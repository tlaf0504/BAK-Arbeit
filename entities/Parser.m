classdef Parser
    
    properties(Constant)
        number_of_triangle_nodes = containers.Map([1,2,3], [3,6,9]);
        number_of_curve_nodes = containers.Map([1,2,3], [2,3,4]);
    end
    
    methods(Static)
        function success = parse(problem_location)
            
            success = 1;
            tmp = pwd;
            cd(problem_location);
            
            % Plot section separator to console
            Misc.print_message(sprintf('%s\n%s\n%s', ...
                Misc.console_section_separator, ...
                'Parser', Misc.console_section_separator));
            
            
            if ~Misc.check_file_existence(problem_location, 'problem_setup.mat')
                msg = sprinft(['Error. Setup-file "problem setup.mat" in folder ', ...
                    '%s not found.'], problem_location);
                Misc.print_error_message(msg);
            end
            load('problem_setup.mat', 'problem_setup');
            
            % Check if the meshing procedure is completed.
            % problem_setup.state must be higher than
            % Setup.setup_state_meshing_finished
            if problem_setup.state < Setup.setup_state_meshing_finished
                msg = sprintf(['Error. Meshing procedure was not ', ...
                    'finished. Maybe the setup procedure was interrupted.', ...
                    'Please rerun the complete meshing procedure.']);
                Misc.print_error_message(msg);
                success = 0;
                return
            end
            
            % ===== Parse mesh file
            tic
            
            Misc.print_message('Parsing mesh file ...')
            [node_data, group_data, triangle_data, curve_data] = Parser. ...
                parse_mesh_file(problem_setup);
            
            time = toc;
            Misc.print_message(sprintf('Done\nElapsed time is %d seconds.\n', time));
            
            
            % ===== Parse settings file
            Misc.print_message('Parsing settings file ...')
            [dirichlet_boundary_conditions, neumann_boundary_conditions, ...
                region_data] = Parser.parse_settings_file(problem_setup);

            Misc.print_message('Done.');
            
            
            % ===== Data prostprocessing
            Misc.print_message('Postprocessing of parsed data....')
            
            dirichlet_boundary_data = ...
                Parser.get_nodes_on_boundary(...
                dirichlet_boundary_conditions, group_data, curve_data, node_data);
            
            neumann_boundary_data = Parser.get_triangle_edges_on_neumann_boundary(...
                neumann_boundary_conditions, ...
                group_data, curve_data, node_data, triangle_data);
            
            material_and_source_properties = ...
                Parser.get_material_and_source_properties_of_elements(region_data, ...
                triangle_data, group_data);
            
            Misc.print_message('Done.');
            
            % ===== Save data to file
            
            data_file = fullfile('internal', 'mesh_data.mat');
            Misc.print_message(sprintf('Writing data to "%s"...', ...
                fullfile(problem_setup.problem_location, data_file)));
            
            
            mesh_data = struct();
            mesh_data.node_data = node_data;
            mesh_data.triangle_data = triangle_data;
            mesh_data.dirichlet_boundary_data = dirichlet_boundary_data;
            mesh_data.neumann_boundary_data = neumann_boundary_data;
            mesh_data.material_and_source_properties = material_and_source_properties;
            
            save(data_file, 'mesh_data');
            Misc.print_message('Done.')
            
            Parser.print_statistics(node_data, triangle_data);
 
            Misc.print_message(Misc.console_section_separator);
            Misc.print_message('\n');
            
            problem_setup.state = Setup.setup_state_parsing_finished;
            Setup.update_problem_setup_file(problem_setup);
            
            cd(tmp);
        end
        
        
        
        
        
        function [node_data, group_data, triangle_data, curve_data] = ...
                parse_mesh_file(problem_setup)
            
            if ~Misc.check_file_existence(problem_setup.problem_location, ...
                    problem_setup.mesh_file)
                error(['Mesh file "%s" does not exist.'], ...
                    fullfile(problem_setup.problem_location, problem_setup.mesh_file))
            end
            
            content = fileread(problem_setup.mesh_file);
            
            [~, node_data] = Parser.parse_nodes(content);
            [~, element_data] = Parser.parse_elements(content);
            [~, group_data] = Parser.parse_physical_groups(content);
            
            clear content;
            
            [~, triangle_data] = Parser. ...
                extract_triangle_elements(...
                Parser.number_of_triangle_nodes(problem_setup.mesh_order), ...
                element_data);
            
            [~, curve_data] = Parser.extract_curve_elements( ...
                Parser.number_of_curve_nodes(problem_setup.mesh_order), ...
                element_data);
        end
        
        
        function [dirichlet_boundary_conditions, neumann_boundary_conditions, ...
                region_data] = ...
                parse_settings_file(problem_setup)
            
            if ~Misc.check_file_existence(problem_setup.problem_location, ...
                    problem_setup.settings_file)
                error(['Settings file "%s" does not exist.'], ...
                    fullfile(problem_setup.problem_location, problem_setup.settings_file))
            end
            
            content = fileread(problem_setup.settings_file);
            
            % ===== Read and parse dirichlet and boundary conditions
            [~, boundary_condition_data] = ...
                Parser.parse_boundary_contidions(content);
            
            % Extract dirichlet boundary conditions
            %
            % The cell elements in the first column contain the name of the boundary
            % condition, the elements of the fourth column the actual values of the
            % conditions.
            idx = strcmpi(boundary_condition_data{1, 3}, 'DRB');
            dirichlet_boundary_conditions = ...
                {boundary_condition_data{1,1}(idx), boundary_condition_data{1,4}(idx)};
            
            % Extract neumann boundary conditions
            idx = strcmpi(boundary_condition_data{1, 3}, 'NRB');
            neumann_boundary_conditions = ...
                {boundary_condition_data{1,1}(idx), boundary_condition_data{1,4}(idx)};
            
            % ===== Read and parse regions
            % A region is a physical surface containing the same material
            % (roh/epsilon/ etc.) and charges/current density etc.
            [~, region_data] = Parser.parse_regions(content);
            
        end
        
        
        function [number_of_boundary_conditions, boundary_condition_data] = ...
                parse_boundary_contidions(file_content)
            
            file_content = Parser.extract_file_section(file_content, '\$Edges', ...
                '\$EndEdges');
            
            % Remove first '$Edges' from section
            file_content = replace(file_content, '$Edges', '');
            
            % First line after '$Edges' contains the number of boundary conditions
            [number_of_boundary_conditions, ~, ~, idx] = sscanf(file_content, ...
                '%d', 1);
            
            % Each further line contains one boundary condition of the following form:
            %
            % <name_of_boundary_condition:string> <Type ID:int> ...
            % <Type of boundary condition:string> <Condition value:float>
            %
            % The Type ID can have two different values: 1 for edge boundary
            % conditions and 2 for face boundary conditions
            %
            % Each boundary condition ca be of the following type (specified by the
            % third entry in each line):
            %
            % 'DRB': Dirichlet boundary condition
            % 'NRB': Neumann boundary condition
            % 'source': Source condition
            boundary_condition_data = textscan(file_content(idx+1:end), ...
                '%s %d %s %f', number_of_boundary_conditions);
        end
        
        function [number_of_regions, region_data] = ...
                parse_regions( file_content)
            
            file_content = Parser.extract_file_section(file_content, '\$Faces', ...
                '\$EndFaces');
            file_content = replace(file_content, '$Faces', '');
            
            [number_of_regions, ~, ~, idx] = sscanf(file_content, ...
                '%d', 1);
            
            % Each region is of the following form: (Interpretation of last two
            % parameters depend on the problem type)
            %
            % <Name:string> <ID:int> <mu_r or epsilon_r:float> <roh or J:float>
            region_data = textscan(file_content(idx+1:end), ...
                '%s %d %f %f', number_of_regions);
        end
        
        function [number_of_nodes, node_data] = parse_nodes(file_content)
            
            file_content = Parser.extract_file_section(file_content, '\$Nodes', '\$EndNodes');
            file_content = replace(file_content, '$Nodes', '');
            
            [number_of_nodes, ~, ~, idx] = sscanf(file_content, '%d', 1);
            
            node_data = textscan(file_content(idx+1:end), '%d %f %f %f', ...
                number_of_nodes);
            
            % First element contains vector of node indices, second column contains
            % matrix of node coordinates.
            node_data = {node_data{1}, [node_data{2}, node_data{3}, node_data{4}]};
        end
        
        function [number_of_elements, element_data] = parse_elements(file_content)
            
            
            file_content = Parser.extract_file_section(file_content, '\$Elements', ...
                '\$EndElements');
            file_content = replace(file_content, '$Elements', '');
            
            [number_of_elements, ~, ~, idx] = sscanf(file_content, '%d', 1);
            
            element_data = textscan(file_content(idx:end), repmat('%d ', 1, 30), ...
                number_of_elements);
            
            element_data = cell2mat(element_data);
        end
        
        function [number_of_groups, group_data] = parse_physical_groups(file_content)
            
            file_content = Parser.extract_file_section(file_content, '\$PhysicalNames', ...
                '\$EndPhysicalNames');
            file_content = replace(file_content, '$PhysicalNames', '');
            
            [number_of_groups, ~, ~, idx] = sscanf(file_content, '%d', 1);
            
            group_data = textscan(file_content(idx:end), '%d %d %s', ...
                number_of_groups);
            
            % Removes the double quotes before and after the names of the physical
            % groups.
            group_names = group_data{1,3};
            for k = 1 : length(group_names)
                group_names{k} = group_names{k}(2:end-1);
            end
            
            group_data{1,3} = group_names;
        end
        
        function [number_of_triangles, triangle_data] = extract_triangle_elements( ...
                number_of_nodes, element_data)
            
            % Depending on the calling function, number of nodes can either be input
            % as an integer ir double, although it always contains an integer value.
            number_of_nodes = int32(number_of_nodes);
            switch number_of_nodes
                case 3
                    element_tag = 2;
                case 6
                    element_tag = 9;
                case 9
                    element_tag = 20;
            end
            
            triangles = element_data(element_data(:, 2) == element_tag, :);
            [number_of_triangles, ~] = size(triangles);
            number_of_tags = triangles(:, 3);
            
            % Column 1: Element id, Column 2: number of tags --> Nodes start at column
            % 3 + number_of_tags
            node_columns_start = 4 + number_of_tags;
            
            
            
            % For each triangle, get column indices of nodes. The column  indices are
            % calculated seperately, which allows the algorithm to extract the values
            % from the correct columns even if some (or all) rows have a different
            % number of tags.
            %
            % The resulting matrix has the following form:
            %[ col_idx_elem_1_node_1, col_idx_elem_1_node_2, ... ;
            %  col_idx_elem_2_node_1, col_idx_elem_1_node_2, ... ;
            % .
            % .
            % .
            % ]
            col_idx = (repmat(node_columns_start, 1, number_of_nodes) + ...
                repmat(0 : number_of_nodes - 1, number_of_triangles, 1));
            
            % The resulting matrix triangle_data should also contain information
            % about the id of the elements and physical group the element belongs to.
            % The id of the the elements is containes in the first column, the id of
            % the pyhsical group in the fourth column.
            %
            % The matrix col_idx is therefore appended to extract also the element
            % ids and the physical group ids.
            %
            % The resulting matrix has the following form:
            %[ 1, 4, col_idx_elem_1_node_1, col_idx_elem_1_node_2, ... ;
            %  1, 4, col_idx_elem_2_node_1, col_idx_elem_1_node_2, ... ;
            % .
            % .
            % .
            % ]
            %
            col_idx = [ones(number_of_triangles, 1), ...
                4*ones(number_of_triangles, 1), col_idx];
            
            % Now the matrix is transformed to a vector to be used in the Matlab
            % buit-in function sub2ind() (See linear matrix indexing of matlab for
            % further information)
            %
            % The resulting (column) vector has the following form:
            %
            % [1; 4; col_idx_elem_1_node_1; ... col_idx_elem_1_node_n; 1; 4; ...
            % col_idx_elem_2_node_1, ... ];
            col_idx = col_idx'; col_idx = col_idx(:);
            
            
            number_of_columns_per_element = number_of_nodes + 2;
            
            % To generate the linear indices of the matrix elements to extract, two
            % vectors are needed. One specifying the columns of the elements (see
            % col_idx) and one specifying the rows. (row_idx)
            %
            % As we want to extract
            % number_of_columns_per_element := number_of_nodes + 2
            % elements from the matrix (id of element, id of physical group, nodes)
            % the vector specifying the rows contains number_of_columns_per_element *
            % number_of_triangles entries. The vector is of the following form:
            %
            %   (number_of_columns_per_element 1s),
            % [ 1, 1, ............................., 1,
            %
            %  (number_of_columns_per_element 2s)
            %  2, 2, .............................., 2,
            %
            %  (number_of_columns_per_element 3s)
            %  3, ........]
            row_idx = repmat(1:number_of_triangles, number_of_columns_per_element, 1);
            row_idx = row_idx(:);
            
            triangle_data_idx = sub2ind(size(triangles), row_idx, col_idx);
            triangle_data = triangles(triangle_data_idx);
            triangle_data = reshape(triangle_data, number_of_columns_per_element, [])';
        end
        
        function [number_of_curves, curve_data] = extract_curve_elements( ...
                number_of_nodes, element_data)
            
            % function [number_of_curves, curve_data] =  extract_curve_elements( ...
            %             number_of_nodes, element_data)
            %
            % DESCRIPTION:
            % The functionality of this function matches withc the functionality of
            % extract_triangle_elements(). Take a look at this funciton for further
            % description of the algorithm.
            %
            %
            
            
            number_of_nodes = int32(number_of_nodes);
            switch number_of_nodes
                case 2
                    element_tag = 1;
                case 3
                    element_tag = 8;
                case 4
                    element_tag = 26;
            end
            
            curves = element_data(element_data(:, 2) == element_tag, :);
            [number_of_curves, ~] = size(curves);
            number_of_tags = curves(:, 3);
            node_columns_start = 4 + number_of_tags;
            
            col_idx = (repmat(node_columns_start, 1, number_of_nodes) + ...
                repmat(0 : number_of_nodes - 1, number_of_curves, 1));
            
            col_idx = [ones(number_of_curves, 1), ...
                4*ones(number_of_curves, 1), col_idx];
            
            col_idx = col_idx'; col_idx = col_idx(:);
            
            
            number_of_columns_per_element = number_of_nodes + 2;
            
            row_idx = repmat(1:number_of_curves, number_of_columns_per_element, 1);
            row_idx = row_idx(:);
            
            curve_data_idx = sub2ind(size(curves), row_idx, col_idx);
            curve_data = curves(curve_data_idx);
            curve_data = reshape(curve_data, number_of_columns_per_element, [])';
            
            
        end
        
        function section = extract_file_section(file_content, section_start_marker, ...
                section_end_marker)
            regex = [section_start_marker, '(.*)', section_end_marker];
            [start_idx, end_idx] = regexp(file_content, regex, 'once');
            
            if isempty(start_idx) || isempty(end_idx)
                Misc.print_error_message(sprintf(['Error: could not extract file section ' ...
                    'between markers "%s" and "%s". Check file and markers.'], ...
                    section_start_marker, section_end_marker));
            end
            
            section = file_content(start_idx : end_idx);
        end
        
        function node_values_on_boundary = get_nodes_on_boundary(...
                boundary_conditions, group_data, curve_data, node_data)
            
            number_of_nodes = length(node_data{1});
            number_of_boundary_conditions = length(boundary_conditions{1});
            
            node_values_on_boundary = -ones(number_of_nodes, 1);
            
            for k = 1 : number_of_boundary_conditions
                name_of_current_boundary_condition = ...
                    boundary_conditions{1,1}{k};
                
                value_of_current_boundary_condition = ...
                    boundary_conditions{1,2}(k);
                
                % Get the index of the physical group corresponding to the current
                % boundary condition within the group_data cell array
                idx_of_boundary_condition_group = find(...
                    strcmpi(name_of_current_boundary_condition, group_data{1,3}));
                
                % Get the id of the physical group of the current boundary condition
                id_of_boundary_condition_group = ...
                    group_data{1,2}(idx_of_boundary_condition_group);
                
                % For the 2D cases handled in this software, the boundary
                % conditions are always curves.
                curves_of_current_boundary_condition = ...
                    curve_data(:, 2) == id_of_boundary_condition_group;
                
                
                nodes_of_boundary_curves = ...
                    curve_data(curves_of_current_boundary_condition, 3:end);
                
                nodes_of_boundary_curves = unique(nodes_of_boundary_curves(:));
                
                node_values_on_boundary(nodes_of_boundary_curves) = ...
                    value_of_current_boundary_condition;
                
                msg = sprintf('Nodes on boundary "%s": \n\t%s', ...
                    name_of_current_boundary_condition, ...
                    sprintf('%d ', nodes_of_boundary_curves));
                Misc.print_message(msg);
            end
            
            ids_of_nodes_on_boundary = int32(find(node_values_on_boundary ~= -1));
            node_values_on_boundary = {ids_of_nodes_on_boundary, ...
                node_values_on_boundary(ids_of_nodes_on_boundary)};
            
        end
        
        function neumann_boundary_data = get_triangle_edges_on_neumann_boundary(...
                neumann_boundary_conditions, ...
                group_data, curve_data, node_data, triangle_data)
            
            
            neumann_boundary_data = cell(1,3);
            
            number_of_boundary_conditions = length(neumann_boundary_conditions{1});
            for k = 1 : number_of_boundary_conditions
                
                current_boundary_condition = ...
                    {neumann_boundary_conditions{1}(k), neumann_boundary_conditions{2}(k)};
                
                
                tmp = Parser.get_nodes_on_boundary( ...
                    current_boundary_condition, group_data, curve_data, node_data);
                nodes_of_boundary = tmp{1};
                
                values = tmp{2};
                clear tmp;
                
                % Get matrix of element nodes
                nodes_of_elements = triangle_data(:, 3:end);
                
                % Matrix of binary values, containing a logic 1 at each element on
                % the current neumann boundary;
                % Matrix has a dimension of number_of_triangles x nodes_per_element
                %
                % Example:
                %
                % 4 elements, where 2 of them are located on the neumann boundary
                % [ 0, 0, 0;
                %   1, 0, 1;
                %   0, 0, 0;
                %   0, 1, 1]
                %
                element_nodes_on_boundary = ...
                    ismember(nodes_of_elements, nodes_of_boundary);
                
                % Vector of binary values containing a logic 1 at each element whose
                % corresponding triangle belongs to the neumann boundary
                %
                % Example for elements from example above:
                % [0; 1; 0; 1]
                elements_on_boundary = sum(element_nodes_on_boundary,2) > 0;
                
                % Remove elements from the matrix which are not on the boundary
                element_nodes_on_boundary = ...
                    element_nodes_on_boundary(elements_on_boundary, :);
                
                edges_on_boundary = ...
                    Parser.assign_nodes_on_neumann_boundary_to_triangle_edges(...
                    triangle_data(elements_on_boundary, 1), element_nodes_on_boundary);
                
                [N, ~] = size(edges_on_boundary);
                
                % Append ids of elements
                neumann_boundary_data{1} = [neumann_boundary_data{1}; ...
                    edges_on_boundary(:, 1)];
                
                % Append element sides
                neumann_boundary_data{2} = [neumann_boundary_data{2}; ...
                    edges_on_boundary(:, 2)];
                
                % Append boundary condition values
                neumann_boundary_data{3} = [neumann_boundary_data{3}; ...
                    values(1) * ones(N, 1)];
            end
        end
        
        function edges_on_boundary = ...
                assign_nodes_on_neumann_boundary_to_triangle_edges(...
                boundary_element_ids, element_nodes_on_boundary)
            
            edges_on_boundary = [];
            
            [~, number_of_nodes] = size(element_nodes_on_boundary);
            
            if number_of_nodes == 3
                nodes_of_triangle_edges = ...
                    FirstOrderTriangleElement.nodes_of_triangle_edges;
            elseif number_of_nodes == 6
                nodes_of_triangle_edges = ...
                    SecondOrderTriangleElement.nodes_of_triangle_edges;
            elseif number_of_nodes == 9
                nodes_of_triangle_edges = ...
                    ThirdOrderTriangleElemement.nodes_of_triangle_edges;
            end
            
            
            for k = 1 : 3
                nodes_of_current_edge = nodes_of_triangle_edges(k, :);
                
                idx = sum(element_nodes_on_boundary == nodes_of_current_edge, 2) == ...
                    number_of_nodes;
                
                ids = boundary_element_ids(idx);
                edge = k;
                
                N = sum(idx);
                
                edges_on_boundary = [edges_on_boundary; [ids(:), edge * ones(N, 1)]];
            end
        end
        
        function material_and_source_properties = ...
                get_material_and_source_properties_of_elements(region_data, ...
                triangle_data, group_data)
            
            [number_of_elements, ~] = size(triangle_data);
            
            % First column: material properties, second column: source value
            material_and_source_properties = zeros(number_of_elements, 2);
            
            % Names of physical groups
            group_names = group_data{3};
            
            % Names of regions. A 'region' is a physical group with constant material
            % properties and a constant source value.
            region_names = region_data{1};
            
            idx = ismember(group_names, region_names);
            group_ids_of_regions = group_data{2}(idx);
            
            
            for k = 1 : length(region_names)
                current_group_id = group_ids_of_regions(k);
                current_material_property = region_data{3}(k);
                current_source_value = region_data{4}(k);
                
                ids_of_elements = triangle_data(:, 2) == current_group_id;
                
                % Number of elements (triangles) belonging to the current group
                N = sum(ids_of_elements);
                
                
                material_and_source_properties(ids_of_elements, :) = ...
                    repmat([current_material_property, current_source_value], ...
                    N, 1);
            end
        end
        
        function print_statistics(node_data, triangle_data)
            
            [number_of_nodes, ~] = size(node_data{2});
            [number_of_elements, ~] = size(triangle_data);
            
            msg = sprintf(...
            ['\nStatistics:\n', ...
                Misc.console_subsection_separator, '\n', ...
                'Number of elements: %d\n', ...
                'Number of nodes: %d\n', ...
                Misc.console_subsection_separator, '\n'], ...
                number_of_elements, number_of_nodes);
            Misc.print_message(msg);
        end
        
    end
end
