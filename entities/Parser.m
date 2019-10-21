classdef Parser
    
    properties(Constant)
        number_of_triangle_nodes = containers.Map(int32([1,2,3]), ...
            int32([3,6,10]));
        number_of_curve_nodes = containers.Map(int32([1,2,3]), ...
            int32([2,3,4]));
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
                msg = sprintf(['Error. Setup-file "problem setup.mat" in folder ', ...
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
            problem_setup.mesh_data_file = 'mesh_data.mat';
            problem_setup.mesh_data_file_hash = DataHash(data_file);
            
            Setup.update_problem_setup_file(problem_setup);
            
            cd(tmp);
        end
        
        
        
        
        
        function [node_data, group_data, triangle_data, curve_data] = ...
                parse_mesh_file(problem_setup)
            
        % function [node_data, group_data, triangle_data, curve_data] = ...
        %        parse_mesh_file(problem_setup)
        %
        % DESCRIPTION: Function for parsing data from a Gmsh version 2.2 mesh file.
        %
        % INPUT:
        %     problem_setup.....Structure containing all necessary information about
        %         the problem.
        %
        % OUTPUT:
        %     node_data.....Structure containing all information about the nodes. 
        %         Returned by parse_nodes() .
        %     group_data.....Structure containing all information about the physical
        %         groups. Returned by parse_physical_groups().
        %     triangle_data.....Structure containing all information about the
        %         FEM triangles. Returned by extract_triangle_elements()
        %     curve_data.....Structure containing all information about the FEM
        %         curves. Returned by extract_curve_elements()
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
        
            % Check if the mesh file exists.
            if ~Misc.check_file_existence(problem_setup.problem_location, ...
                    problem_setup.mesh_file)
                error(['Mesh file "%s" does not exist.'], ...
                    fullfile(problem_setup.problem_location, problem_setup.mesh_file))
            end
            
            % Get content of mesh file as character array.
            content = fileread(problem_setup.mesh_file);
            
            % Parse nodes, elements and physical groups
            node_data = Parser.parse_nodes(content);
            element_data = Parser.parse_elements(content);
            group_data = Parser.parse_physical_groups(content);
            
            clear content;
            
            % Extract triangle elements from element_data
            triangle_data = Parser. ...
                extract_triangle_elements(...
                Parser.number_of_triangle_nodes(problem_setup.mesh_order), ...
                element_data);
            
            curve_data = Parser.extract_curve_elements( ...
                Parser.number_of_curve_nodes(problem_setup.mesh_order), ...
                element_data);
        end
        
        
        function [dirichlet_boundary_conditions, neumann_boundary_conditions, ...
                region_data] = ...
                parse_settings_file(problem_setup)
            
        % function [dirichlet_boundary_conditions, neumann_boundary_conditions, ...
        %        region_data] = ...
        %        parse_settings_file(problem_setup)
        %
        % DESCRIPTION: Function for parsing data from the settings file. 
        %
        % INPUT:
        %     problem_setup.....Structure containing all necessary information about
        %         the problem.
        %
        % OUTPUT:
        %     dirichlet_boundary_conditions.....Structure containing all information
        %         about the dirichlet boundary conditions. It contains the following
        %         members:
        %         -) number_of_boundary_conditions
        %         -) names.....Names of the physical curves representing the boundary
        %             conditions. Each of the names must correspond to a name of a
        %             physical curve in the mesh file.
        %         -) values.....Boundary condition values
        %     neumann_boundary_conditions.....Same structure as 
        %         dirichlet_boundary_conditions, but for neumann boundary conditions.
        %     region_data.....Structure containing all information about the problem
        %         regions. Each region is represented by a physical surface in the
        %         mesh file. Returned by parse_regions().
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            % Check if settings file exists
            if ~Misc.check_file_existence(problem_setup.problem_location, ...
                    problem_setup.settings_file)
                error(['Settings file "%s" does not exist.'], ...
                    fullfile(problem_setup.problem_location, ...
                    problem_setup.settings_file))
            end
            
            % Get complete settings file content as character array
            content = fileread(problem_setup.settings_file);
            
            % ===== Read and parse dirichlet and neumann boundary conditions
            boundary_condition_data = ...
                Parser.parse_boundary_contidions(content);
            
            % Extract dirichlet boundary conditions by searching for boundary
            % conditions of type 'DRB'
            idx = strcmpi(boundary_condition_data.types, 'DRB');
            
            dirichlet_boundary_conditions = struct( ...
                'number_of_boundary_conditions', sum(idx), ...
                'values', boundary_condition_data.values(idx)...
                );
            
            % Cell-arrays containing strings must be added in this way, otherwise
            % matlab creates multidimensional structures. 
            % (Matlab R2018a, Linux version)
            dirichlet_boundary_conditions.names = boundary_condition_data.names(idx);
           
            
            % Extract neumann boundary conditions
            idx = strcmpi(boundary_condition_data.types, 'NRB');
            neumann_boundary_conditions = struct( ...
                'number_of_boundary_conditions', sum(idx), ...
                'values', boundary_condition_data.values(idx));
            
            neumann_boundary_conditions.names =  boundary_condition_data.names(idx);
            
            % ===== Read and parse regions
            % A region is a physical surface containing the same material
            % (roh/epsilon/ etc.) and charges/current density etc.
            region_data = Parser.parse_regions(content);
            
        end
        
        
        function boundary_condition_data = parse_boundary_contidions(file_content)
            
        % function boundary_condition_data = parse_boundary_contidions(file_content)
        %
        % DESCRIPTION: Function for parsing the boundary conditions from the settings
        %     file.
        %
        %     As this program only supports 2D problems, the boundary conditions are
        %     always one-dimensional.
        %     The boundary conditions are contained in the '$Edges'-section of the
        %     settings file.
        %
        %     Each name of a boundary condition must corredpond to the name of a
        %     physical curve in the mesh file.
        %
        % INPUT:
        %     file_content.....Character array containing the content of the settings
        %         file.
        %
        % OUTPUT:
        %     boundary_condition_data.....Structure containing all information about 
        %         the parsed boundary conditions.
        %         It contains the following members:
        %         -) number_of_boundary_conditions.....Number of parsed boundary
        %             conditions
        %         -) names.....Cell array containing the names of the boundary
        %             conditions. This names must match one of the physical curves in
        %             the mesh file.
        %         -) types.....Cell array containing the types of the boundary
        %             conditions. Can be either 'NRB' for neumann boundary conditions
        %             or 'DRB' for dirichlet boundary conditions.
        %         -) values.....Vector of double values containing the actual values
        %             of the boundary conditions.
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            % Extract file content between '%Edges' and '$EndEdges'
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
            % The Type ID is always 1 for entries in the '$Edges' - section of the
            % file.
            %
            % Each boundary condition can be of one of the following types 
            % (specified by the third entry in each line):
            %
            % 'DRB': Dirichlet boundary condition
            % 'NRB': Neumann boundary condition
            boundary_condition_data_tmp = textscan(file_content(idx+1:end), ...
                '%s %d %s %f', number_of_boundary_conditions);
            
            % Store data to output structure
            boundary_condition_data = struct(...
                'number_of_boundary_conditions', number_of_boundary_conditions, ...
                'values', boundary_condition_data_tmp{4});
            
            % Cell-arrays containing strings must be added in this way, otherwise
            % matlab creates multidimensional structures. 
            % (Matlab R2018a, Linux version)
            boundary_condition_data.names = boundary_condition_data_tmp{1};
            boundary_condition_data.types = boundary_condition_data_tmp{3};
            
        end
 
        function region_data = parse_regions(file_content)
            
        % function region_data = parse_regions(file_content)
        %
        % DESCRIPTION: Function for parsing the regions defined in the '$Faces'
        %     section of the settings file.
        %
        %     As this program only supports 2D problems, regions are always
        %     two-dimensional. Their definition is contained in the '$Faces' section
        %     of the settings file.
        %    
        %     Each region name must correspond to the name of a physical surface in 
        %     the mesh file.
        %
        % INPUT:
        %     file_content.....Character array containing the content of the settings
        %         file.
        %
        % OUTPUT:
        %     region_data.....Structure containing all information about the parsed 
        %         boundary regions.
        %         It contains the following members:
        %         -) number_of_regions.....Number of parsed regions
        %         -) names.....Cell array containing the names of the regions.
        %             This names must match one of the physical surfaces in the mesh 
        %             file.
        %         -) material_values.....Vector of float-values containing the
        %             material property. (e.g. the relative permettivity)
        %             This program only supports isotropic materials. 
        %         -) source_values.....Source values of each region. (e.g. current
        %             density or free colume charges)
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            file_content = Parser.extract_file_section(file_content, '\$Faces', ...
                '\$EndFaces');
            file_content = replace(file_content, '$Faces', '');
            
            [number_of_regions, ~, ~, idx] = sscanf(file_content, ...
                '%d', 1);
            
            % Each region is of the following form: (Interpretation of last two
            % parameters depend on the problem type)
            %
            % <Name:string> <Type ID:int> <mu_r or epsilon_r:float> <roh or J:float>
            %
            % Type ID specifies the dimension of the region. In the '$Faces' section,
            % this parameter is always set to one. (As regions are 2 dimensional ;P )
            region_data_tmp = textscan(file_content(idx+1:end), ...
                '%s %d %f %f', number_of_regions);
            
            region_data = struct( ...
                'number_of_regions', number_of_regions, ...
                'material_values', region_data_tmp{3}, ...
                'source_values', region_data_tmp{4}...
                );
            
            % Cell-arrays containing strings must be added in this way, otherwise
            % matlab creates multidimensional structures. 
            % (Matlab R2018a, Linux version)
            region_data.names = region_data_tmp{1};
        end
        
        function node_data = parse_nodes(file_content)
        
        % function [number_of_nodes, node_data] = parse_nodes(file_content)
        %
        % DESCRIPTION: Function for parsing the nodes from a Gmsh version 2.2 mesh
        %     file.
        %
        % INPUT:
        %     file_content.....Character array containing the content of the mesh
        %         file.
        %
        % OUTPUT:
        %     node_data.....Structure containing all information about the parsed
        %         nodes. It contains the following members:
        %         -) number_of_nodes.....Number of parsed nodes
        %         -) ids.....IDs of the nodes
        %         -) coordinates.....number_of_nodes x 3 matrix containing the
        %             coordinates of the nodes. The first column contains the x-, the
        %             second column the y- and the third column the z-coordinates of
        %             the node. 
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            
            % Extract the file content between '$Nodes' and '$EndNodes' from the mesh
            % file.
            file_content = Parser.extract_file_section(file_content, '\$Nodes', ...
                '\$EndNodes');
            
            % Remove '$Nodes' from the content
            file_content = replace(file_content, '$Nodes', '');
            
            % After '$Nodes', the file contains the number of nodes
            [number_of_nodes, ~, ~, idx] = sscanf(file_content, '%d', 1);
            
            % The next lines contain each an integer number (representing the ID of
            % the node) and three floating point numbers representing the coordinates
            % of the node. First value: x-coordinate, second value: y-coordinate and
            % third value: z-coordinate.
            node_data_tmp = textscan(file_content(idx+1:end), '%d %f %f %f', ...
                number_of_nodes);
           
            node_data = struct(...
                'number_of_nodes', number_of_nodes, ...
                'ids', node_data_tmp{1}, ...
                'coordinates', [node_data_tmp{2}, node_data_tmp{3}, node_data_tmp{4}]);
            
        end
        
        function element_data = parse_elements(file_content)
            
        % function [number_of_elements, element_data] = parse_elements(file_content)
        %
        % DESCRIPTION: Function for parsing the elements from a Gmsh version 2.2 mesh
        %     file.
        %
        % INPUT:
        %     file_content.....Character array containing the content of the mesh
        %         file.
        %
        % OUTPUT:
        %     element_data.....Matrix containing all information about the parsed
        %         elements. The structure of the matrix is described in the code
        %         below. Use e.g. extract_triangle_elements() or
        %         extract_curve_elements() to extract the corresponding elements from
        %         the matrix.
        % 
        %
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            % Extract file content between '$Elements' and '$EndElements' markers.
            file_content = Parser.extract_file_section(file_content, '\$Elements', ...
                '\$EndElements');
            
            % Remove first line containing '$Elements'
            file_content = replace(file_content, '$Elements', '');
            
            % Next line contains number of elements
            [number_of_elements, ~, ~, idx] = sscanf(file_content, '%d', 1);
            
            
            % Parsing the data of all elements is a little bit more complicated as
            % each element definition can consist of a different number of values.
            % e.g. a 2D line is defined by less entries than a 3D triangle. ( See
            % Gmsh manual for more information)
            %
            % Whats more is that the number of tags can vary from element to element.
            % 
            % The algorithm always tries to parse 30 integer values, resulting in 
            % a 1x30 cell array as a result for element_data_tmp, where each cell
            % contains a vector of integer numbers with number_of_elements entries.
            %
            % The first cell contains the elements number.
            % The second cell contains the elements type.
            % The third cell contains the number of the 'tags'. 
            % The next number_of_tags cells contain the elements tags.
            % The next number_of_elements_nodes cells contain the ids of the elements
            % nodes.
            % The remaining cells contain zeros.
            element_data = textscan(file_content(idx:end), repmat('%d ', 1, 30), ...
                number_of_elements);

            % As element_data contains all integer values, it is transformed to a
            % matrix.
            element_data = cell2mat(element_data);
        end
        
        function group_data = parse_physical_groups(file_content)
            
        % function group_data = parse_physical_groups(file_content)
        %
        % DESCRIPTION: Function for parsing the physical groups from a version 2.2
        %     mesh file.
        %
        % INPUT:
        %     file_content.....Character array containing the content of the mesh
        %         file.
        %
        % OUTPUT:
        %     group_data.....Structure containing the members described below, 
        %         holding all information about the physical groups contained in the 
        %         problem.
        %         Member description:
        %
        %         -) number_of_groups.....Number of physical gorups
        %         -) dimensions.....Vector of number_of_groups integer numbers
        %             containing the dimension of the physical group. e.g. 1 for a
        %             curve or 2 for a 2D surface. 
        %         -) tags.....Tags of the physical groups. Each element in the
        %             'Elements'-section of the mesh file contains beside its ID,
        %             type and nodes also some 'tags' specifying e.g. the membership 
        %             of a physical group of the current element.
        %             The integer values contained in this vector are therefore used
        %             to assign a element to a physical group.
        %        -) names.....Cell array containing the names of the physical groups.
        % 
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
        
            file_content = Parser.extract_file_section(file_content, '\$PhysicalNames', ...
                '\$EndPhysicalNames');
            
            % Remove first line containing '$PhysicalNames'
            file_content = replace(file_content, '$PhysicalNames', '');
            
            % Next line contains the number of physical groups
            [number_of_groups, ~, ~, idx] = sscanf(file_content, '%d', 1);
            
            % Every group definition contains of two integer values (dimension and
            % tag) and a string value (name of the physical group)
            %
            % Group data_tmp is a 1x3 cell array containg in the first and second
            % cell an integer vector with the dimension and tag values, and in the
            % third cell the names of the physical groups.
            group_data_tmp = textscan(file_content(idx:end), '%d %d %s', ...
                number_of_groups);
            
            % Removes the double quotes before and after the names of the physical
            % groups.
            group_names = group_data_tmp{1,3};
            for k = 1 : length(group_names)
                group_names{k} = group_names{k}(2:end-1);
            end
            
            
            group_data = struct(...
                'number_of_groups', number_of_groups, ...
                'dimensions', group_data_tmp{1,1}, ...
                'tags', group_data_tmp{1,2});
            
            % Has to be added now, because when adding in the expression above,
            % Matlab creates a 5x1 struct due to cell array of strings.
            % (Behaviour of Matlab R2018a, Linux version)
            group_data.names = group_names;
            
        end
        
        function triangle_data = extract_triangle_elements( ...
                number_of_element_nodes, element_data)
                 
        % function [number_of_triangles, triangle_data] = extract_triangle_elements( ...
        %        number_of_nodes, element_data)
        %
        % DESCRIPTION: Function for extracting the triangle elements from the
        %     element_data matrix. 
        %
        %     The following procedure is used to extract the required data from the
        %     element_data matrix:
        %
        %     1.) Determination of the matrix rows that contain data of triangles by
        %         using the element-type id contained in the second column of each 
        %         row.
        %     2.) Calculation of the indices of columns containing ellement ID,
        %         physical group and nodes. Column indices are stored in variable
        %         'col_idx'.
        %     3.) Extraction of data using linear indexing.
        %
        % INPUT:
        %     number_of_element_nodes.....Number of nodes per element. (Integer value)
        %     element_data.....Matrix of integer values containing the element
        %     definitions. Returned by parse_elements().
        %     See parse_elements()-header for definition of this matrix.
        %
        % OUTPUT:
        %     triangle_data.....Struct containing all information about the FEM
        %         triangles of the problem area
        %         The structure contains the following members:
        %         -) number_of_triangles.....Number of triangles in the problem area.
        %         -) IDs.....Integer vector containing the IDs of the triangles
        %         -) physical_groups.....Integer vector containing the ids of the
        %             physical groups the triangles are assigned to. Each triangle
        %             can only be assigned to one physical group.
        %         -) nodes.....number_of_triangles x number_of_element_nodes matrix
        %             of integer values containing the ids of the global nodes 
        %             corresponding to the elements local nodes. 
        %             For the node ordering, see 
        %             http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
        %             (05.07.2019)
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
        
        
            
        % See http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format for definition
        % of the element types.
        %
        % Element-type IDs of supported triangles:
        %
        % 2: First order (3-node) triangle
        % 9: Second order (6-node) triangle
        % 21: Thrid order (10-node) triangle
            switch number_of_element_nodes
                case 3
                    element_type_ID = 2;
                    finite_element = FirstOrderTriangleElement;
                case 6
                    element_type_ID = 9;
                    finite_element = SecondOrderTriangleElement;
                case 10
                    element_type_ID = 21;
                    finite_element = ThirdOrderTriangleElement;
            end
            
            % Second column of element_data contains element-type IDs.
            % Variable 'triangles' contains complete rows of element_data for all 
            % triangles. 
            triangles = element_data(element_data(:, 2) == element_type_ID, :);
            
            [number_of_triangles, ~] = size(triangles);
            
            % The third column of each row contains the number of tags for each
            % triangle.
            number_of_tags = triangles(:, 3);
            
            % Column 1: Element id, Column 2: number of tags --> Nodes start at column
            % 3 + number_of_tags
            node_columns_start = 4 + number_of_tags;
            
            
            % Get column indices of nodes for each triangle. The column indices of
            % the node ids are calculated seperately because triangles can have
            % different number of tags.
            %
            % The resulting matrix has the following form and dimention:
            %[ col_idx_elem_1_node_1, col_idx_elem_1_node_2, ... ;
            %  col_idx_elem_2_node_1, col_idx_elem_1_node_2, ... ;
            % .
            % .
            % .
            % ]
            % Dimension: number_of_triangles x number_of_nodes
            %
            col_idx = (repmat(node_columns_start, 1, number_of_element_nodes) + ...
                repmat(0 : number_of_element_nodes - 1, number_of_triangles, 1));
            
            % The resulting structure 'triangle_data'should also contain information
            % about the id of the elements and physical group the elements belong to.
            % The ids of the the elements are contained in the first column, the ids 
            % of the pyhsical groups in the fourth column.
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
            
            % The matrix is now transformed to a vector to allow linear index
            % calculation by the built-in function sub2ind. (See corresponding
            % manpage for further information)
            %
            % The resulting (column) vector has the following form:
            %
            % [1; 4; col_idx_elem_1_node_1; ... ; col_idx_elem_1_node_n; 1; 4; ...
            % col_idx_elem_2_node_1, ... ];
            col_idx = col_idx'; col_idx = col_idx(:);
            
            % Calculate how many values should be extracted from each row.
            number_of_columns_per_element = number_of_element_nodes + 2;
            
            % To generate the linear indices of the matrix elements to extract, two
            % vectors are needed. One specifying the columns of the elements (see
            % col_idx) and one specifying the rows. (row_idx)
            %
            % The vector specifying the rows is calculated as following:
            %
            % We want to extract
            % number_of_columns_per_element := number_of_nodes + 2
            % elements from the matrix (id of element, id of physical group, nodes).
            %
            % Therefore the vector specifying the rows contains
            %
            %     number_of_columns_per_element * number_of_triangles 
            %
            % entries. 
            % The resulting vector is of the following form:
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
            
            % Calculation of linear matrix indices of the required elements
            triangle_data_idx = sub2ind(size(triangles), row_idx, col_idx);
            
            
            % Linear indexing delivers vector results -> triangle_data_tmp must be
            % reshaped. 
            triangle_data_tmp = triangles(triangle_data_idx);
            
            % Reshaping vector to number_of_triangles x number_of_columns_per_element
            % matrix.
            % The first column contains the element IDs, the second column the
            % physical group tags the the remaining columns the element nodes.
            triangle_data_tmp = reshape(triangle_data_tmp, ...
                number_of_columns_per_element, [])';
            
            % For second and third order elements, the node ordering must be
            % corrected as Gmsh uses a different ordering than this program.
            % See http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering for further
            % information
            %
            % The node orderings used in this program are defined in the
            % corresponding classes as comment headers.
            %
            if number_of_element_nodes > 3
                node_ordering = finite_element.triangle_node_ordering;
                triangle_data_tmp(:, 2 + node_ordering(:, 2)) = ...
                    triangle_data_tmp(:, 2 + node_ordering(:, 1));
            end
            
            
            % Store data to output structure
            triangle_data = struct( ...
                'number_of_triangles', number_of_triangles, ...
                'IDs', triangle_data_tmp(:, 1), ...
                'physical_groups', triangle_data_tmp(:, 2), ...
                'nodes', triangle_data_tmp(:, 3:end));
        end
        
        function curve_data = extract_curve_elements( ...
                number_of_element_nodes, element_data)
            
        % function [number_of_curves, curve_data] =  extract_curve_elements( ...
        %             number_of_nodes, element_data)
        %
        % DESCRIPTION: Function for extracting the curve elements from element_data.
        %
        % The calculations of this function matches in large parts with the
        % calculations of extract_curve_elements(). 
        % See extract_curve_elements() for deeper information about the algorithm.
        %
        %
        % INPUT:
        %     number_of_element_nodes.....Number of nodes per element. (Integer value)
        %     element_data.....Matrix of integer values containing the element
        %     definitions. Returned by parse_elements().
        %     See parse_elements()-header for definition of this matrix.
        %
        % OUTPUT:
        %     curve_data.....Struct containing all information about the FEM
        %         curves of the problem area
        %         The structure contains the following members:
        %         -) number_of_curves.....Number of curves in the problem area.
        %         -) IDs.....Integer vector containing the IDs of the curves
        %         -) physical_groups.....Integer vector containing the ids of the
        %             physical groups the curves are assigned to. Each curve
        %             can only be assigned to one physical group.
        %         -) nodes.....number_of_curves x number_of_element_nodes matrix
        %             of integer values containing the ids of the global nodes 
        %             corresponding to the elements local nodes. 
        %             For the node ordering, see 
        %             http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
        %             (05.07.2019)
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            
            switch number_of_element_nodes
                case 2
                    element_type_ID = 1;
                    finite_element = FirstOrderTriangleElement;
                case 3
                    element_type_ID = 8;
                    finite_element = SecondOrderTriangleElement;
                case 4
                    element_type_ID = 26;
                    finite_element = ThirdOrderTriangleElement;
            end
            
            curves = element_data(element_data(:, 2) == element_type_ID, :);
            [number_of_curves, ~] = size(curves);
            number_of_tags = curves(:, 3);
            node_columns_start = 4 + number_of_tags;
            
            col_idx = (repmat(node_columns_start, 1, number_of_element_nodes) + ...
                repmat(0 : number_of_element_nodes - 1, number_of_curves, 1));
            
            col_idx = [ones(number_of_curves, 1), ...
                4*ones(number_of_curves, 1), col_idx];
            
            col_idx = col_idx'; col_idx = col_idx(:);
            
            
            number_of_columns_per_element = number_of_element_nodes + 2;
            
            row_idx = repmat(1:number_of_curves, number_of_columns_per_element, 1);
            row_idx = row_idx(:);
            
            curve_data_idx = sub2ind(size(curves), row_idx, col_idx);
            curve_data_tmp = curves(curve_data_idx);
            curve_data_tmp = reshape(curve_data_tmp, ...
                number_of_columns_per_element, [])';
            
            
            % For second and third order elements, the node ordering must be
            % corrected as Gmsh uses a different ordering than this program.
            % See http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering for further
            % information
            %
            % The node orderings used in this program are defined in the
            % corresponding classes as comment headers.
            %
            if number_of_element_nodes > 2
                node_ordering = finite_element.curve_node_ordering;
                curve_data_tmp(:, 2 + node_ordering(:, 2)) = ...
                    curve_data_tmp(:, 2 + node_ordering(:, 1));
            end
            
            curve_data = struct( ...
                'number_of_curves', number_of_curves, ...
                'IDs', curve_data_tmp(:, 1), ...
                'physical_groups', curve_data_tmp(:, 2), ...
                'nodes', curve_data_tmp(:, 3:end));
            
            
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
        
        function nodes_on_boundary = get_nodes_on_boundary(boundary_conditions, ...
                group_data, curve_data, node_data)
            
        % function nodes_on_boundary = get_nodes_on_boundary(boundary_conditions, ...
        %         group_data, curve_data, node_data)
        %
        % DESCRIPTION: Function for determining the nodes and their corresponding
        %     boundary values on a set of boundary conditions.
        %
        %
        %
        % INPUT:
        %     boundary_conditions.....Structure containing information about a set of
        %         boundary conditions. Must have an equal form to the structures
        %         returned by parse_settings_file()
        %     group_data.....Structure containing information about the physical
        %         groups. Returned by parse_physical_groups().
        %     curve_data.....Structure containing information about the curve
        %         elements of the FEM mesh. Returned by extract_curve_elements().
        %     node_data.....Structure containing information about the nodes. 
        %         Returned by parse_nodes(). 
        %
        % OUTPUT:
        %     nodes_on_boundary.....Structure containing information about the nodes
        %         on the current boundary and their boundary values. The structure
        %         contains the following members:
        %         -) IDs.....IDs of nodes on the boundary
        %         -) values.....Boundary condition values of nodes
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            
            number_of_nodes = node_data.number_of_nodes;
            number_of_boundary_conditions = ...
                boundary_conditions.number_of_boundary_conditions;
            
            boundary_condition_names = boundary_conditions.names;
            boundary_condition_values = boundary_conditions.values;
            
            group_names = group_data.names;
            group_ids = group_data.tags;
            
            node_values_on_boundary = -ones(number_of_nodes, 1);
            
            for k = 1 : number_of_boundary_conditions
                
                % Get name and value of current boundary condition
                name_of_current_boundary_condition = ...
                    boundary_condition_names{k};
                
                value_of_current_boundary_condition = ...
                    boundary_condition_values(k);
                
                % Get the index of the physical group corresponding to the current
                % boundary condition
                idx_of_boundary_condition_group = find(...
                    strcmpi(name_of_current_boundary_condition, group_names));
                
                % Get the id of the physical group of the current boundary condition
                id_of_boundary_condition_group = ...
                    group_ids(idx_of_boundary_condition_group);
                
                % Get the physical curves assigned to the current boundary condition
                curves_of_current_boundary_condition = ...
                    curve_data.physical_groups == id_of_boundary_condition_group;
                
                % Get the nodes of the curves on the current boundary condition
                nodes_of_boundary_curves = ...
                    curve_data.nodes(curves_of_current_boundary_condition, :);
                
                % Remove redundant nodes
                nodes_of_boundary_curves = unique(nodes_of_boundary_curves(:));
                
                node_values_on_boundary(nodes_of_boundary_curves) = ...
                    value_of_current_boundary_condition;
                
                msg = sprintf('Nodes on boundary "%s": \n\t%s', ...
                    name_of_current_boundary_condition, ...
                    sprintf('%d ', nodes_of_boundary_curves));
                Misc.print_message(msg);
            end
            
            ids_of_nodes_on_boundary = int32(find(node_values_on_boundary ~= -1));
            
            nodes_on_boundary = struct(...
                'IDs', ids_of_nodes_on_boundary, ...
                'values', node_values_on_boundary(ids_of_nodes_on_boundary));  
        end
        
        function neumann_boundary_data = get_triangle_edges_on_neumann_boundary(...
                neumann_boundary_conditions, ...
                group_data, curve_data, node_data, triangle_data)
            
        % function neumann_boundary_data = get_triangle_edges_on_neumann_boundary(...
        %        neumann_boundary_conditions, ...
        %        group_data, curve_data, node_data, triangle_data)
        %
        % DESCRIPTION: Function for determining the edges of a triangle on the
        %     neumann boundary.
        %
        %     Neumann boundary conditions lead to curve integrals when calculating
        %     the element equation systems. Therefore it is requires to know which
        %     edge of the triangle belongs to the neumann boundary.
        %
        %     The function uses get_nodes_on_boundary() to determine the global nodes
        %     on the neumann boundary and searches for triangles whose local nodes
        %     are assigned to the boundary nodes. 
        %     
        %     Depending on which nodes of the triangles are located on the boundary,
        %     the algorithm can determine the correct triangle edge on the boundary.
        % 
        %     The algorithm numbers the triangle edges as following:
        %         -) The edge on the zeta-axis (the horizontal axis) is edge '1'.
        %         -) The edge on the eta-axis (the vertical axis) is edge '2'.
        %         -) The remaining edge is edge '3'.
        %
        %
        %
        % INPUT:
        %     neumann_boundary_conditions.....Structure containing information about
        %         the neumann boundary conditions.
        %     group_data.....Structure containing information about the physical
        %         groups. Returned by parse_physical_groups().
        %     curve_data.....Structure containing information about the curve
        %         elements of the FEM mesh. Returned by extract_curve_elements().
        %     node_data.....Structure containing information about the nodes. 
        %         Returned by parse_nodes(). 
        %
        % OUTPUT:
        %     neumann_boundary_data.....Structure containing the IDs and edges of
        %         triangles on the neumann boundary and the corresponding 
        %         boundary condition values. The structure has following members:
        %         -) IDs.....IDs of triangles on the neumann boundary
        %         -) edges.....Edges of the corresponding triangles on the boundary
        %         -) values.....Boundary condition values
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            
            neumann_boundary_data_tmp = cell(1,3);
            
            boundary_condition_names = neumann_boundary_conditions.names;
            boundary_condition_values = neumann_boundary_conditions.values;
            
            number_of_boundary_conditions = ...
                neumann_boundary_conditions.number_of_boundary_conditions;
            
            % Get matrix of element nodes
                nodes_of_elements = triangle_data.nodes;
            
            for k = 1 : number_of_boundary_conditions
                
                % Create temporary structure to get nodes of current boundary
                % condition
                current_boundary_condition = struct( ... 
                    'number_of_boundary_conditions', 1, ...
                    'names', {boundary_condition_names(k)}, ...
                    'values', boundary_condition_values(k));
                
                
                % 'tmp' contains now IDs and values of the nodes on the current 
                % neumann boundary. 
                tmp = Parser.get_nodes_on_boundary( ...
                    current_boundary_condition, group_data, curve_data, node_data);
                
                % Vector containing the nodes of the current boundary
                nodes_of_boundary = tmp.IDs;
                
                % Vector containing the values of the current boundary (all the
                % same)
                values = tmp.values;
                clear tmp;
                
                
                % ===== Determination of triangle edges on boundary
                % From the algorithm above we now know the nodes on the current 
                % boundary. We know search for triangles which contain the same nodes. 
                %
                
                % 'nodes_of_elements' contain the nodes of each triangle as its rows.
                % The algorithm now checks for each element of the matrix if it is
                % equal to one node on the boundary.
                %
                % The result is a binary matrix indicating which triangle has nodes
                % local nodes on the current boundary.
                % 
                % E.g. in the matrix below the triangles 1 and 3 have nodes on the 
                % boundary.
                %
                %[ 0, 1, 1;
                %  0, 0, 0;
                %  1, 0, 1;
                %  0, 0, 0 ]
                %
                element_nodes_on_boundary = ...
                    ismember(nodes_of_elements, nodes_of_boundary);
                
                % 'elements_on_boundary' is a vector of binary values with logic 1s
                % for triangles on the current boundary.
                % 
                % E.g. for the matrix from above: [1; 0; 1; 0]
                % (Elements 1 and 3 have nodes on the boundary)
                elements_on_boundary = sum(element_nodes_on_boundary,2) > 0;
                
                % Remove triangles from the matrix without nodes on the boundary
                element_nodes_on_boundary = ...
                    element_nodes_on_boundary(elements_on_boundary, :);
                
                % Using the node ordering specified in the Gmsh manual, the algorithm
                % can determine the edge of a triangle on the boundary.
                % 
                % E.g. for a first order triangle whose nodes 2 and 3 are on the
                % boundary, the algorithm assigns edge 3 as edge on the boundary.
                edges_on_boundary = ...
                    Parser.assign_nodes_on_neumann_boundary_to_triangle_edges(...
                    triangle_data.IDs(elements_on_boundary), element_nodes_on_boundary);
                
                [N, ~] = size(edges_on_boundary);
                
                % Append ids of elements
                neumann_boundary_data_tmp{1} = [neumann_boundary_data_tmp{1}; ...
                    edges_on_boundary(:, 1)];
                
                % Append element sides
                neumann_boundary_data_tmp{2} = [neumann_boundary_data_tmp{2}; ...
                    edges_on_boundary(:, 2)];
                
                % Append boundary condition values
                neumann_boundary_data_tmp{3} = [neumann_boundary_data_tmp{3}; ...
                    values(1) * ones(N, 1)];
            end
            
            neumann_boundary_data = struct(...
                'IDs', neumann_boundary_data_tmp{1}, ...
                'edges', neumann_boundary_data_tmp{2}, ...
                'values', neumann_boundary_data_tmp{3});
        end
        
        function edges_on_boundary = ...
                assign_nodes_on_neumann_boundary_to_triangle_edges(...
                boundary_element_ids, element_nodes_on_boundary)
            
            edges_on_boundary = [];
            
            [~, number_of_element_nodes] = size(element_nodes_on_boundary);
            
            if number_of_element_nodes == 3
                nodes_of_triangle_edges = ...
                    FirstOrderTriangleElement.nodes_of_triangle_edges;
            elseif number_of_element_nodes == 6
                nodes_of_triangle_edges = ...
                    SecondOrderTriangleElement.nodes_of_triangle_edges;
            elseif number_of_element_nodes == 10
                nodes_of_triangle_edges = ...
                    ThirdOrderTriangleElement.nodes_of_triangle_edges;
            end
            
            
            for k = 1 : 3
                nodes_of_current_edge = nodes_of_triangle_edges(k, :);
                
                idx = sum(element_nodes_on_boundary == nodes_of_current_edge, 2) == ...
                    number_of_element_nodes;
                
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
            material_and_source_properties_tmp = zeros(number_of_elements, 2);
            
            % Names of physical groups
            group_names = group_data.names;
            
            % Names of regions. A 'region' is a physical group with constant material
            % properties and a constant source value.
            region_names = region_data.names;
            
            idx = ismember(group_names, region_names);
            group_ids_of_regions = group_data.tags(idx);
            
            
            for k = 1 : length(region_names)
                current_group_id = group_ids_of_regions(k);
                current_material_property = region_data.material_values(k);
                current_source_value = region_data.source_values(k);
                
                ids_of_elements = triangle_data.physical_groups == current_group_id;
                
                % Number of elements (triangles) belonging to the current group
                N = sum(ids_of_elements);
                
                
                material_and_source_properties_tmp(ids_of_elements, :) = ...
                    repmat([current_material_property, current_source_value], ...
                    N, 1);
            end
            
            % IDs of triangles
            material_IDs = triangle_data.IDs;
            
            % Not each triangle contains a source. Extract triangles with sources and
            % store triangle ids and source values to struct.
            triangles_with_sources_bin = material_and_source_properties_tmp(:, 2) ~= 0;
            
            triangles_with_sources = triangle_data.IDs(triangles_with_sources_bin);
            source_values = ...
                material_and_source_properties_tmp(triangles_with_sources_bin, 2);
            
            material_and_source_properties = struct( ...
                'material_IDs', material_IDs, ...
                'material_values', material_and_source_properties_tmp(:, 1), ...
                'triangles_with_sources', triangles_with_sources, ...
                'source_values', source_values);
        end
        
        function print_statistics(node_data, triangle_data)
            
            number_of_nodes = node_data.number_of_nodes;
            number_of_elements = triangle_data.number_of_triangles;
            
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
