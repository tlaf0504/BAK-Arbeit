classdef Parser
    
    properties
        Property1
    end
    
    methods(Static)
        function [success, problem_setup] = parse(location, geometry_filename)
            
            success = 1;

            
            %fprintf('Parsing mesh file ...\n')
            Misc.print_message('Parsing mesh file ...')
            
            tmp = pwd;
            cd(location);
            problem_setup = Setup.new_problem_setup(location, geometry_filename);
            Parser.parse_mesh_file(problem_setup);
            cd(tmp);
            
            
            
        end
        
        
        
        
        function [number_of_nodes, node_data, number_of_triangles, triangle_data, ...
                number_of_groups, group_data] = parse_mesh_file(problem_setup)
            file = fopen(problem_setup.mesh_file, 'r');
            
            content = fileread(problem_setup.mesh_file);
            fclose(file);
            
            
            [number_of_nodes, node_data] = Parser.parse_nodes(content);
            [~, element_data] = Parser.parse_elements(content);
            [number_of_groups, group_data] = Parser.parse_physical_groups(content);
            
            clear content;
            
            % Hier morgen weiter.
            [number_of_triangles, triangle_data] = Parser.extract_triangle_elements(3, ...
                element_data);
            
        end
        
        
        function [number_of_nodes, node_data] = parse_nodes(file_content)
            
            regex = '\$Nodes(.*)\$EndNodes';
            [start_idx, end_idx] = regexp(file_content, regex, 'once');
            
            file_content = file_content(start_idx : end_idx);
            file_content = replace(file_content, '$Nodes', '');
            
            [number_of_nodes, ~, ~, idx] = sscanf(file_content, '%d', 1);
            
            node_data = textscan(file_content(idx+1:end), '%d %f %f %f', ...
                number_of_nodes);
            
            % First element contains vector of node indices, second column contains
            % matrix of node coordinates. ( Format of each matrix row: [x,y,z])
            node_data = {node_data{1}, [node_data{2}, node_data{3}, node_data{4}]};
        end
        
        
        
        function [number_of_elements, element_data] = parse_elements(file_content)
            
            
            regex = '\$Elements(.*)\$EndElements';
            [start_idx, end_idx] = regexp(file_content, regex, 'once');
            
            file_content = file_content(start_idx : end_idx);
            file_content = replace(file_content, '$Elements', '');
            
            [number_of_elements, ~, ~, idx] = sscanf(file_content, '%d', 1);
            
            element_data = textscan(file_content(idx:end), repmat('%d ', 1, 30), ...
                number_of_elements);
            
            element_data = cell2mat(element_data);
        end
        
        function [number_of_groups, group_data] = parse_physical_groups(file_content)
            regex = '\$PhysicalNames(.*)\$EndPhysicalNames';
            [start_idx, end_idx] = regexp(file_content, regex, 'once');
            
            if isempty(start_idx) && isempty(end_idx)
                fprintf('Mesh file contains no physical groups.\n\r')
                group_data = 0;
                number_of_groups = 0;
                return
            end
            
            file_content = file_content(start_idx : end_idx);
            file_content = replace(file_content, '$PhysicalNames', '');
            
            [number_of_groups, ~, ~, idx] = sscanf(file_content, '%d', 1);
            
            group_data = textscan(file_content(idx:end), '%d %d %s', ...
                number_of_groups);
            
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
            
            
            col_idx = (repmat(node_columns_start, 1, number_of_nodes) + ...
                repmat(0 : number_of_nodes - 1, number_of_triangles, 1))';
            col_idx = col_idx(:);
            row_idx = repmat(1:number_of_triangles, number_of_nodes, 1);
            row_idx = row_idx(:);
            
            triangle_data_idx = sub2ind(size(triangles), row_idx, col_idx);
            triangle_data = triangles(triangle_data_idx);
            triangle_data = reshape(triangle_data, number_of_nodes, [])';
        end
    end
end

