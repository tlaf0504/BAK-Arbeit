classdef Solver
    methods(Static)
        
        function vec = get_node_index_to_equation_index_mapping(...
                number_of_nodes, nodes_with_dirichlet_boundary_conditions)
            
            vec = -ones(number_of_nodes,1);
            cnt = 1;
            
            for k = 1 : number_of_nodes
                if ~ismember(k, nodes_with_dirichlet_boundary_conditions)
                    vec(k) = cnt;
                    cnt = cnt + 1;
                end
            end
        end
        
        
        
        function success= solve(problem_location)
            tmp = pwd;
            cd(problem_location);
            
            % Plot section separator to console
            Misc.print_message(sprintf('%s\n%s\n%s', ...
                Misc.console_section_separator, ...
                'Solver', Misc.console_section_separator));
            
            % Check if file 'problem_setup.mat' exists
            % This contains all neccessary information about the problem.
            if ~Misc.check_file_existence(problem_location, 'problem_setup.mat')
                msg = sprinft(['Error. Setup-file "problem setup.mat" in folder ', ...
                    '%s not found.'], problem_location);
                Misc.print_error_message(msg);
            end
            load('problem_setup.mat', 'problem_setup');
            

            % Check if the meshing parsing procedure was completed.
            % problem_setup.state must be higher than
            % Setup.setup_state_meshing_finished
            if problem_setup.state < Setup.setup_state_parsing_finished
                msg = sprintf(['Error. Parsing procedure was not ', ...
                    'finished. Maybe the setup procedure was interrupted.', ...
                    'Please rerun the complete parsing procedure.']);
                Misc.print_error_message(msg);
                success = 0;
                return
            end
            
             % Check if file 'internal/mesh_data.mat' exists. This file contains all
            % information about the FEM mesh
            mesh_data_location = fullfile(problem_location, 'internal');
            if ~Misc.check_file_existence(mesh_data_location, 'mesh_data.mat')
                msg = sprinft(['Error. Mesh data file "mesh_data.mat" in folder ', ...
                    '%s not found.'], mesh_data_location);
                Misc.print_error_message(msg);
            end
            load(fullfile(mesh_data_location, 'mesh_data.mat'), 'mesh_data');
            

            % Assemble equation system
            Misc.print_message('Assembling equation system...')
            tic
            
            [A, r, result_to_global_node_mapping, success] = Solver.assemble(...
                mesh_data.node_data, ...
                mesh_data.triangle_data, ...
                mesh_data.material_and_source_properties, ...
                mesh_data.dirichlet_boundary_data, ...
                mesh_data.neumann_boundary_data, ...
                problem_setup.problem_type, ...
                7,... % Plane integration order
                3, ...% Curve integration order
                problem_setup.mesh_order);
            time = toc;
            Misc.print_message(sprintf('Done\nElapsed time is %d seconds.\n', time));
            
            if ~success
                return
            end
            
            
            Misc.print_message('Solving equation system...')
            tic
            unknowns = A\r;
            time = toc;
            
            Misc.print_message(sprintf('Done\nElapsed time is %d seconds.\n', time));
            
             Misc.print_message(sprintf('Writing data to "%s"...\n', ...
                fullfile(problem_setup.problem_location, 'results', 'results.mat')));
            
            save(fullfile(problem_location, 'results', 'results.mat'), 'unknowns', ...
                'result_to_global_node_mapping')
            
            Misc.print_message(Misc.console_section_separator);
            Misc.print_message('\n');
            
            problem_setup.state = Setup.setup_state_solving_finished;
            problem_setup.result_file = 'results.mat';
            problem_setup.result_file_hash = DataHash(fullfile( ...
                'results', 'results.mat'));
            Setup.update_problem_setup_file(problem_setup);
            
            
           
            
            cd(tmp);
        end
        
        
        function [A, r, result_to_global_node_mapping, success] = assemble(...
                node_data, triangle_data, ...
                material_and_source_properties, dirichlet_boundary_data, ...
                neumann_boundary_data, problem_type_number, ...
                plane_integration_order, curve_integration_order, ...
                mesh_order)
            
            % function [A, r, result_to_global_node_mapping] = assemble(...
            %    node_data, triangle_data, ...
            %    material_and_source_properties, dirichlet_boundary_data, ...
            %    neumann_boundary_data, problem_type, ...
            %    plane_integration_order, curve_integration_order, ...
            %    finite_element_type)
            %
            % DESCRIPTION:
            %
            % Will be added later
            
            
            
            % Assignment of finite element class
            finite_element_type = ...
                Misc.get_finite_element_class_from_mesh_order(mesh_order);
            
            % Assignment of problem type
            [problem_type, success] = Misc. ...
                get_problem_type_class_from_problem_type_number( ...
                problem_type_number);
            if ~success
                return
            end
            
            
            % Column 1 contains id of element and column 2 the id of the
            % corresponding physical group. This information can be ignored.
            node_mapping = triangle_data.nodes;
            N_glob_nodes = node_data.number_of_nodes;
            
            % Get number of elements (N_e) and number of element nodes (N)
            [N_e, N] = size(node_mapping);
            
            % Second column contains node coordinates
            global_node_coordinates = node_data.coordinates;
           
            global_nodes_on_dirichlet_boundary = dirichlet_boundary_data.IDs;
            
            % Number of nodes on dirichlet boundary
            N_dir = length(global_nodes_on_dirichlet_boundary);
            
            % Number of unknowns
            N_unknown = N_glob_nodes - N_dir;
            
            % Coefficient matrix and right-side vector of global equation system
            A = sparse(N_unknown, N_unknown);
            r = zeros(N_unknown, 1);
            
            % Material constants
%             vacuum_material = Misc.get_vacuum_material(problem_type_number);
            vacuum_material = 1;
            
            % Extract material and source properties
            triangles_with_sources = ...
                material_and_source_properties.triangles_with_sources;
            source_values = material_and_source_properties.source_values;

            
            for element_id = 1 : N_e
                local_node_mapping = node_mapping(element_id, :);
                local_node_mapping = local_node_mapping(:);
                
                xe = global_node_coordinates(local_node_mapping, 1);
                ye = global_node_coordinates(local_node_mapping, 2);
               
                material_x = vacuum_material * material_and_source_properties.material_values(element_id);
                material_y = vacuum_material * material_and_source_properties.material_values(element_id);
                
                % Local coefficient matrix and right-side vector
                A_loc = zeros(N, N);
                r_loc = zeros(N, 1);
                
                % ===== DEBUG, REMOVE LATER =====
                A_loc_tmp = zeros(N, N);
                
                row_idx = 1 : N;
                row_idx = row_idx(:);
                
                col_idx = 1 : N;
                col_idx = col_idx(:);
                
                % ===== Dirichlet boundary conditions
                
                % Vector of logical values indicating which global node(s) on the
                % dirichlet boundary match with local nodes of the current element.
                % Dimenison : N_dir (Number of dirichlet boundary conditions) x 1
                %
                % e.g. nodes_on_dirichlet_boundary = [1;5;3;9] and 
                % local_node mapping = [1,10,15] then 
                % dirichlet_conditions_of_current_element = [1;0;0;0].
                dirichlet_conditions_of_currnet_element_bin_idx = ...
                    ismember(global_nodes_on_dirichlet_boundary, local_node_mapping);
                
                % Vector of logical values indicating which local node is located on
                % the dirichlet boundary.
                % Dimension: number_of_element_nodes x 1
                %
                % Using the example from above: 
                % local_nodes_with_dirichlet_conditions = [1,0,0];
                local_nodes_on_dirichlet_boundary_bin_idx = ...
                    ismember(local_node_mapping, global_nodes_on_dirichlet_boundary);
                
                % IDs of local nodes on the dirichlet boundary.
                % On entry for each local node on the dirichlet boundary
                %
                % E.g. if node 1 and 3 lie on the dirichlet boundary:
                % local_nodes_on_dirichlet_boundary = [1;3];
                %
                local_nodes_on_dirichlet_boundary = ...
                    row_idx(local_nodes_on_dirichlet_boundary_bin_idx);
                
                % Given values for local nodes on dirichlet boundary.
                % One entry for each local node on the dirichlet boundary.
                %
                % E.g. using the example from above with 10V for node 1 and 0V for
                % node 3:
                % values_of_local_nodes_on_dirichlet_boundary = [10;3]
                %
                values_of_local_nodes_on_dirichlet_boundary = ...
                    dirichlet_boundary_data.values(dirichlet_conditions_of_currnet_element_bin_idx);
                
                % ===== Neumann boundary conditions
                
                % Entries for current element in members of neumenn_boundary_data
                % struct.
                % 
                % The struct contains 3 members of the same size. Each member
                % contains one entry for each neumann boundary condition. 
                % (See Parser.get_triangle_edges_on_neumann_boundary() for further
                % information)
                % 
                % e.g. 
                % .IDs = [1,5,6,1];
                % .edges = [1,1,1,2];
                % .values = [5,5,5,3];
                %
                % The 'IDs' member contains the IDs of the triangles having edges on
                % the neumenn boundary.
                %
                % The 'edges' member contains the edge of the triangle on the neumann 
                % boundary. See Parser.get_triangle_edges_on_neumann_boundary() how
                % the edges are defined.)
                %
                % The 'values' member contain the values for the single boundary
                % conditions as defined in the settings file.
                %
                % As a triangle contains three sides, each triangle can possibly have
                % 3 neumann boundary conditions. As shown in the example, multiple
                % entries for the same triangle indicate multiple neumann boundary
                % conditions for each element.
                
                % Get indices of the entries for the current element
                neumann_boundary_data_entries_for_current_element = ...
                    neumann_boundary_data.IDs == element_id;
                
                % Get the edges of the current element on the neumann boundary
                triangle_edges_of_current_element_on_neumann_boundary = ...
                    neumann_boundary_data.edges( ...
                    neumann_boundary_data_entries_for_current_element);
                
                % Get the boundary values for the current element
                neumann_boundary_values_for_current_element = ...
                    neumann_boundary_data.values( ...
                    neumann_boundary_data_entries_for_current_element);
                
                number_of_neumann_conditions_for_current_element = ...
                    length(neumann_boundary_values_for_current_element);
                    
                

                % A dirichlet boundary condition at local node k removes the k-th 
                % row(the k-th equation) of the element-equation system. The elements 
                % of A_loc will be zero in this row.
                row_idx(local_nodes_on_dirichlet_boundary) = [];
                
                % Due to line above, the algorithm will only iterate above the 
                % elements without dirichlet boundary conditions.
                for l = 1 : length(row_idx)
                    
                    row = row_idx(l);
                    
                    % ===== Calculation of right-side vector elements
                    
                    % === Plane integral
                    % Only needed for electrostatic problems, not for static current
                    % problemes
                    %
                    % Integration also not needed if no free charges occurr in the
                    % current element
                    if  ismember(element_id, triangles_with_sources)
                        fun = @(zeta, eta) problem_type. ...
                            right_side_plane_integrant(finite_element_type, xe, ...
                            ye, row, zeta, eta, ...
                            source_values(triangles_with_sources == element_id));
                        
                        r_loc(row) = GaussianQuadrature. ...
                            integrate_2D_normalized_triangle_region(fun,...
                            plane_integration_order);
                    end
                    
                    
                   
                    % === Neumann boundary integrals
                    for b = 1: number_of_neumann_conditions_for_current_element
                        edge = triangle_edges_of_current_element_on_neumann_boundary(b);
                        value = neumann_boundary_values_for_current_element(b);
                        
                        % Integrant for curve integration over triangle edge
                        fun = @(t) problem_type. ...
                            right_side_neumann_integrant(finite_element_type, ...
                            xe, ye, row, t, edge, value);
                        
                        r_loc(row)  = GaussianQuadrature.integrate_1d_line(fun, ...
                            curve_integration_order, 1, 0) + r_loc(row);
                    end
                    
                    
                     % ===== Calculation of coefficient matrix of element equation 
                     % system =====
                    
                    for m = 1 : length(col_idx)
                        
                        col = col_idx(m);
                        

                        fun = @(zeta, eta) problem_type. ...
                            element_matrix_integrant(...
                            finite_element_type, xe, ye, row , col, zeta, ...
                            eta, material_x, material_y);
                        
                        coef = GaussianQuadrature. ...
                            integrate_2D_normalized_triangle_region(fun,...
                            plane_integration_order);
                        
                        % If the k-th local node contains a dirichlet boundary 
                        % condition, the coeffient in the k-th column multiplied with 
                        % the value of the boundary condition is brought to the right 
                        % side of the current equation
                        if ismember(col, local_nodes_on_dirichlet_boundary)
                            idx = (local_nodes_on_dirichlet_boundary == col);
                            r_loc(row) = r_loc(row) - coef * ...
                                values_of_local_nodes_on_dirichlet_boundary(idx);
                        else
                            A_loc(row, col) = coef;
                        end
                        
                        % ===== DEBUG, REMOVE LATER =====
                        A_loc_tmp(row, col) = coef;
                        
                        
                    end
                end
                
                % ===== Add current local equation system to global equation 
                %   system =====
                
                global_node_to_equation_system_index_mapping = ...
                    Solver.get_node_index_to_equation_index_mapping( ...
                    N_glob_nodes, dirichlet_boundary_data.IDs);
                
                result_to_global_node_mapping = find( ...
                    global_node_to_equation_system_index_mapping ~= -1);
                
                
                % If the k-th local node has a dirichlet boundary condition, the 
                % elements from the k-th column should not be included in the global 
                % equation system as they are already condidered in the right-side 
                % vector.
                col_idx(local_nodes_on_dirichlet_boundary) = [];
                
                % global_node_to_equadion_system_mapping maps the global nodes to 
                % their corresponding row and column within the global equation 
                % system.
                % E.g. the unknown potential of global node 6 maps to the 5th 
                % row/column of the global equation system as, for example, the 3rd 
                % global node has a dirichlet boundary condition.
                %
                % loca_node_mapping maps the local node to the global nodes and 
                % row_idx and col_idx contains the indices of the equation system to 
                % copy.
                global_row_idx = global_node_to_equation_system_index_mapping(...
                    local_node_mapping(row_idx));
                global_col_idx = global_node_to_equation_system_index_mapping(...
                    local_node_mapping(col_idx));
                
 
                A(global_row_idx, global_col_idx) = A(global_row_idx, global_col_idx) + ...
                    A_loc(row_idx, col_idx);
                
                r(global_row_idx) = r(global_row_idx) + r_loc(row_idx);
                
               
            end
        end
       

    end
end