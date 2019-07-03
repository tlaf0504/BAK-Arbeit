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
        
        function [A, r, result_to_global_node_mapping] = assemble(...
                node_data, triangle_data, ...
                material_properties, dirichlet_boundary_conditions, ...
                neumann_boundary_conditions, sources, problem_type, ...
                plane_integration_order, curve_integration_order, ...
                finite_element_type)
            
            % function [A, r] = assemble(global_node_coordinates, node_mapping, ...
            %    material_properties, dirichlet_boundary_conditions, ...
            %    neumann_boundary_conditions, volume_charges, problem_type, ...
            %    plane_integration_order, curve_integration_order)
            %
            % DESCRIPTION:
            % Function for assembling the global coordinate system for an
            % electrostatic or a static current problem.
            %
            % The following definitions are used in the further descripiton:
            % N_nodes.....Number of global nodes
            % N_e.....Number of finite elements
            % N.....Number of element nodes
            % N_dir.....Number of dirichlet boundary conditions
            % N_neu.....Number of neumann boundary conditions
            % N_roh.....Number of elements with free volume charges
            % 
            % INPUT:
            % global_node_coordinates.....N_nodes x 2 matrix containing the
            %   x-coordinates of the globla nodes in the first column and the
            %   y-coordinates in the second column.
            % node_mapping.....N_e x N matrix containing the 
            % material_properties.....N_e x 2 matrix contianing the material
            %   propertinode_plot_data.dirichlet_boundary_conditions(:, 1)es of each element. The first column contains the material
            %   property in x-direction, the second column for the y-direction.
            % dirichlet_boundary_conditions.....1x2 cell array where the first
            %   elements contains a nx1 matrix of integer values, where n corresponds
            %   containing the ids of the nodes on the dirichlet boundary.
            %   The second entry in the cell array conains an nx1 matrix of double
            %   values whith the corresponding boundary condition values.
            % neumann_boundary_conditions.....1 x 3 cell array containing the neumann
            %   boundary conditions of the problem.
            %   The first entry contains the ids of elements on the neumann boundary,
            %   the second entry contains the sides of the triangles on the boundary 
            %   and the thrid elements contains the actual boundary condition value.
            %   Multiple entries are used if one element contains multiple boundary
            %   conditions. (e.g. the two sides of the triangle lie on the boundary)
            % sources.....N_roh x 2 matrix containing information about the
            %   sources (volume charges/ current densities etc.) within the problem 
            %   area. Each row stands for one finite element containing a source. 
            %   The first column contains the number of the finite element and the 
            %   second column contains the source value.
            % problem_type.....Either 'E-Static' or 'Static current' (not case
            %   sensitive) specifiying if either an electrostatic problem or a static
            %   current problem should be solved.
            % plane_integration_order.....The integration order (number of sampling 
            %   points) used for gaussian plane integration. Possible vaues: 3, 7, 15
            %  curve_integration_order.....The integration order used for curve
            %  integration. Values from 1 to 6 are possible.
            % finite_element_type.....Type of the finite elements used for the
            %   problem. Either 'First order', 'Second order', 'Third order' (not 
            %   case sensitive) are possible.
            %
            %
            
            if strcmpi(finite_element_type, 'First order')
                error('First order triangles are not implemented yet.')
                
            elseif strcmpi(finite_element_type, 'Second order')
                finite_element_type = SecondOrderTriangleElement;
                
            elseif strcmpi(finite_element_type, 'Third order')
                error('Third order triangles are not implemented yet.')
                
            else
                error(['Invalid finite element type specified. Please use ', ...
                    'either "First order", "Second order" or "Third order"'])
            end
            
            % Column 1 contains id of element and column 2 the id of the
            % corresponding physical group. This information can be ignored.
            node_mapping = triangle_data(:, 3:end);
            
            % Get number of elements (N_e) and number of element nodes (N)
            [N_e, N] = size(node_mapping);
            
            % Second column contains node coordinates
            global_node_coordinates = node_data{2};
            % Get number of global nodes
            [N_glob_nodes, ~] = size(global_node_coordinates);
            
            nodes_on_dirichlet_boundary = dirichlet_boundary_conditions{1};
            dirichlet_boundary_values = dirichlet_boundary_conditions{2};
            N_dir = length(nodes_on_dirichlet_boundary);
            
            % Number of unknown potentials
            N_unknown = N_glob_nodes - N_dir;
            
            % Coefficient matrix and right-side vector of global equation system
            A = sparse(N_unknown, N_unknown);
            r = zeros(N_unknown, 1);

            
            for e = 1 : N_e
                local_node_mapping = node_mapping(e, :);
                local_node_mapping = local_node_mapping(:);
                
                xe = global_node_coordinates(local_node_mapping, 1);
                ye = global_node_coordinates(local_node_mapping, 2);
                
                if strcmpi(problem_type, 'E-Static')
                    % Vacuum permittivity
                    vacuum_material = 8.854187812813e-12;
                elseif strcmpi(problem_type, 'Static current')
                    % Magnetic vacuum permeability
                    vacuum_material = 4*pie-7;
                end
                material_x = vacuum_material * material_properties(e, 1);
                material_y = vacuum_material * material_properties(e, 2);
                
                % Local coefficient matrix and right-side vector
                A_loc = zeros(N, N);
                r_loc = zeros(N, 1);
                
                row_idx = 1 : N;
                row_idx = row_idx(:);
                
                col_idx = 1 : N;
                col_idx = col_idx(:);
                
                % Vector containing logical ones at dirichlet boundary conditions of 
                % current element.
                % Dimenison : N_dir (Number of dirichlet boundary conditions) x 1
                dirichlet_conditions_of_currnet_element_bin_idx = ...
                    ismember(nodes_on_dirichlet_boundary, local_node_mapping);
                
                % Vector containing logical ones at local nodes with dirichlet 
                % boundary conditions.
                % Dimenstion N (Number of element nodes) x 1
                local_nodes_with_dirichlet_conditions_bin_idx = ...
                    ismember(local_node_mapping, nodes_on_dirichlet_boundary);
                
                % First column: Local nodes with dirichlet boundary conditions
                % Second and third column: Global node + boundary conditions 
                % corresponding to local node in first column.
                %
                % E.g.:
                %
                % 1 25 10.0
                % 3 73 0.0
                %
                % Local node 1 corresponds to global node 25 with dirichlet boudnary 
                % condition of 10 (V, A, etc.). Local node 3 corresponds to global 
                % node 73 with a boundary condition of 0.
                %
                dirichlet_boundary_data = {...
                    row_idx(local_nodes_with_dirichlet_conditions_bin_idx), ...
                    dirichlet_boundary_conditions{1}(dirichlet_conditions_of_currnet_element_bin_idx), ...
                    dirichlet_boundary_conditions{2}(dirichlet_conditions_of_currnet_element_bin_idx)};
                
                
                % Get Neumann boundary conditions of element. It is a 1x3 cell array,
                % where each 'row' of the internal cells stands for one neumann 
                % boudnary condition of the current element. The rows are of the 
                % following form:
                %
                % Column 1: Side of the current element with a boundary condition.
                % e.g. for a first order triangle, the side
                % between local nodes 1 and 3 is side 1, between local nodes 2 and 3
                % side 2 etc.
                %
                % Column 2: Boundary condition value.
                
                boundary_condition_idx = neumann_boundary_conditions{1} ...
                    == e;
                if ~isempty(boundary_condition_idx)
                    
                    neumann_boundary_flag = 1; 
                    neumann_boundary_conditions_of_element = ...
                        {neumann_boundary_conditions{2}(boundary_condition_idx) ...
                        neumann_boundary_conditions{3}(boundary_condition_idx)};
                else
                    neumann_boundary_flag = 0;
                end

                
                
                % A dirichlet boundary condition at local node k removes the k-th 
                % row(the k-th equation) of the element-equation system. The elements 
                % of A_loc will be zero in this row.
                row_idx(local_nodes_with_dirichlet_conditions_bin_idx) = [];
                
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
                    if ~isempty(sources) && strcmpi(problem_type, 'E-Static') && ...
                            sources(e) ~= 0
                        fun = @(zeta, eta) ElectrostaticProblem. ...
                            right_side_plane_integrant(finite_element_type, xe, ...
                            ye, row, zeta, eta, sources(e));
                        
                        r_loc(row) = GaussianQuadrature. ...
                            integrate_2D_normalized_triangle_region(fun,...
                            plane_integration_order);
                    end
                    
                    
                    % === Neumann integral
                    if neumann_boundary_flag
                        
                        number_of_neumann_boundary_conditions = ...
                            length(neumann_boundary_conditions_of_element{1});
                        
                        for b = 1: number_of_neumann_boundary_conditions
                            side = neumann_boundary_conditions_of_element{1}(b);
                            value = neumann_boundary_conditions_of_element{2}(b);
                            
                            fun = @(t) ElectrostaticProblem. ...
                                right_side_neumann_integrant(finite_element_type, ...
                                xe, ye, row, t, side, value);
                            
                            r_loc(row)  = GaussianQuadrature.integrate_1d_line(fun, ...
                                curve_integration_order, 1, 0) + r_loc(row);
                        end
                    end
                    
                    
                    
                    
                    for m = 1 : length(col_idx)
                        
                        col = col_idx(m);
                        
                        % ===== Calculation of local coefficient matrix =====
                        
                        fun = @(zeta, eta) ElectrostaticProblem. ...
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
                        if ismember(col, dirichlet_boundary_data{1})
                            idx = (dirichlet_boundary_data{1} == col);
                            r_loc(row) = r_loc(row) - coef * ...
                                dirichlet_boundary_data{3}(idx);
                        else
                            A_loc(row, col) = coef;
                        end
                    end
                end
                
                % ===== Add current local equation system to global equation 
                %   system =====
                
                global_node_to_equation_system_index_mapping = ...
                    Solver.get_node_index_to_equation_index_mapping( ...
                    N_glob_nodes, dirichlet_boundary_conditions{1});
                
                result_to_global_node_mapping = find( ...
                    global_node_to_equation_system_index_mapping ~= -1);
                
                % If the k-th local node has a dirichlet boundary condition, the 
                % elements from the k-th column should not be included in the global 
                % equation system as they are already condidered in the right-side 
                % vector.
                col_idx(local_nodes_with_dirichlet_conditions_bin_idx) = [];
                
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
        
        function success= solve(problem_location)
            
            success = 1;
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
            
            
            material_and_source_properties = mesh_data.material_and_source_properties;
            
            
            % Assembler need N_elements x 2 matrix where each row contains the
            % material porperites in x- and y direction.
            % 
            % In this version of the program, only isotropic materials are supported.
            material_properties = [material_and_source_properties(:, 1), ...
                material_and_source_properties(:, 1)];
            
            elements_with_sources = find(material_and_source_properties(:, 2));
            
            source_data = [elements_with_sources, ...
                material_and_source_properties(elements_with_sources)];
            
            % Assemble equation system
            Misc.print_message('Assembling equation system...')
            tic
            
            % May be changed in future when other element typed are supported
            switch(problem_setup.mesh_order)
                case 1
                    order_string = 'First order';
                case 2
                    order_string = 'Second order';
                case 3
                    order_string = 'Third order';
            end
            
            [A, r, result_to_global_node_mapping] = Solver.assemble(...
                mesh_data.node_data, mesh_data.triangle_data, ...
                material_properties, mesh_data.dirichlet_boundary_data, ...
                mesh_data.neumann_boundary_data, ...
                source_data, 'E-Static', 7,3, order_string);
            time = toc;
            Misc.print_message(sprintf('Done\nElapsed time is %d seconds.\n', time));
            
            
            Misc.print_message('Solving equation system...')
            tic
            unknowns = A\r;
            time = toc;
            Misc.print_message(sprintf('Done\nElapsed time is %d seconds.\n', time));
            
            save(fullfile(problem_location, 'results', 'results.mat'), 'unknowns', ...
                'result_to_global_node_mapping')
            
            cd(tmp);
        end
    end
end