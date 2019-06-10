classdef Solver
    methods(Static)
        
        function vec = get_node_index_to_equation_index_mapping(...
                number_of_nodes, nodes_with_direchlet_boundary_conditions)
            
            vec = -ones(number_of_nodes,1);
            cnt = 1;
            
            for k = 1 : number_of_nodes
                if ~ismember(k, nodes_with_direchlet_boundary_conditions)
                    vec(k) = cnt;
                    cnt = cnt + 1;
                end
            end
        end
        
        function [A, r] = assemble(global_node_coordinates, node_mapping, ...
                material_properties, direchlet_boundary_conditions, ...
                neumann_boundary_conditions, volume_charges, problem_type, ...
                plane_integration_order, curve_integration_order, ...
                finite_element_type)
            
            % function [A, r] = assemble(global_node_coordinates, node_mapping, ...
            %    material_properties, direchlet_boundary_conditions, ...
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
            % N_dir.....Number of direchlet boundary conditions
            % N_neu.....Number of neumann boundary conditions
            % N_roh.....Number of elements with free volume charges
            % 
            % INPUT:
            % global_node_coordinates.....N_nodes x 2 matrix containing the
            %   x-coordinates of the globla nodes in the first column and the
            %   y-coordinates in the second column.
            % node_mapping.....N_e x N matrix containing the 
            % material_properties.....N_e x 2 matrix contianing the material
            %   properties of each element. The first column contains the material
            %   property in x-direction, the second column for the y-direction.
            % direchlet_boundary_conditions.....N_dir x 2 matrix containing the
            %   direchlet boundary conditions of the problem. Each row stands for one
            %   conditions, where the first column denotes the number of the global
            %   node and the second column the actual value. 
            % neumann_boundary_conditions.....N_neu x 3 matrix containing the neumann
            %   boundary conditions of the problem. Each row stands for one condition
            %   where the first column denotes the number of the finite element with
            %   the condition, the second column the side of the triangle and the
            %   third column the actual value of the condition.
            %   Multiple entries are used if one element contains multiple boundary
            %   conditions. (e.g. the two bottom sidef of a second order triangle)
            % colume_charges.....N_roh x 2 matrix containing information about the
            %   free volume charges within the problem area. Each row stands for one
            %   finite element with free volume charges. The first column contains
            %   the number of the finite element and the second column contains the
            %   CONSTANT charge density within the finite element.
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
            
            % Get number of elements (N_e) and number of element nodes (N)
            [N_e, N] = size(node_mapping);
            
            % Get number of global nodes
            [N_glob_nodes, ~] = size(global_node_coordinates);
            
            % Get number of direchlet boundary conditions
            [N_dir, ~] = size(direchlet_boundary_conditions);
            
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
                
                
                material_x = material_properties(e, 1);
                material_y = material_properties(e, 2);
                
                % Local coefficient matrix and right-side vector
                A_loc = zeros(N, N);
                r_loc = zeros(N, 1);
                
                row_idx = 1 : N;
                row_idx = row_idx(:);
                
                col_idx = 1 : N;
                col_idx = col_idx(:);
                
                % Vector containing logical ones at direchlet boundary conditions of 
                % current element.
                % Dimenison : N_dir (Number of direchlet boudnary conditions) x 1
                direchlet_conditions_of_currnet_element_bin_idx = ...
                    ismember(direchlet_boundary_conditions(:,1), local_node_mapping);
                
                % Vector containing logical ones at local nodes with direchlet 
                % boundary conditions.
                % Dimenstion N (Number of element nodes) x 1
                local_nodes_with_direchlet_conditions_bin_idx = ...
                    ismember(local_node_mapping, direchlet_boundary_conditions(:,1));
                
                % First column: Local nodes with direchlet boundary conditions
                % Second and third column: Global node + boundary conditions 
                % corresponding to local node in first column.
                %
                % E.g.:
                %
                % 1 25 10
                % 3 73 0
                %
                % Local node 1 corresponds to global node 25 with direchlet boudnary 
                % condition of 10 (V, A, etc.). Local node 3 corresponds to global 
                % node 73 with a boundary condition of 0.
                %
                direchlet_boundary_data = [...
                    row_idx(local_nodes_with_direchlet_conditions_bin_idx), ...
                    direchlet_boundary_conditions(...
                    direchlet_conditions_of_currnet_element_bin_idx, :)];
                
                
                % Get Neumann boundary conditions of element. It is a nx3 matrix,
                % where each row stands for one neumann boudnary condition of the
                % current element. The rows are of the following form:
                %
                % Column 1: Side of the current element with a boundary condition.
                % e.g. for a first order triangle, the side
                % between local nodes 1 and 3 is side 1, between local nodes 2 and 3
                % side 2 etc.
                %
                % Column 2: Boundary condition value.
                neumann_boundary_conditions_of_element = ...
                    neumann_boundary_conditions(neumann_boundary_conditions(:, 1) ...
                    == e, 2:3);
                
                
                % A direchlet boundary condition at local node k removes the k-th 
                % row(the k-th equation) of the element-equation system. The elements 
                % of A_loc will be zero in this row.
                row_idx(local_nodes_with_direchlet_conditions_bin_idx) = [];
                
                % Due to line above, the algorithm will only iterate above the 
                % elements without direchlet boundary conditions.
                for l = 1 : length(row_idx)
                    
                    row = row_idx(l);
                    
                    % ===== Calculation of right-side vector elements
                    
                    % === Plane integral
                    % Only needed for electrostatic problems, not for static current
                    % problemes
                    %
                    % Integration also not needed if no free charges occurr in the
                    % current element
                    if strcmpi(problem_type, 'E-Static') && volume_charges(e) ~= 0
                        fun = @(zeta, eta) ElectrostaticProblem. ...
                            right_side_plane_integrant(finite_element_type, xe, ...
                            ye, row, zeta, eta, volume_charges(e));
                        
                        r_loc(row) = GaussianQuadrature. ...
                            integrate_2D_normalized_triangle_region(fun,...
                            plane_integration_order);
                    end
                    
                    
                    % === Neumann integral
                    
                    [number_of_neumann_boundary_conditions, ~] = ...
                        size(neumann_boundary_conditions_of_element);
                    
                    for b = 1: number_of_neumann_boundary_conditions
                        side = neumann_boundary_conditions_of_element(b, 1);
                        value = neumann_boundary_conditions_of_element(b, 2);
                        
                        fun = @(t) ElectrostaticProblem. ...
                            right_side_neumann_integrant(finite_element_type, ...
                            xe, ye, row, t, side, value);
                        
                        r_loc(row)  = GaussianQuadrature.integrate_1d_line(fun, ...
                            curve_integration_order, 1, 0) + r_loc(row);
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
                        
                        % If the k-th local node contains a direchlet boundary 
                        % condition, the coeffient in the k-th column multiplied with 
                        % the value of the boundary condition is brought to the right 
                        % side of the current equation
                        if ismember(col, direchlet_boundary_data(:,1))
                            idx = (direchlet_boundary_data(:, 1) == col);
                            r_loc(row) = r_loc(row) - coef * ...
                                direchlet_boundary_data(idx, 3);
                        else
                            A_loc(row, col) = coef;
                        end
                    end
                end
                
                % ===== Add current local equation system to global equation 
                %   system =====
                
                global_node_to_equation_system_index_mapping = ...
                    Solver.get_node_index_to_equation_index_mapping( ...
                    N_glob_nodes, direchlet_boundary_conditions(:,1));
                
                % If the k-th local node has a direchlet boundary condition, the 
                % elements from the k-th column should not be included in the global 
                % equation system as they are already condidered in the right-side 
                % vector.
                col_idx(local_nodes_with_direchlet_conditions_bin_idx) = [];
                
                % global_node_to_equadion_system_mapping maps the global nodes to 
                % their corresponding row and column within the global equation 
                % system.
                % E.g. the unknown potential of global node 6 maps to the 5th 
                % row/column of the global equation system as, for example, the 3rd 
                % global node has a direchlet boundary condition.
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