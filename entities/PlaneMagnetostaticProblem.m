classdef PlaneMagnetostaticProblem
    
    methods(Static)
        % ------------ Functions for element matrix calculation
        function val = element_matrix_integrant(finite_element_class, xe, ye, k , l, ...
                zeta, eta, mu_x, mu_y)
            
            if ~iscolumn(xe)
                xe = xe';
            end
            
            if ~iscolumn(ye)
                ye = ye';
            end
            
            
            % Get vairables needed for calculation
            jacobi_determinant = finite_element_class.jacobi_determinant( ...
                zeta, eta, xe, ye)';
            
            dN_dZeta = finite_element_class. ...
                get_zeta_shape_function_derivative_matrix(zeta, eta);
            
            dN_dEta = ...
                finite_element_class. ...
                get_eta_shape_function_derivative_matrix(zeta, eta);
            
          
            % Partial derivatives of N_k and N_l after r and z each
            dNk_dx = (1 ./jacobi_determinant)' .* ...
                ( ye' * dN_dEta .* dN_dZeta(k,:) - ...
                  ye' * dN_dZeta .* dN_dEta(k,:));
              
            dNk_dy = (1 ./jacobi_determinant)' .* ...
                ( - xe' * dN_dEta .* dN_dZeta(k,:) + ...
                  xe' * dN_dZeta .* dN_dEta(k,:));
              
              
            dNl_dx = (1 ./jacobi_determinant)' .* ...
                ( ye' * dN_dEta .* dN_dZeta(l,:) - ...
                  ye' * dN_dZeta .* dN_dEta(l,:));
              
            dNl_dy = (1 ./jacobi_determinant)' .* ...
                ( - xe' * dN_dEta .* dN_dZeta(l,:) + ...
                  xe' * dN_dZeta .* dN_dEta(l,:));

            
            val = (mu_y * (dNk_dx .* dNl_dx) + ...
                  mu_x * (dNk_dy .* dNl_dy)) .* jacobi_determinant';
        end
        
        function val = right_side_plane_integrant(finite_element_class, xe, ye, k, zeta, ...
                eta, J)
            
            if ~iscolumn(xe)
                xe = xe';
            end
            
            if ~iscolumn(ye)
                ye = ye';
            end

            shape_function = finite_element_class.get_shape_function(k, zeta, ...
                eta);
            
            val = shape_function .* J .* ...
                finite_element_class.jacobi_determinant(zeta, eta, xe, ye);
        end
        
        
        function val = right_side_neumann_integrant(finite_element_class, xe, ye, k, ...
                t, side, alpha)
            
            if ~iscolumn(xe)
                xe = xe';
            end
            
            if ~iscolumn(ye)
                ye = ye';
            end
            
            

            zero_vec = zeros(1,length(t));
            
            if side == 1
                
                
                
                dN_dZeta = finite_element_class. ...
                    get_zeta_shape_function_derivative_matrix(t, zero_vec);
                
                val = finite_element_class.get_shape_function(k, t, zero_vec) .* ...
                    alpha .* ...
                    sqrt((xe' * dN_dZeta).^2 + (ye' * dN_dZeta).^2);
                
            elseif side == 2
                
                
                
                dN_dEta = finite_element_class. ...
                    get_eta_shape_function_derivative_matrix(zero_vec, t);
                
                val = finite_element_class.get_shape_function(k, zero_vec, t).* ...
                    alpha .* ...
                    sqrt((xe' * dN_dEta).^2 + (ye' * dN_dEta).^2);
                
            elseif side == 3
                zeta = 1-t;
                eta = t;
                
                
                
                dN_dZeta = finite_element_class. ...
                    get_zeta_shape_function_derivative_matrix(zeta, eta);
                
                dN_dEta = finite_element_class. ...
                    get_eta_shape_function_derivative_matrix(zeta, eta);
                
                val = finite_element_class.get_shape_function(k, zeta, eta).* ...
                    alpha .*...
                    sqrt((xe' * dN_dEta - xe' * dN_dZeta).^2 + ...
                    (ye' * dN_dEta - ye' * dN_dZeta).^2);
                
           end
                    
        end
        
        function val = get_energy_integrant(finite_element_class, xe, ye, ...
                node_potentials, zeta, eta, epsilon_x, epsilon_y)
            error('Not implemented!')
            
            if ~iscolumn(xe)
                xe = xe';
            end
            
            if ~iscolumn(ye)
                ye = ye';
            end
            
            if ~iscolumn(node_potentials)
                node_potentials = node_potentials';
            end
            
            
            % Get vairables needed for calculation
            jacobi_determinant = finite_element_class.jacobi_determinant( ...
                zeta, eta, xe, ye);
            
            dN_dZeta = finite_element_class. ...
                get_zeta_shape_function_derivative_matrix(zeta, eta);
            
            dN_dEta = ...
                finite_element_class. ...
                get_eta_shape_function_derivative_matrix(zeta, eta);
            
            part_1 = epsilon_x ./ jacobi_determinant .* ( ...
                (node_potentials' * dN_dZeta) .* (ye' * dN_dEta) - ...
                (node_potentials' * dN_dEta) .* (ye' * dN_dZeta)).^2;
            
            part_2 = epsilon_y ./ jacobi_determinant .* ( ...
                -(node_potentials' * dN_dZeta) .* (xe' * dN_dEta) + ...
                (node_potentials' * dN_dEta) .* (xe' * dN_dZeta)).^2;
            
            val = part_1 + part_2;
            
        end
    end
    
end