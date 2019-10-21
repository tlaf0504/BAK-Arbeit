classdef RotationSymmetricMagnetostaticProblem
    
    methods(Static)
        % ------------ Functions for element matrix calculation
        function val = element_matrix_integrant(finite_element_class, re, ze, k , l, ...
                zeta, eta, mu_r, mu_z)
            
            if ~iscolumn(re)
                re = re';
            end
            
            if ~iscolumn(ze)
                ze = ze';
            end
            
            
            % Get vairables needed for calculation
            jacobi_determinant = finite_element_class.jacobi_determinant( ...
                zeta, eta, re, ze)';
            
            dN_dZeta = finite_element_class. ...
                get_zeta_shape_function_derivative_matrix(zeta, eta);
            
            dN_dEta = ...
                finite_element_class. ...
                get_eta_shape_function_derivative_matrix(zeta, eta);
            
            % Shape functions of element nodes k and l
            N_k = finite_element_class.get_shape_function(k, zeta, eta);
            N_l = finite_element_class.get_shape_function(l, zeta, eta);
            
            % Partial derivatives of N_k and N_l after r and z each
            dNk_dr = (1 ./jacobi_determinant)' .* ...
                ( ze' * dN_dEta .* dN_dZeta(k,:) - ...
                  ze' * dN_dZeta .* dN_dEta(k,:));
              
            dNk_dz = (1 ./jacobi_determinant)' .* ...
                ( - re' * dN_dEta .* dN_dZeta(k,:) + ...
                  re' * dN_dZeta .* dN_dEta(k,:));
              
              
            dNl_dr = (1 ./jacobi_determinant)' .* ...
                ( ze' * dN_dEta .* dN_dZeta(l,:) - ...
                  ze' * dN_dZeta .* dN_dEta(l,:));
              
            dNl_dz = (1 ./jacobi_determinant)' .* ...
                ( - re' * dN_dEta .* dN_dZeta(l,:) + ...
                  re' * dN_dZeta .* dN_dEta(l,:));
              
            
            r = re' * finite_element_class.get_shape_function_matrix(zeta, eta);
            
            
            val = (mu_z * ((N_k.*N_l)./r + N_k.*dNk_dr + N_l.*dNl_dr + r.*dNk_dr.*dNl_dr) + ...
                  mu_r * r .* dNk_dz .* dNl_dz) .* jacobi_determinant';
        end
        
        function val = right_side_plane_integrant(finite_element_class, re, ze, k, zeta, ...
                eta, J)
            
            if ~iscolumn(re)
                re = re';
            end
            
            if ~iscolumn(ze)
                ze = ze';
            end
            
            r = re' * finite_element_class.get_shape_function_matrix(zeta, eta);
            
            shape_function = finite_element_class.get_shape_function(k, zeta, ...
                eta);
            
            val = shape_function .* J .* r .*...
                finite_element_class.jacobi_determinant(zeta, eta, re, ze);
        end
        
        
        function val = right_side_neumann_integrant(finite_element_class, re, ze, k, ...
                t, side, alpha)
            
            if ~iscolumn(re)
                re = re';
            end
            
            if ~iscolumn(ze)
                ze = ze';
            end
            
            

            zero_vec = zeros(1,length(t));
            
            if side == 1
                
                r = re' * finite_element_class.get_shape_function_matrix(t, zero_vec);
                
                dN_dZeta = finite_element_class. ...
                    get_zeta_shape_function_derivative_matrix(t, zero_vec);
                
                val = finite_element_class.get_shape_function(k, t, zero_vec) .* ...
                    alpha .* r .*...
                    sqrt((re' * dN_dZeta).^2 + (ze' * dN_dZeta).^2);
                
            elseif side == 2
                
                r = re' * finite_element_class.get_shape_function_matrix(zero_vec, t);
                
                dN_dEta = finite_element_class. ...
                    get_eta_shape_function_derivative_matrix(zero_vec, t);
                
                val = finite_element_class.get_shape_function(k, zero_vec, t).* ...
                    alpha .* r .*...
                    sqrt((re' * dN_dEta).^2 + (ze' * dN_dEta).^2);
                
            elseif side == 3
                zeta = 1-t;
                eta = t;
                
                r = re' * finite_element_class.get_shape_function_matrix(zeta, eta);
                
                dN_dZeta = finite_element_class. ...
                    get_zeta_shape_function_derivative_matrix(zeta, eta);
                
                dN_dEta = finite_element_class. ...
                    get_eta_shape_function_derivative_matrix(zeta, eta);
                
                val = finite_element_class.get_shape_function(k, zeta, eta).* ...
                    alpha .* r .*...
                    sqrt((re' * dN_dEta - re' * dN_dZeta).^2 + ...
                    (ze' * dN_dEta - ze' * dN_dZeta).^2);
                
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