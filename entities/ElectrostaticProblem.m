classdef ElectrostaticProblem
    
    methods(Static)
        % ------------ Functions for element matrix calculation
        function val = element_matrix_integrant(finite_element, k , l, ...
                zeta, eta, epsilon_x, epsilon_y)
            
            assert(isa(finite_element, 'SecondOrderTriangleElement'))
            
            % Get vairables needed for calculation
            jacobi_determinant = finite_element.jacobi_determinant( ...
                zeta, eta);
            
            dN_dZeta = ...
                finite_element.get_zeta_shape_function_derivatives(...
                zeta, eta);
            
            dN_dEta = ...
                finite_element.get_eta_shape_function_derivatives(...
                zeta, eta);
            
            xe = finite_element.xe;
            ye = finite_element.ye;
            
            
            % Split calculation into two parts and calculate each part
            % separately
            part_1 = epsilon_x / jacobi_determinant * ...
                (ye' * dN_dEta * dN_dZeta(k) - ye' * dN_dZeta * dN_dEta(k)) ...
                * ...
                (ye' * dN_dEta * dN_dZeta(l) - ye' * dN_dZeta * dN_dEta(l));
            
            
            part_2 = epsilon_y / jacobi_determinant * ...
                (- xe' * dN_dEta * dN_dZeta(k) + xe' * dN_dZeta * dN_dEta(k)) ...
                * ...
                (-xe' * dN_dEta * dN_dZeta(l) + xe' * dN_dZeta * dN_dEta(l));
            
            
            val = part_1 + part_2;
        end
        
    end
    
end