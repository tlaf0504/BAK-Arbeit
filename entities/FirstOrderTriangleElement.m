classdef FirstOrderTriangleElement
    %FIRSTORDERTRIANGLEELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        nodes_of_triangle_edges = [1,1,0;
                                   1,0,1;
                                   0,1,1];
    end
    
    methods(Static)
        % --------- Element shape functions and their derivatives
        
        % ----- Shape functions
        % Shape function of node 1
        function val = shape_function_1(zeta, eta)
           val = 1 - eta - zeta;
        end
        
        % Shape function of node 2
        function val = shape_function_2(zeta, eta)
            val = zeta;
        end
        
        % Shape function of node 3
        function val = shape_function_3(zeta, eta)
            val = eta;
        end
        
        
        % ----- Derivatives of shape functions
        
        % --- Zeta derivatives
        function val = dN1_dZeta(zeta, eta)
            [a,b] = size(zeta);
            val = -ones(a,b);
        end
        
        function val = dN2_dZeta(zeta, eta)
            [a,b] = size(zeta);
            val = ones(a,b);
        end 
        
        function val = dN3_dZeta(zeta, eta)
            [a,b] = size(zeta);
            val = zeros(a,b);
        end 
        
        % --- Eta derivatives
        function val = dN1_dEta(zeta, eta)
            [a,b] = size(zeta);
            val = -ones(a,b);
        end
        
        function val = dN2_dEta(zeta, eta)
            [a,b] = size(zeta);
            val = -zeros(a,b);
        end 
        
        function val = dN3_dEta(zeta, eta)
            [a,b] = size(zeta);
            val = ones(a,b);
        end 
        
        
        % ------- Misc functions
        
        function val = get_shape_function(k, zeta, eta)
            
            assert(k > 0 & k <= FirstOrderTriangleElement.number_of_nodes)
            assert(all(zeta >= 0 & zeta <= 1 & eta >= 0 & eta <= 1))
            
            if k == 1
                val = SecondOrderTriangleElement.shape_function_1( ...
                    zeta, eta);
                
            elseif k == 2
                val = SecondOrderTriangleElement.shape_function_2( ...
                    zeta, eta);
                
                
            elseif k == 3
                val = SecondOrderTriangleElement.shape_function_3( ...
                    zeta, eta);
            end
            
        end
        
        function det_J = jacobi_determinant(eta, zeta, xe, ye)
            if ~iscolumn(xe)
                xe = xe';
            end
            
            if ~iscolumn(ye)
                ye = ye';
            end
            
            
            
            det_J = 0;
            
            
        end
        
        function J = jacobi_matrix(zeta, eta, xe, ye)
            j11 = (ye' * SecondOrderTriangleElement.v21) * zeta + ...
                (ye' * SecondOrderTriangleElement.v22) * eta + ...
                (ye' * SecondOrderTriangleElement.v23);
           
           j12 = -((ye' * SecondOrderTriangleElement.v11) * zeta + ...
               (ye' * SecondOrderTriangleElement.v12) * eta + ...
               (ye' * SecondOrderTriangleElement.v13) );
           
           j21 = -( (xe' * SecondOrderTriangleElement.v21) * zeta + ...
               (xe' * SecondOrderTriangleElement.v22) * eta + ...
               (xe' * SecondOrderTriangleElement.v23) );
           
           j22 = (xe' * SecondOrderTriangleElement.v11) * zeta + ...
               (xe' * SecondOrderTriangleElement.v12) * eta + ...
               (xe' * SecondOrderTriangleElement.v13);
           
           J = [j11, j12; j21, j22];
        end
        
        function J_inv = inverse_jacobi_matrix(zeta, eta, xe, ye)
           
           j11 = (ye' * SecondOrderTriangleElement.v21) * zeta + ...
                (ye' * SecondOrderTriangleElement.v22) * eta + ...
                (ye' * SecondOrderTriangleElement.v23);
           
           j12 = -((ye' * SecondOrderTriangleElement.v11) * zeta + ...
               (ye' * SecondOrderTriangleElement.v12) * eta + ...
               (ye' * SecondOrderTriangleElement.v13) );
           
           j21 = -( (xe' * SecondOrderTriangleElement.v21) * zeta + ...
               (xe' * SecondOrderTriangleElement.v22) * eta + ...
               (xe' * SecondOrderTriangleElement.v23) );
           
           j22 = (xe' * SecondOrderTriangleElement.v11) * zeta + ...
               (xe' * SecondOrderTriangleElement.v12) * eta + ...
               (xe' * SecondOrderTriangleElement.v13);
            
           J_inv = (1 ./ obj.jacobi_determinant(zeta, eta, xe, ye)) .* ...
               [j22, -j12; -j21, j11];
        end
        
        function val = get_zeta_shape_function_derivative(k , zeta, ...
                eta)
            
            assert(k > 0 & k <= SecondOrderTriangleElement.number_of_nodes)
            assert(all(zeta >= 0 & zeta <= 1 & eta >= 0 & eta <= 1))
            
            if k == 1
                val = SecondOrderTriangleElement.dN1_dZeta(zeta, eta);
                
            elseif k == 2
                val = SecondOrderTriangleElement.dN2_dZeta(zeta, eta);
                
            elseif k == 3
                val = SecondOrderTriangleElement.dN3_dZeta(zeta, eta);
                
            elseif k == 4
                val = SecondOrderTriangleElement.dN4_dZeta(zeta, eta);
                
            elseif k == 5
                val = SecondOrderTriangleElement.dN5_dZeta(zeta, eta);
                
            elseif k == 6
                val = SecondOrderTriangleElement.dN6_dZeta(zeta, eta);
            end
        end
     
        function mat = get_zeta_shape_function_derivative_matrix(...
                zeta, eta)
            
            assert(all(zeta >= 0 & zeta <= 1 & eta >= 0 & eta <= 1))
            
            zeta = zeta(:); zeta = zeta'; % Ensure that zeta is a row vector
            eta = eta(:); eta = eta'; % Ensure that eta is a row vector
            
            % mat contains derivates for one zeta/eta pair in wach column
            % mat = 
            % [dN1_dEta(zeta(1), eta(1)), dN1_dEta(zeta(2), eta(2)), ...
            %  dN2_dEta(zeta(1), eta(1)), dN2_dEta(zeta(2), eta(2)), ...
            %  .
            %  .
            %  .
            % ]
            mat = [SecondOrderTriangleElement.dN1_dZeta(zeta, eta);
                   SecondOrderTriangleElement.dN2_dZeta(zeta, eta);
                   SecondOrderTriangleElement.dN3_dZeta(zeta, eta);
                   ];
            
        end
        
        function val = get_eta_shape_function_derivative(k , zeta, ...
                eta)
            
            assert(k > 0 & k <= SecondOrderTriangleElement.number_of_nodes)
            assert(all(zeta >= 0 & zeta <= 1 & eta >= 0 & eta <= 1))
            
            if k == 1
                val = SecondOrderTriangleElement.dN1_dEta(zeta, eta);
                
            elseif k == 2
                val = SecondOrderTriangleElement.dN2_dEta(zeta, eta);
                
            elseif k == 3
                val = SecondOrderTriangleElement.dN3_dEta(zeta, eta);
                
            elseif k == 4
                val = SecondOrderTriangleElement.dN4_dEta(zeta, eta);
                
            elseif k == 5
                val = SecondOrderTriangleElement.dN5_dEta(zeta, eta);
                
            elseif k == 6
                val = SecondOrderTriangleElement.dN6_dEta(zeta, eta);
            end
        end
        
        
        function mat = get_eta_shape_function_derivative_matrix( ...
                zeta, eta)
            
            assert(all(zeta >= 0 & zeta <= 1 & eta >= 0 & eta <= 1))
            
            zeta = zeta(:); zeta = zeta'; % Ensure that zeta is a row vector
            eta = eta(:); eta = eta'; % Ensure that eta is a row vector
            
            % mat contains derivates for one zeta/eta pair in wach column
            % mat = 
            % [dN1_dEta(zeta(1), eta(1)), dN1_dEta(zeta(2), eta(2)), ...
            %  dN2_dEta(zeta(1), eta(1)), dN2_dEta(zeta(2), eta(2)), ...
            %  .
            %  .
            %  .
            % ]
            mat = [SecondOrderTriangleElement.dN1_dEta(zeta, eta);
                   SecondOrderTriangleElement.dN2_dEta(zeta, eta);
                   SecondOrderTriangleElement.dN3_dEta(zeta, eta);
                   SecondOrderTriangleElement.dN4_dEta(zeta, eta);
                   SecondOrderTriangleElement.dN5_dEta(zeta, eta);
                   SecondOrderTriangleElement.dN6_dEta(zeta, eta)];
            
        end
        
    end
end

