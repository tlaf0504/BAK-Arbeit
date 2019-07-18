classdef FirstOrderTriangleElement
    %FIRSTORDERTRIANGLEELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        nodes_of_triangle_edges = [1,1,0;
                                   0,1,1;
                                   1,0,1];
                               
        number_of_nodes = 3;
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
                val = FirstOrderTriangleElement.shape_function_1( ...
                    zeta, eta);
                
            elseif k == 2
                val = FirstOrderTriangleElement.shape_function_2( ...
                    zeta, eta);
                
                
            elseif k == 3
                val = FirstOrderTriangleElement.shape_function_3( ...
                    zeta, eta);
            end
            
        end
        
        function det_J = jacobi_determinant(eta, zeta, xe, ye)
            xe = xe(:);
            ye = ye(:);
            
            number_of_required_output_values = length(eta);
            
            % For linear elements, the jacobi determinant is independent from zeta
            % and eta.
            % Output is row-vector
            det_J = ((-xe(1) + xe(2)) * (-ye(1) + ye(3)) - ...
                (-ye(1) + ye(2)) * (-xe(1) + xe(3))) + ...
                ones(1, number_of_required_output_values);     
        end
        
        function J = jacobi_matrix(zeta, eta, xe, ye)
            j11 = (ye' * FirstOrderTriangleElement.v21) * zeta + ...
                (ye' * FirstOrderTriangleElement.v22) * eta + ...
                (ye' * FirstOrderTriangleElement.v23);
           
           j12 = -((ye' * FirstOrderTriangleElement.v11) * zeta + ...
               (ye' * FirstOrderTriangleElement.v12) * eta + ...
               (ye' * FirstOrderTriangleElement.v13) );
           
           j21 = -( (xe' * FirstOrderTriangleElement.v21) * zeta + ...
               (xe' * FirstOrderTriangleElement.v22) * eta + ...
               (xe' * FirstOrderTriangleElement.v23) );
           
           j22 = (xe' * FirstOrderTriangleElement.v11) * zeta + ...
               (xe' * FirstOrderTriangleElement.v12) * eta + ...
               (xe' * FirstOrderTriangleElement.v13);
           
           J = [j11, j12; j21, j22];
        end
        
        function J_inv = inverse_jacobi_matrix(zeta, eta, xe, ye)
           
           j11 = (-xe(1) + xe(2));
           
           j12 = (-ye(1) + ye(2));
           
           j21 = (-xe(1) + xe(3));
           
           j22 = (-ye(1) + ye(3));
            
           J_inv = (1 ./ FirstOrderTriangleElement.jacobi_determinant(zeta, eta, xe, ye)) .* ...
               [j22, -j12; -j21, j11];
        end
        
        function val = get_zeta_shape_function_derivative(k , zeta, ...
                eta)
            
            assert(k > 0 & k <= FristOrderTriangleElement.number_of_nodes)
            assert(all(zeta >= 0 & zeta <= 1 & eta >= 0 & eta <= 1))
            
            if k == 1
                val = FristOrderTriangleElement.dN1_dZeta(zeta, eta);
                
            elseif k == 2
                val = FristOrderTriangleElement.dN2_dZeta(zeta, eta);
                
            elseif k == 3
                val = FristOrderTriangleElement.dN3_dZeta(zeta, eta);
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
            mat = [FirstOrderTriangleElement.dN1_dZeta(zeta, eta);
                   FirstOrderTriangleElement.dN2_dZeta(zeta, eta);
                   FirstOrderTriangleElement.dN3_dZeta(zeta, eta);
                   ];
            
        end
        
        function val = get_eta_shape_function_derivative(k , zeta, ...
                eta)
            
            assert(k > 0 & k <= FristOrderTriangleElement.number_of_nodes)
            assert(all(zeta >= 0 & zeta <= 1 & eta >= 0 & eta <= 1))
            
            if k == 1
                val = FristOrderTriangleElement.dN1_dEta(zeta, eta);
                
            elseif k == 2
                val = FristOrderTriangleElement.dN2_dEta(zeta, eta);
                
            elseif k == 3
                val = FristOrderTriangleElement.dN3_dEta(zeta, eta);
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
            mat = [FirstOrderTriangleElement.dN1_dEta(zeta, eta);
                   FirstOrderTriangleElement.dN2_dEta(zeta, eta);
                   FirstOrderTriangleElement.dN3_dEta(zeta, eta);
                   ];
            
        end
        
        
        function mat = get_shape_function_matrix(zeta, eta)
            
             mat = [FirstOrderTriangleElement.shape_function_1(zeta, eta);
                   FirstOrderTriangleElement.shape_function_2(zeta, eta);
                   FirstOrderTriangleElement.shape_function_3(zeta, eta);
                    ];
        end
        
    end
end

