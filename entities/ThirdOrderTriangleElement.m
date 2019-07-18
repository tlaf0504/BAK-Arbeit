classdef ThirdOrderTriangleElement
    
    properties(Constant)
        number_of_nodes = 10;
        
        % Assignment of element nodes to the corresponding triangle edges.
        %
        % This assignment is needed for correct
        % calculation of the right side coefficients when the current triangle is
        % located on the neumann boundary, as the line integral is depended of the
        % triangles edge on the boundary.
        %
        % Triangle edge '1' is defined as the edge on the zeta-axis
        % Triangle edge '2' is defined as the edge on the eta-axis
        % Triangle edge '3' is defined as the edge between the points (1,0) and (0,1)
        %
        % Each row of the matrix below corresponds to one triangle edge.
        nodes_of_triangle_edges = [1,1,1,1,0,0,0,0,0,0;
                                   0,0,0,1,1,1,1,0,0,0;
                                   0,0,0,0,0,0,1,1,1,1];
                               
                               
       % Correction of node ordering of triangles and curves.
       % Gmsh has a different node ordering than used in this program. See 
       % http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering (15.07.2019) for
       % further information.
       % Attention: Gmsh starts node numbering with 0. As matlab starts indexing with
       % 1, the gmsh node numbers also start with 1.
       
                              % Gmsh  Program
       triangle_node_ordering = [ 1,      1;
                                  2,      4;
                                  3,      7;
                                  4,      2;
                                  5,      3;
                                  6,      5;
                                  7,      6;
                                  8,      8;
                                  9,      9;
                                  10,    10];
                              
                              
                        % Gmsh  Program 
      curve_node_ordering = [1,    1;
                             2,    4;
                             3,    2;
                             4,    3];
    end
    
    
    methods(Static)
        
        % --------- Element shape functions and their derivatives
        
        % ----- Shape functions
        % Shape function of node 1
        function val = shape_function_1(zeta, eta)
           val = - (9*eta.^3)/2 - (27*eta.^2.*zeta)/2 + 9*eta.^2 - (27*eta.*zeta.^2)/2 + ...
               18*eta.*zeta - (11*eta)/2 - (9*zeta.^3)/2 + 9*zeta.^2 - (11*zeta)/2 + 1;
        end
        
        % Shape function of node 2
        function val = shape_function_2(zeta, eta)
            val = (27*eta.^2.*zeta)/2 + 27*eta.*zeta.^2 - (45*eta.*zeta)/2 + ...
                (27*zeta.^3)/2 - (45*zeta.^2)/2 + 9*zeta;
        end
        
        % Shape function of node 3
        function val = shape_function_3(zeta, eta)
            val = (9*eta.*zeta)/2 - (9*zeta)/2 - (27*eta.*zeta.^2)/2 + 18*zeta.^2 - ...
                (27*zeta.^3)/2;
        end
        
        % Shape function of node 4
        function val = shape_function_4(zeta, eta)
            val = (9*zeta.^3)/2 - (9*zeta.^2)/2 + zeta;
        end
        
        % Shape function of node 5
        function val = shape_function_5(zeta, eta)
            val = (27*eta.*zeta.^2)/2 - (9*eta.*zeta)/2;
        end
        
        % Shape function of node 6
        function val = shape_function_6(zeta, eta)
            val = (27*zeta.*eta.^2)/2 - (9*zeta.*eta)/2;
        end
        
        
        % Shape function of node 7
        function val = shape_function_7(zeta, eta)
            val = (9*eta.^3)/2 - (9*eta.^2)/2 + eta;
        end
        
        
        % Shape function of node 8
        function val = shape_function_8(zeta, eta)
            val = (9*eta.*zeta)/2 - (9*eta)/2 - (27*eta.^2.*zeta)/2 + 18*eta.^2 - ...
                (27*eta.^3)/2;
        end
        
        
        % Shape function of node 9
        function val = shape_function_9(zeta, eta)
            val = (27*eta.^3)/2 + 27*eta.^2.*zeta - (45*eta.^2)/2 + (27*eta.*zeta.^2)/2 - ...
                (45*eta.*zeta)/2 + 9*eta;
        end
        
        
        % Shape function of node 10
        function val = shape_function_10(zeta, eta)
            val =  - 27*eta.^2.*zeta - 27*eta.*zeta.^2 + 27*eta.*zeta;
        end
        
        
        % ----- Derivatives of shape functions
        
        % --- Zeta derivatives
        function val = dN1_dZeta(zeta, eta)
            val = 18*eta + 18*zeta - 27*eta.*zeta - (27*eta.^2)/2 - (27*zeta.^2)/2 - ...
                11/2;
        end
        
        function val = dN2_dZeta(zeta, eta)
            val = (27*eta.^2)/2 + 54*eta.*zeta - (45*eta)/2 + (81*zeta.^2)/2 - ...
                45*zeta + 9;
        end 
        
        function val = dN3_dZeta(zeta, eta)
            val = (9*eta)/2 + 36*zeta - 27*eta.*zeta - (81*zeta.^2)/2 - 9/2;
        end 
        
        function val = dN4_dZeta(zeta, eta)
            val = (27*zeta.^2)/2 - 9*zeta + 1;
        end 
        
        function val = dN5_dZeta(zeta, eta)
            val = 27*eta.*zeta - (9*eta)/2;
        end 
        
        function val = dN6_dZeta(zeta, eta)
            val = (27*eta.^2)/2 - (9*eta)/2;
        end
        
        function val = dN7_dZeta(zeta, eta)
            [a,b] = size(zeta);
            val = zeros(a,b);
        end
        
        function val = dN8_dZeta(zeta, eta)
            val = (9*eta)/2 - (27*eta.^2)/2;
        end
        
        function val = dN9_dZeta(zeta, eta)
            val = 27*eta.*zeta - (45*eta)/2 + 27*eta.^2;
        end
        
        function val = dN10_dZeta(zeta, eta)
            val = 27*eta - 54*eta.*zeta - 27*eta.^2;
        end
        
        % --- Eta derivatives
        function val = dN1_dEta(zeta, eta)
            val = 18*eta + 18*zeta - 27*eta.*zeta - (27*eta.^2)/2 - (27*zeta.^2)/2 - ...
                11/2;
        end
        
        function val = dN2_dEta(zeta, eta)
            val = 27*eta.*zeta - (45*zeta)/2 + 27*zeta.^2;
        end 
        
        function val = dN3_dEta(zeta, eta)
            val = (9*zeta)/2 - (27*zeta.^2)/2;
        end 
        
        function val = dN4_dEta(zeta, eta)
            [a,b] = size(zeta);
            val = zeros(a,b);
        end 
        
        function val = dN5_dEta(zeta, eta)
            val = (27*zeta.^2)/2 - (9*zeta)/2;
        end 
        
        function val = dN6_dEta(zeta, eta)
            val = 27*eta.*zeta - (9*zeta)/2;
        end
        
        function val = dN7_dEta(zeta, eta)
            val = (27*eta.^2)/2 - 9*eta + 1;
        end
        
        function val = dN8_dEta(zeta, eta)
            val = 36*eta + (9*zeta)/2 - 27*eta.*zeta - (81*eta.^2)/2 - 9/2;
        end
        
        function val = dN9_dEta(zeta, eta)
            val = (81*eta.^2)/2 + 54*eta.*zeta - 45*eta + (27*zeta.^2)/2 - ...
                (45*zeta)/2 + 9;
        end
        
        function val = dN10_dEta(zeta, eta)
            val = 27*zeta - 54*eta.*zeta - 27*zeta.^2;
        end
        
        
        % ------- Misc functions
        
        function val = get_shape_function(k, zeta, eta)
            
            assert(k > 0 & k <= ThirdOrderTriangleElement.number_of_nodes)
            assert(all(zeta >= 0 & zeta <= 1 & eta >= 0 & eta <= 1))
            
            if k == 1
                val = ThirdOrderTriangleElement.shape_function_1( ...
                    zeta, eta);
                
            elseif k == 2
                val = ThirdOrderTriangleElement.shape_function_2( ...
                    zeta, eta);
                
                
            elseif k == 3
                val = ThirdOrderTriangleElement.shape_function_3( ...
                    zeta, eta);
                
            elseif k == 4
                val = ThirdOrderTriangleElement.shape_function_4( ...
                    zeta, eta);
                
            elseif k == 5
                val = ThirdOrderTriangleElement.shape_function_5( ...
                    zeta, eta);
                
            elseif k == 6
                val = ThirdOrderTriangleElement.shape_function_6( ...
                    zeta, eta);
                
            elseif k == 7
                val = ThirdOrderTriangleElement.shape_function_7( ...
                    zeta, eta);
                
            elseif k == 8
                val = ThirdOrderTriangleElement.shape_function_8( ...
                    zeta, eta);
                
            elseif k == 9
                val = ThirdOrderTriangleElement.shape_function_9( ...
                    zeta, eta);
                
            elseif k == 10
                val = ThirdOrderTriangleElement.shape_function_10( ...
                    zeta, eta);
            end
            
            
            
        end
        
        function val = get_zeta_shape_function_derivative(k , zeta, ...
                eta)
            
            assert(k > 0 & k <= ThirdOrderTriangleElement.number_of_nodes)
            assert(all(zeta >= 0 & zeta <= 1 & eta >= 0 & eta <= 1))
            
            if k == 1
                val = ThirdOrderTriangleElement.dN1_dZeta(zeta, eta);
                
            elseif k == 2
                val = ThirdOrderTriangleElement.dN2_dZeta(zeta, eta);
                
            elseif k == 3
                val = ThirdOrderTriangleElement.dN3_dZeta(zeta, eta);
                
            elseif k == 4
                val = ThirdOrderTriangleElement.dN4_dZeta(zeta, eta);
                
            elseif k == 5
                val = ThirdOrderTriangleElement.dN5_dZeta(zeta, eta);
                
            elseif k == 6
                val = ThirdOrderTriangleElement.dN6_dZeta(zeta, eta);
                
            elseif k == 7
                val = ThirdOrderTriangleElement.dN7_dZeta(zeta, eta);
                
            elseif k == 8
                val = ThirdOrderTriangleElement.dN8_dZeta(zeta, eta);
                
            elseif k == 9
                val = ThirdOrderTriangleElement.dN9_dZeta(zeta, eta);
                
            elseif k == 10
                val = ThirdOrderTriangleElement.dN10_dZeta(zeta, eta);
                
            end
        end
        
        function val = get_eta_shape_function_derivative(k , zeta, ...
                eta)
            
            assert(k > 0 & k <= ThirdOrderTriangleElement.number_of_nodes)
            assert(all(zeta >= 0 & zeta <= 1 & eta >= 0 & eta <= 1))
            
            if k == 1
                val = ThirdOrderTriangleElement.dN1_dEta(zeta, eta);
                
            elseif k == 2
                val = ThirdOrderTriangleElement.dN2_dEta(zeta, eta);
                
            elseif k == 3
                val = ThirdOrderTriangleElement.dN3_dEta(zeta, eta);
                
            elseif k == 4
                val = ThirdOrderTriangleElement.dN4_dEta(zeta, eta);
                
            elseif k == 5
                val = ThirdOrderTriangleElement.dN5_dEta(zeta, eta);
                
            elseif k == 6
                val = ThirdOrderTriangleElement.dN6_dEta(zeta, eta);
                
            elseif k == 7
                val = ThirdOrderTriangleElement.dN7_dEta(zeta, eta);
                
            elseif k == 8
                val = ThirdOrderTriangleElement.dN8_dEta(zeta, eta);
                
            elseif k == 9
                val = ThirdOrderTriangleElement.dN9_dEta(zeta, eta);
                
            elseif k == 10
                val = ThirdOrderTriangleElement.dN10_dEta(zeta, eta);
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
            mat = [ThirdOrderTriangleElement.dN1_dZeta(zeta, eta);
                   ThirdOrderTriangleElement.dN2_dZeta(zeta, eta);
                   ThirdOrderTriangleElement.dN3_dZeta(zeta, eta);
                   ThirdOrderTriangleElement.dN4_dZeta(zeta, eta);
                   ThirdOrderTriangleElement.dN5_dZeta(zeta, eta);
                   ThirdOrderTriangleElement.dN6_dZeta(zeta, eta);
                   ThirdOrderTriangleElement.dN7_dZeta(zeta, eta);
                   ThirdOrderTriangleElement.dN8_dZeta(zeta, eta);
                   ThirdOrderTriangleElement.dN9_dZeta(zeta, eta);
                   ThirdOrderTriangleElement.dN10_dZeta(zeta, eta)];
            
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
            mat = [ThirdOrderTriangleElement.dN1_dEta(zeta, eta);
                   ThirdOrderTriangleElement.dN2_dEta(zeta, eta);
                   ThirdOrderTriangleElement.dN3_dEta(zeta, eta);
                   ThirdOrderTriangleElement.dN4_dEta(zeta, eta);
                   ThirdOrderTriangleElement.dN5_dEta(zeta, eta);
                   ThirdOrderTriangleElement.dN6_dEta(zeta, eta);
                   ThirdOrderTriangleElement.dN7_dEta(zeta, eta);
                   ThirdOrderTriangleElement.dN8_dEta(zeta, eta);
                   ThirdOrderTriangleElement.dN9_dEta(zeta, eta);
                   ThirdOrderTriangleElement.dN10_dEta(zeta, eta)];
            
        end
        
        function [j11, j12, j21, j22] = jacobi_matrix(zeta, eta, xe, ye)
            
            xe = xe(:); xe = xe';
            ye = ye(:); ye = ye';
            
            zeta = zeta(:);
            eta = eta(:);
            
            zeta_derivatives = ...
                ThirdOrderTriangleElement.get_zeta_shape_function_derivative_matrix(...
                zeta, eta);
            
            eta_derivatives = ...
                ThirdOrderTriangleElement.get_eta_shape_function_derivative_matrix(...
                zeta, eta);
            
            
            j11 = xe * zeta_derivatives;
            j12 = ye * zeta_derivatives;
            j21 = xe * eta_derivatives;
            j22 = ye * eta_derivatives;
            
        end
        
        function det_J = jacobi_determinant(zeta, eta, xe, ye)
            [j11, j12, j21, j22] = ...
                ThirdOrderTriangleElement.jacobi_matrix(zeta, eta, xe, ye);
            
            det_J = (j11 .* j22) - (j12 .* j21); 
        end
        
        function mat = get_shape_function_matrix(zeta, eta)
             mat = [ThirdOrderTriangleElement.shape_function_1(zeta, eta);
                   ThirdOrderTriangleElement.shape_function_2(zeta, eta);
                   ThirdOrderTriangleElement.shape_function_3(zeta, eta);
                   ThirdOrderTriangleElement.shape_function_4(zeta, eta);
                   ThirdOrderTriangleElement.shape_function_5(zeta, eta);
                   ThirdOrderTriangleElement.shape_function_6(zeta, eta);
                   ThirdOrderTriangleElement.shape_function_7(zeta, eta);
                   ThirdOrderTriangleElement.shape_function_8(zeta, eta);
                   ThirdOrderTriangleElement.shape_function_9(zeta, eta);
                   ThirdOrderTriangleElement.shape_function_10(zeta, eta)];
        end
        
        
    end
    
    
end

