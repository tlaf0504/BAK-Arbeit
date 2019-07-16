classdef SecondOrderTriangleElement
    
    properties(Constant)
        v11 = [4; -8; 4; 0; 0; 0];
        v12 = [4; -4; 0; 4; 0; -4];
        v13 = [-3; 4; -1; 0; 0; 0];
        v21 = [4; -4; 0; 4; 0; -4];
        v22 = [4; 0; 0; 0; 4; -8];
        v23 = [-3; 0; 0; 0; -1; 4];
        
        number_of_nodes = 6;
        
        % Assignment of element nodes to the corresponding triangle edges.
        % According to http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
        % (21.06.2019) the nodes 0, 3 and 1 belong to the bottom triangle edge.
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
        nodes_of_triangle_edges = [1,1,1,0,0,0;
                                   0,0,1,1,1,0;
                                   0,0,0,1,1,1];
                               
       % Correction of node ordering of triangles and curves.
       % Gmsh has a different node ordering than used in this program. See 
       % http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering (15.07.2019) for
       % further information.
       % Attention: Gmsh starts node numbering with 0. As matlab starts indexing with
       % 1, the gmsh node numbers also start with 1.
       
                               % Gmsh  Program
       triangle_node_ordering = [ 1,      1;
                                  2,      3;
                                  3,      5;
                                  4,      2;
                                  5,      4;
                                  6,      6];
                         % Gmsh  Program 
      curve_node_ordering = [1,    1;
                             2,    3;
                             3,    2];
    end
 
    methods(Static)
        
        % --------- Element shape functions and their derivatives
        
        % ----- Shape functions
        % Shape function of node 1
        function val = shape_function_1(zeta, eta)
           val = (1 - zeta - eta) .* (1 - 2*zeta - 2*eta);
        end
        
        % Shape function of node 2
        function val = shape_function_2(zeta, eta)
            val = 4*zeta .* (1 - zeta - eta);
        end
        
        % Shape function of node 3
        function val = shape_function_3(zeta, eta)
            val = zeta .* (2 * zeta - 1);
        end
        
        % Shape function of node 4
        function val = shape_function_4(zeta, eta)
            val = 4 .* zeta .* eta;
        end
        
        % Shape function of node 5
        function val = shape_function_5(zeta, eta)
            val = eta .* (2 * eta - 1);
        end
        
        % Shape function of node 6
        function val = shape_function_6(zeta, eta)
            val = 4*eta .* (1 - zeta - eta);
        end
        
        
        % ----- Derivatives of shape functions
        
        % --- Zeta derivatives
        function val = dN1_dZeta(zeta, eta)
            val = 4*eta + 4*zeta - 3;
        end
        
        function val = dN2_dZeta(zeta, eta)
            val = 4 - 4*eta - 8*zeta;
        end 
        
        function val = dN3_dZeta(zeta, eta)
            val = 4*zeta - 1;
        end 
        
        function val = dN4_dZeta(zeta, eta)
            val = 4*eta;
        end 
        
        function val = dN5_dZeta(zeta, eta)
            [a,b] = size(zeta);
            val = zeros(a,b);
        end 
        
        function val = dN6_dZeta(zeta, eta)
            val = -4*eta;
        end
        
        % --- Eta derivatives
        function val = dN1_dEta(zeta, eta)
            val = 4*eta + 4*zeta - 3;
        end
        
        function val = dN2_dEta(zeta, eta)
            val = -4*zeta;
        end 
        
        function val = dN3_dEta(zeta, eta)
            [a,b] = size(zeta);
            val = zeros(a,b);
        end 
        
        function val = dN4_dEta(zeta, eta)
            val = 4*zeta;
        end 
        
        function val = dN5_dEta(zeta, eta)
            val = 4*eta - 1;
        end 
        
        function val = dN6_dEta(zeta, eta)
            val = 4 - 4*zeta - 8*eta;
        end
        
        
        % ------- Misc functions
        
        function val = get_shape_function(k, zeta, eta)
            
            assert(k > 0 & k <= SecondOrderTriangleElement.number_of_nodes)
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
                
            elseif k == 4
                val = SecondOrderTriangleElement.shape_function_4( ...
                    zeta, eta);
                
            elseif k == 5
                val = SecondOrderTriangleElement.shape_function_5( ...
                    zeta, eta);
                
            elseif k == 6
                val = SecondOrderTriangleElement.shape_function_6( ...
                    zeta, eta);
            end
            
        end
        
        function det_J = jacobi_determinant(zeta, eta, xe, ye)
            if ~iscolumn(xe)
                xe = xe';
            end
            
            if ~iscolumn(ye)
                ye = ye';
            end
            
            v11 = SecondOrderTriangleElement.v11;
            v12 = SecondOrderTriangleElement.v12;
            v13 = SecondOrderTriangleElement.v13;
            v21 = SecondOrderTriangleElement.v21;
            v22 = SecondOrderTriangleElement.v22;
            v23 = SecondOrderTriangleElement.v23;

            
            a1 = v11' * xe; a2 = v12' * xe; a3 = v13' * xe;
            b1 = v21' * ye; b2 = v22' * ye; b3 = v23' * ye;
            c1 = v11' * ye; c2 = v12' * ye; c3 = v13' * ye;
            d1 = v21' * xe; d2 = v22' * xe; d3 = v23' * xe;
            
            
            a1s = a1 * b1;
            a2s = a2 * b2;
            a3s = a1*b2 + a2*b2;
            a4s = a1*b3 + a3*b1;
            a5s = a2*b3 + b3*a2;
            a6s = a3*b3;
            
            b1s = c1 * d1;
            b2s = c2 * d2;
            b3s = c1*d2 + c2*d2;
            b4s = c1*d3 + c3*d1;
            b5s = c2*d3 + c3*d2;
            b6s = c3*d3;
            
            det_J = (a1s - b1s) .* zeta.^2 + (a2s - b2s) .* eta.^2 + ...
                (a3s - b3s) .* zeta .* eta + (a4s - b4s) .* zeta + ...
                (a5s - b5s) .* eta + a6s - b6s;
            
            
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
                   SecondOrderTriangleElement.dN4_dZeta(zeta, eta);
                   SecondOrderTriangleElement.dN5_dZeta(zeta, eta);
                   SecondOrderTriangleElement.dN6_dZeta(zeta, eta)];
            
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