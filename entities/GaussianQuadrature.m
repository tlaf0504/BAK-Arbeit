classdef GaussianQuadrature
    
    properties(Constant)
        
        % Three-point and seven-point integration is implemented
        possible_numbers_of_integration_points = [3,7];
        
        % Coordinates and weights of integration points for 7 point
        % gaussian quadrature of a trangle region with normalized
        % corner-coordinates
        %
        % Table was taken from K.-J.Bathe "Finite Elemente Methoden", 2nd
        % edition, Springer 2002; S.547
        
        
        % Vectors for three-point integration
        %
        % Coordinate 'r' in book corresponds to coordinate 'zeta' in this
        % algorithm
        
        zeta_3_point = [1/6;
                        2/3;
                        1/6];
          
        eta_3_point = [GaussianQuadrature.zeta_3_point(1);
                       GaussianQuadrature.zeta_3_point(1);
                       GaussianQuadrature.zeta_3_point(2)];
                   
       w_3_point = [1/3;
                    1/3;
                    1/3];
        
        
        % Vectors for seven-point integration
        
        zeta_7_point = [0.1012865073235;
                        0.7974269853531;
                        0.1012865073235;
                        0.4701420641051;
                        0.4701420641051;
                        0.0597158717898;
                        1/3];
                
        eta_7_point = [GaussianQuadrature.zeta_7_point(1);
                       GaussianQuadrature.zeta_7_point(1);
                       GaussianQuadrature.zeta_7_point(2);
                       GaussianQuadrature.zeta_7_point(6);
                       GaussianQuadrature.zeta_7_point(4);
                       GaussianQuadrature.zeta_7_point(4);
                       GaussianQuadrature.zeta_7_point(7)];
                   
        w_7_point = [ 0.1259391805448;
                      0.125939180544;
                      0.125939180544;
                      0.1323941527885;
                      0.1323941527885;
                      0.1323941527885;
                      0.225];
                 
                 
    end
    
    
    methods(Static)
        
        function val = integrate_2D_normalized_triangle_region(...
                integrant, number_of_integration_points)
            
            assert(ismember(number_of_integration_points, [3,7]));
            
            N = number_of_integration_points;
            if N == 7
                zeta_vec = GaussianQuadrature.zeta_7_point;
                eta_vec = GaussianQuadrature.eta_7_point;
                w_vec = GaussianQuadrature.w_7_point;
                
            else
                zeta_vec = GaussianQuadrature.zeta_3_point;
                eta_vec = GaussianQuadrature.eta_3_point;
                w_vec = GaussianQuadrature.w_3_point;
                
            end
            
            val = 0;
            for k = 1 : N
                val = val + 1/2 * w_vec(k) * ...
                    integrant(zeta_vec(k), eta_vec(k));
            end
            
        end
        
    end
end

