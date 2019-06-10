classdef GaussianQuadrature
    
    properties(Constant)
        
        % Three-point, seven-point and fiveteen-point integration for
        % normalized triangular regions (fixed interval limits of [0,1]) 
        % and line integration with one to six points for arbitrary limits 
        % are implemented.
        
        possible_numbers_of_triangle_integration_points = [3,7,15];
        possible_numbers_of_line_integration_points = 1:6;
        
        
        
        % ============ Constants for triangle integration =================
        
        % Coordinates and weights of integration points for 7 and 3 point
        % gaussian quadrature of a trangle region with normalized
        % corner-coordinates. (Triangle corners habe coordinated of (0,0),
        % (0,1) and (1,0)
        %
        % Due to the normalization the integraion has fixes limits of [0,1]
        % in both directions.
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
                 
                  
                  
        % ============== Constants for line integration ===================
        
        % Coordinates and weights of integration points for 1- to 6-point
        % gaussian quadrature of a line with arbitrary integraion limits
        %
        % Table was taken from K.-J.Bathe "Finite Elemente Methoden", 2nd
        % edition, Springer 2002; S.542
        % 
        
        % One-point integration
        r_line_1_point = 0;
        w_line_1_point = 2;
        
        % Two-point integration
        r_line_2_point = [0.577350269189626; -0.577350269189626];
        w_line_2_point = [1;1];
        
        % Three-point integration
        r_line_3_point = [0.774596669241483; 0; -0.774596669241483];
        w_line_3_point = [5/9; 8/9; 5/9];
        
        % Four-point integration
        r_line_4_point = [0.861136311594053; 0.339981043584856; ...
            -0.339981043584856; -0.861136311594053];
        w_line_4_point = [0.347854845137454; 0.652145154862546; ...
            0.652145154862546; 0.347854845137454];
        
        % Five-point integration
        r_line_5_point = [0.906179845938664; 0.538469310105683; ...
            0; -0.538469310105683; -0.906179845938664];
        w_line_5_point = [0.236926885056189; 0.478628670499366; ...
            128/225; 0.478628670499366; 0.236926885056189];
        
        % Six-point integration
        r_line_6_point = [0.932469514203152; 0.661209386466265; ...
            0.238619186083197; -0.238619186083197; -0.661209386466265; ...
            -0.932469514203152];
        w_line_6_point = [0.171324492379170; 0.360761573048139; ...
            0.467913934572691; 0.467913934572691; 0.360761573048139; ...
            0.171324492379170];
        
                 
    end
    
    
    methods(Static)
        
        function val = integrate_2D_normalized_triangle_region(...
                integrant, number_of_integration_points)
            
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
            
            val = 1/2 * (integrant(zeta_vec', eta_vec') * w_vec);     
        end
        
        
        
        function val = integrate_1d_line(integrant, ...
                number_of_integration_points, upper_limit, lower_limit)
            
            N_allowed = ...
                GaussianQuadrature. ...
                possible_numbers_of_line_integration_points;
            
            assert(ismember(number_of_integration_points, N_allowed), ...
                ['Assertion failed.\nNumber of integration points must ' ...
                'be one of the following values: %s\n', ...
                'You entered %d.'], ...
                sprintf('%d ', N_allowed), number_of_integration_points);
            
            
            % Coordinates and weights are specified for a nomalized
            % integration interval of [-1, +1].
            % For arbitray limits the corrdinates and weights must be
            % transformed to the input limits.
            %
            % According to K.-J.Bathe "Finite Elemente Methoden", 2nd
            % edition, Springer 2002; S.541 the the transformation can be 
            % archieved as following:
            %
            % r_trans = (a+b)/2 + (b-a)/2 * r_i
            % w_trans = (b-a)/2 * w_i
            %
            % where r_i is the coordinate and w_i is the weight in the
            % normalized interval, r_trans and w_trans are the transformed
            % parameters in the arbitrary interval and 'a' and 'b' are the
            % limits of the arbitrary interval.
            %
            coordinate_offset = (upper_limit + lower_limit) / 2;
            coordinate_scaler = (upper_limit - lower_limit) / 2;
            weight_scaler = coordinate_scaler;
            
            
            switch number_of_integration_points
                case 1
                    r_vec = GaussianQuadrature.r_line_1_point;
                    w_vec = GaussianQuadrature.w_line_1_point;
                    
                case 2
                    r_vec = GaussianQuadrature.r_line_2_point;
                    w_vec = GaussianQuadrature.w_line_2_point;
                    
                    
                case 3
                    r_vec = GaussianQuadrature.r_line_3_point;
                    w_vec = GaussianQuadrature.w_line_3_point;
                    
                    
                case 4
                    r_vec = GaussianQuadrature.r_line_4_point;
                    w_vec = GaussianQuadrature.w_line_4_point;
                    
                case 5
                    r_vec = GaussianQuadrature.r_line_5_point;
                    w_vec = GaussianQuadrature.w_line_5_point;
                    
                case 6
                    r_vec = GaussianQuadrature.r_line_6_point;
                    w_vec = GaussianQuadrature.w_line_6_point;
            end
            
            coordinates = r_vec * coordinate_scaler + coordinate_offset;
            coordinates = coordinates(:); % Make column vector
            weights = w_vec * weight_scaler;
            
            val = integrant(coordinates') * weights;
        end
        
    end
end

