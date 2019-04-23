classdef SecondOrderTriangleElement
    
    properties(Constant)
        v11 = [4; -8; 4; 0; 0; 0];
        v12 = [4; -4; 0; 4; 0; -4];
        v13 = [-3; 4; -1; 0; 0; 0];
        v21 = [4; -4; 0; 4; 0; -4];
        v22 = [4; 0; 0; 0; 4; -8];
        v23 = [-3; 0; 0; 0; -1; 4];
        
        number_of_nodes = 6;
    end
    
    properties(Access = public)
        % Node coordinates
        xe
        ye
        
    end
    
    properties(Access = private)    
    end
    
    methods
        % Getter and Setter
        function obj = set.xe(obj, values)
            assert(length(values) == obj.number_of_nodes, ...
                ['Assetion failed. \nWanted to set variable "xe" to a ',...
                'vector of length %d, but only vectors of length %d ', ...
                'are allowed'], length(values), obj.number_of_nodes);
            obj.xe = values(:);
        end
        
        function obj = set.ye(obj, values)
            assert(length(values) == obj.number_of_nodes, ...
                ['Assetion failed. \nWanted to set variable "ye" to a ',...
                'vector of length %d, but only vectors of length %d ', ...
                'are allowed'], length(values), obj.number_of_nodes);
            obj.ye = values(:);
        end
    end
    
    methods(Access = public)
        
        % Construtor of object
        function obj = SecondOrderTriangleElement(xe, ye)
            obj.xe = xe(:);
            obj.ye = ye(:);
        end
        
        function det_J = jacobi_determinant(obj, eta, zeta)
            v11 = obj.v11;
            v12 = obj.v12;
            v13 = obj.v13;
            v21 = obj.v21;
            v22 = obj.v22;
            v23 = obj.v23;
            
            xe = obj.xe;
            ye = obj.ye;
            
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
            
            det_J = (a1s - b1s) * zeta^2 + (a2s - b2s) * eta^2 + ...
                (a3s - b3s) * zeta * eta + (a4s - b4s) * zeta + ...
                (a5s - b5s) * eta + a6s - b6s;
            
            
        end
        
        function J_inv = inverse_jacobi_matrix(obj, zeta, eta)
           j11 = obj.ye' * obj.v21 * zeta + obj.ye' * obj.v22 * eta + ...
               obj.ye' * obj.v23;
           
           j12 = -(obj.ye' * obj.v11 * zeta + obj.ye' * obj.v12 * eta + ...
               obj.ye' * obj.v13);
           
           j21 = -(obj.xe' * obj.v21 * zeta + obj.xe' * obj.v22 * eta + ...
               obj.xe' * obj.v23);
           
           j22 = obj.xe' * obj.v11 * zeta + obj.xe' * obj.v12 * eta + ...
               obj.xe' * obj.v13;
           
           J_inv = (1 / obj.jacobi_determinant(zeta, eta)) * ...
               [j11, j12; j21, j22];
        end
        
        function val = get_zeta_shape_function_derivatives(obj,zeta, eta)
            val = [obj.dN1_dZeta(zeta, eta);
                   obj.dN2_dZeta(zeta, eta);
                   obj.dN3_dZeta(zeta, eta);
                   obj.dN4_dZeta(zeta, eta);
                   obj.dN5_dZeta(zeta, eta);
                   obj.dN6_dZeta(zeta, eta)];
        end
        
        function val = get_eta_shape_function_derivatives(obj,zeta, eta)
            val = [obj.dN1_dEta(zeta, eta);
                   obj.dN2_dEta(zeta, eta);
                   obj.dN3_dEta(zeta, eta);
                   obj.dN4_dEta(zeta, eta);
                   obj.dN5_dEta(zeta, eta);
                   obj.dN6_dEta(zeta, eta)];
        end
        
        
        
    end
    
    methods(Static)
        
        % --------- Element shape functions and their derivatives
        
        % ----- Shape functions
        % Shape function of node 1
        function val = shape_function_1(zeta, eta)
           val = (1 - zeta - eta) * (1 - 2*zeta - 2*eta);
        end
        
        % Shape function of node 2
        function val = shape_function_2(zeta, eta)
            val = 4*zeta * (1 - zeta - eta);
        end
        
        % Shape function of node 3
        function val = shape_function_3(zeta, eta)
            val = zeta * (2 * zeta - 1);
        end
        
        % Shape function of node 4
        function val = shape_function_4(zeta, eta)
            val = 4 * zeta * eta;
        end
        
        % Shape function of node 5
        function val = shape_function_5(zeta, eta)
            val = eta * (2 * eta - 1);
        end
        
        % Shape function of node 6
        function val = shape_function_6(zeta, eta)
            val = 4*eta * (1 - zeta - eta);
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
            val = 0;
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
            val = 0;
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
        
        
    end
    
end