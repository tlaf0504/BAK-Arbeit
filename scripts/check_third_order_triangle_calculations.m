clear all
close all
clc

addpath('../entities')

load('third_order_triangle_shape_functions.mat', 'shape_functions', 'shape_functions_zeta_derivatives', ...
    'shape_functions_eta_derivatives')


syms zeta eta

xe = [0, 1/6, 1/3, 1/2, 1/3, 1/6, 0, 0, 0, 1/2];
ye = [0, 1/6, 1/3, 1/2, 2/3, 5/6, 1, 2/3, 1/3, 1/2];

J = [ shape_functions_zeta_derivatives * xe', shape_functions_zeta_derivatives * ye';
      shape_functions_eta_derivatives * xe', shape_functions_eta_derivatives * ye'];
  
J_inv = inv(J);
               
              
  
  
  
  %%%%%%%
  det_J_num = ThirdOrderTriangleElement.jacobi_determinant(1/3, 2/99,xe, ye)
  
  det_J_sym = subs(det(J), [zeta, eta], [1/3, 2/99])
  
  %%%%%%%
  k = 1; 
  l = 1;
  
  zeta_num = 0;
  eta_num = 0;
  
  kij_integrant_sym = (J_inv(1,1) * shape_functions_zeta_derivatives(k) + J_inv(1,2) * shape_functions_eta_derivatives(k)) * ...
            (J_inv(1,1) * shape_functions_zeta_derivatives(l) + J_inv(1,2) * shape_functions_eta_derivatives(l)) +...
            ...
            (J_inv(2,1) * shape_functions_zeta_derivatives(k) + J_inv(2,2) * shape_functions_eta_derivatives(k)) * ...
            (J_inv(2,1) * shape_functions_zeta_derivatives(l) + J_inv(2,2) * shape_functions_eta_derivatives(l));
        
  kij = subs(kij_integrant_sym, [zeta, eta], [zeta_num, eta_num])
        
  ElectrostaticProblem.element_matrix_integrant(ThirdOrderTriangleElement, xe, ye, k , l, ...
                zeta_num, eta_num, 1, 1)
        
        
            