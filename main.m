clear all
close all
clc

addpath('entities')

a = [0,1/2,1];
b = [0,1/2,1];

xe = [0,1/2,1,1/2,0,0];
ye = [0,0,0,1/2,1,1/2];
% Calculation of one element of element matrix
%
% Compute element k_1_2 of element hatrix using epsilon_x = epsilon_y = 1
fun = @(zeta, eta) ElectrostaticProblem.element_matrix_integrant(...
        SecondOrderTriangleElement, xe, ye, 1 , 1, zeta, eta, 1, 1);
    
fprintf('Element k_12 of element matrix:\nk_12 = %f\n\n', ...
    GaussianQuadrature.integrate_2D_normalized_triangle_region(fun,7))

% Calculation of one element neumann-element of the right-side vector
% 
% Compute line integral of shape function 5 on side 3 with sigma = 1
fun = @(t) ElectrostaticProblem.right_side_neumann_integrant(...
    SecondOrderTriangleElement, xe, ye, 5, t, 3, 1);

fprintf('Neuman part or right-side vector element 5:\n%f\n\n', ...
    GaussianQuadrature.integrate_1d_line(fun, 5, 1, 0))

%arrayfun(@(a) GaussianQuadrature.integrate_1d_line(fun, 5, 1, 0), 0:100000);

% Calculation of one plane-element of right-side vector
%
% Compute plane integral of second node shape function
fun = @(zeta, eta) ElectrostaticProblem.right_side_plane_integrant(...
    SecondOrderTriangleElement, xe, ye, 2, zeta, eta, 1);

fprintf('Plane part or right-side vector element 5:\n%f\n\n', ...
    GaussianQuadrature.integrate_2D_normalized_triangle_region(fun,7))



% Compute complete equation system for one element
fprintf('======= Calculation of element-matrix and right-side vector for one element =========\n\n')
A = ones(6,6);
r = ones(6,1);

eps = 8.854178*10^(-12);
sigma = 1
roh = 1;

fprintf('Elapsed time for calculation of one element:\n')
tic
for k = 1 : 6
    for l = 1 : 6
        fun = @(zeta, eta) ElectrostaticProblem.element_matrix_integrant(...
            SecondOrderTriangleElement, xe, ye, k , l, zeta, eta, eps, eps);
        A(k,l) = GaussianQuadrature.integrate_2D_normalized_triangle_region(fun,7);
    end
    
    fun = @(zeta, eta) ElectrostaticProblem.right_side_plane_integrant(...
        SecondOrderTriangleElement, xe, ye, k, zeta, eta, roh);
    
    r(k) = GaussianQuadrature.integrate_2D_normalized_triangle_region(fun,7);
    
    % Assumung neumann-condition on triangle side 3
    fun = @(t) ElectrostaticProblem.right_side_neumann_integrant(...
    SecondOrderTriangleElement, xe, ye, k, t, 3, sigma);

    r(k) = r(k) + GaussianQuadrature.integrate_1d_line(fun, k, 1, 0);
    
end
toc
fprintf('\n')
fprintf('Element-matrix A is:\n')
disp(A)

fprintf('Right-side vector r is :\n')
disp(r)
        