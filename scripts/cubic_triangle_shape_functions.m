clear all
close all
clc

% Script for calculating the shape funcitons for a third order (10-point) triangular 
% finite element
%
% We use the following Ansatz for the potential:
%
% V = c0 + c1 * zeta + c2 * eta + c3 * zeta * eta + c4 * zeta^2 + c5 * eta^2 + ...
%     c6 * zeta^2 * eta + c7 * zeta * eta^2 + c8 * zeta^3 + c9 * eta^2

% Node numbering of the triangle:
%
% 7\
% |  \
% |    \ 
% 8     6
% |       \
% |         \
% 9     10    5
% |             \
% |               \
% 1-----2-----3-----4
%

% Node coordinates:
%       (zeta, eta)
% P1 =  (0,      0)
% P2 =  (1/3,    0)
% P3 =  (2/3,    0)
% P4 =  (1,      0)
% P5 =  (2/3,  1/3)
% P6 =  (1/3,  2/3)
% P7 =  (0,      1)
% P8 =  (0,    2/3)
% P9 =  (0,    1/3)
% P10 = (1/3,  1/3)


% c0 to c9 are the unknowns. Setting up linear equation system [V] = [A] * [c] and
% solve it for [c]

% Equations: V1 = c0 + c1 * zeta_P1 + c2 * eta_P1 + ....
%            V2 = c0 + c1 * zeta_P2 + c2 * eta_P2 + ....

% zeta coordinates of nodes
zeta = [0, 1/3, 2/3, 1, 2/3, 1/3, 0, 0, 0, 1/3];

% eta coordinates of nodes
eta = [0, 0, 0, 0, 1/3, 2/3, 1, 2/3, 1/3, 1/3];

A = [
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1, zeta(2),  eta(2),  zeta(2) *  eta(2),  zeta(2)^2,  eta(2)^2,  zeta(2)^2 *  eta(2),  zeta(2) *  eta(2)^2,  zeta(2)^3,  eta(2)^3; 
    1, zeta(3),  eta(3),  zeta(3) *  eta(3),  zeta(3)^2,  eta(3)^2,  zeta(3)^2 *  eta(3),  zeta(3) *  eta(3)^2,  zeta(3)^3,  eta(3)^3;
    1, zeta(4),  eta(4),  zeta(4) *  eta(4),  zeta(4)^2,  eta(4)^2,  zeta(4)^2 *  eta(4),  zeta(4) *  eta(4)^2,  zeta(4)^3,  eta(4)^3;
    1, zeta(5),  eta(5),  zeta(5) *  eta(5),  zeta(5)^2,  eta(5)^2,  zeta(5)^2 *  eta(5),  zeta(5) *  eta(5)^2,  zeta(5)^3,  eta(5)^3;
    1, zeta(6),  eta(6),  zeta(6) *  eta(6),  zeta(6)^2,  eta(6)^2,  zeta(6)^2 *  eta(6),  zeta(6) *  eta(6)^2,  zeta(6)^3,  eta(6)^3;
    1, zeta(7),  eta(7),  zeta(7) *  eta(7),  zeta(7)^2,  eta(7)^2,  zeta(7)^2 *  eta(7),  zeta(7) *  eta(7)^2,  zeta(7)^3,  eta(7)^3;
    1, zeta(8),  eta(8),  zeta(8) *  eta(8),  zeta(8)^2,  eta(8)^2,  zeta(8)^2 *  eta(8),  zeta(8) *  eta(8)^2,  zeta(8)^3,  eta(8)^3;
    1, zeta(9),  eta(9),  zeta(9) *  eta(9),  zeta(9)^2,  eta(9)^2,  zeta(9)^2 *  eta(9),  zeta(9) *  eta(9)^2,  zeta(9)^3,  eta(9)^3;
    1,  zeta(10), eta(10), zeta(10) * eta(10), zeta(10)^2, eta(10)^2, zeta(10)^2 * eta(10), zeta(10) * eta(10)^2, zeta(10)^3, eta(10)^3;
    ];

syms zeta
syms eta

polynomials = [1, zeta, eta, zeta * eta, zeta^2, eta^2, zeta^2 * eta, zeta * eta^2, zeta^3, eta^3];

A_inv = inv(A);

% Set 'zero' values unequal to zero due to numerical inaccuracies to zero
% All abs values in matrix are at least >= 1. Only the very small values are below
% this level. --> Use threshold level of 0.5.
A_inv(abs(A_inv) < 0.5) = 0;


% ===== Calculate shape functions
shape_functions = polynomials * A_inv;

% ===== Calculate partial zeta and eta derivatives
shape_functions_zeta_derivatives = diff(shape_functions, zeta);
shape_functions_eta_derivatives = diff(shape_functions, eta);


% ===== Generate latex expressions and store to file

shape_functions_latex_expressions = {};
shape_functions_zeta_derivatives_latex_expressions = {};
shape_functions_eta_derivatives_latex_expressions = {};
for k = 1 : 10
    shape_functions_latex_expressions{k} = latex(shape_functions(k));
    
    shape_functions_zeta_derivatives_latex_expressions{k} = ...
        latex(shape_functions_zeta_derivatives(k));
    
    shape_functions_eta_derivatives_latex_expressions{k} = ...
        latex(shape_functions_eta_derivatives(k));
end

% save('third_order_triangle_shape_functions_latex_expressions', ...
%     'shape_functions_latex_expressions', ...
%     'shape_functions_zeta_derivatives_latex_expressions', ...
%     'shape_functions_eta_derivatives_latex_expressions');

save('third_order_triangle_shape_functions.mat', 'shape_functions', ...
    'shape_functions_zeta_derivatives', 'shape_functions_eta_derivatives')



