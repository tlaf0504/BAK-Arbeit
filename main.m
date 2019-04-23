clear all
close all
clc

addpath('entities')

elem = SecondOrderTriangleElement([0,1/2,1,1/2,0,0], [0,0,0,1/2,1,1/2]);
fun = @(zeta, eta) ElectrostaticProblem.element_matrix_integrant(...
        elem, 1 , 2, zeta, eta, 0.00001, 0.0005);
val = GaussianQuadrature.integrate_2D_normalized_triangle_region(fun,7)

