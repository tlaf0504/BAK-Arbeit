clear all
close all
clc


%% DESCRIPTION
%
% Script for calculation of the coefficients used for caluclation of the
% shape functions
%
%%

syms u1 u2 u3 u4 u5 u6
u = [u1; u2; u3; u4; u5; u6];
A = [
     1, 0,   0,   0,   0,   0   ;
     1, 1/2, 0,   1/4, 0,   0   ;
     1, 1,   0,   1,   0,   0   ;
     1, 1/2, 1/2, 1/4, 1/4, 1/4 ;
     1, 0,   1,   0,   0,   1   ;
     1, 0,   1/2, 0,   0,   1/4 
     ];
 
 
% Coefficients for shape functions 
c = inv(A) * u