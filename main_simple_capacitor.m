clear all
close all
clc

addpath('entities')


V = 5;
Ne = 4; % Number of elements
N = 6; % Number of nodes of each element
N_glob_nodes = 15;

% epsilon_0
eps = 8.854178*10^(-12);

% Simple plate capacitor with 6 triangular finite elements of second order


x = [0; 0; 0; 0.5; 1; 1; 1; 0.5; 0.5; 0; 0; 0.5; 1; 1; 0.5];
y = [1; 0.75; 0.5; 0.5; 0.5; 0.75; 1; 1; 0.75; 0.25; 0; 0; 0; 0.25; 0.25];

node_coordinates = [x, y];
                
              
% Mapping of local nodes to global nodes.
node_mapping = [1,  2,  3,  9,  7, 8;
                3,  4,  5,  6,  7, 9;
                3,  10, 11, 15, 5, 4;
                11, 12, 13, 14, 5, 15];
            
            
% Mapping of direchlet boundaries
%
% Column 1 denote the global nodes and column 2 the corresponding direchlet
% boundary condition
direchlet_boundary_conditions = [1, V;
                                 8, V;
                                 7, V;
                                 11, 0;
                                 12, 0;
                                 13, 0];
                             
[N_dir, ~] = size(direchlet_boundary_conditions);            
N_unknown = N_glob_nodes - N_dir;
      
% Coefficient matrix of assembled equation system            
A = sparse(N_unknown, N_unknown);

% Right-side vector
r = zeros(N_unknown, 1);

global_node_to_equation_system_index_mapping = ...
    solver.get_node_index_to_equation_index_mapping( ...
                N_glob_nodes, direchlet_boundary_conditions(:,1));

                             
for e = 1 : Ne
    local_node_mapping = node_mapping(e, :);
    local_node_mapping = local_node_mapping(:);
    
    xe = node_coordinates(local_node_mapping, 1);
    ye = node_coordinates(local_node_mapping, 2);
    
    
    % Local coefficient matrix and right-side vector
    A_loc = zeros(N, N);
    r_loc = zeros(N, 1);
    
    row_idx = 1 : N;
    row_idx = row_idx(:);
    
    col_idx = 1 : N;
    col_idx = col_idx(:);
    
    % N_dir x 1 vector containing logal 1s at elements corresponding to a local node
    % of the current finite element
    direchlet_conditions_of_currnet_element_bin_idx = ...
        ismember(direchlet_boundary_conditions(:,1), local_node_mapping);
    
    local_nodes_with_direchlet_conditions_bin_idx = ...
        ismember(local_node_mapping, direchlet_boundary_conditions(:,1));
    
    direchlet_boundary_data = [...
        row_idx(local_nodes_with_direchlet_conditions_bin_idx), ...
        direchlet_boundary_conditions(direchlet_conditions_of_currnet_element_bin_idx, :)];
    
    [number_of_direchlet_conditions_within_element, ~] = ...
        size(direchlet_boundary_data);
    
    % A direchlet boundary condition at local node k removes the k-th (the k-th
    % equation) of the element-equation system.
    row_idx(local_nodes_with_direchlet_conditions_bin_idx) = [];
    
    for l = 1 : length(row_idx)
        
        row = row_idx(l);
        
        % Normally the right-side vector-elements would be calculates in here, but we
        % assume no free charges and no inhomogeneous Neumann boundary conditions
        
        for m = 1 : length(col_idx)
            
            col = col_idx(m);
            
            fun = @(zeta, eta) ElectrostaticProblem.element_matrix_integrant(...
                SecondOrderTriangleElement, xe, ye, row , col, zeta, eta, eps, eps);
            coef = GaussianQuadrature.integrate_2D_normalized_triangle_region(fun,7);
            
            % If the k-th local node contains a direchlet boundary condition, the
            % coeffient in the k-th column multiplied with the value of the boundary
            % condition can be brought to the right side of the current equation
            if ismember(col, direchlet_boundary_data(:,1))
                idx = direchlet_boundary_data(:, 1) == col;
                r_loc(row) = r_loc(row) - coef * direchlet_boundary_data(idx, 3);
            else
                A_loc(row, col) = coef;
            end
        end
    end
    
    col_idx(local_nodes_with_direchlet_conditions_bin_idx) = [];
    
    global_row_idx = global_node_to_equation_system_index_mapping(local_node_mapping(row_idx));
    global_col_idx = global_node_to_equation_system_index_mapping(local_node_mapping(col_idx));
    
    A(global_row_idx, global_col_idx) = A(global_row_idx, global_col_idx) + ...
        A_loc(row_idx, col_idx);
    
    r(global_row_idx) = r(global_row_idx) + r_loc(row_idx);
    
end

V_unknown = A\r;
V_unknown

