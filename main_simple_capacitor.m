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

% Permettivities in elemens
permettivities = eps * ones(Ne, 2);

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
                             
% Neumann boundary conditions
% 
% Column 1; Element number
% Column 2: Triangle side
% Column 3: Neumann condition
neumann_boundary_conditions = [1, 1, 0;
                               3, 2, 0];

                             
[N_dir, ~] = size(direchlet_boundary_conditions);            
N_unknown = N_glob_nodes - N_dir;
      
% Coefficient matrix of assembled equation system            
A = sparse(N_unknown, N_unknown);

% Right-side vector
r = zeros(N_unknown, 1);

global_node_to_equation_system_index_mapping = ...
    Solver.get_node_index_to_equation_index_mapping( ...
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
    
    % Vector containing logical ones at direchlet boundary conditions of current
    % element
    % Dimenison : N_dir (Number of direchlet boudnary conditions) x 1
    direchlet_conditions_of_currnet_element_bin_idx = ...
        ismember(direchlet_boundary_conditions(:,1), local_node_mapping);
    
    % Vector containing logical ones at local nodes with direchlet boundary
    % conditions
    % Dimenstion N (Number of element nodes) x 1
    local_nodes_with_direchlet_conditions_bin_idx = ...
        ismember(local_node_mapping, direchlet_boundary_conditions(:,1));
    
    % First column: Local nodes with direchlet boundary conditions
    % Second and third column: Global node + boundary conditions corresponding to
    % local node in first column.
    %
    % E.g.:
    %
    % 1 25 10
    % 3 73 0
    %
    % Local node 1 corresponds to global node 25 with direchlet boudnary condition of
    % 10 (V, A, etc.). Local node 3 corresponds to global node 73 with a boundary
    % condition of 0
    %
    direchlet_boundary_data = [...
        row_idx(local_nodes_with_direchlet_conditions_bin_idx), ...
        direchlet_boundary_conditions(...
        direchlet_conditions_of_currnet_element_bin_idx, :)];
    
    [number_of_direchlet_conditions_within_element, ~] = ...
        size(direchlet_boundary_data);
    
    % A direchlet boundary condition at local node k removes the k-th row(the k-th
    % equation) of the element-equation system. The elements of A_loc will be zero in
    % this row.
    row_idx(local_nodes_with_direchlet_conditions_bin_idx) = [];
    
    % Due to line above, the algorithm will only iterate above the elements without
    % direchlet boundary conditions.
    for l = 1 : length(row_idx)
        
        row = row_idx(l);
        
        % Normally the right-side vector-elements would be calculated in here, but we
        % assume no free charges and no inhomogeneous Neumann boundary conditions in
        % this example.
        
        for m = 1 : length(col_idx)
            
            col = col_idx(m);
            
            % ===== Calculation of local equation system =====
            
            fun = @(zeta, eta) ElectrostaticProblem.element_matrix_integrant(...
                SecondOrderTriangleElement, xe, ye, row , col, zeta, eta, eps, eps);
            coef = GaussianQuadrature.integrate_2D_normalized_triangle_region(fun,7);
            
            % If the k-th local node contains a direchlet boundary condition, the
            % coeffient in the k-th column multiplied with the value of the boundary
            % condition is brought to the right side of the current equation
            if ismember(col, direchlet_boundary_data(:,1))
                idx = direchlet_boundary_data(:, 1) == col;
                r_loc(row) = r_loc(row) - coef * direchlet_boundary_data(idx, 3);
            else
                A_loc(row, col) = coef;
            end
        end
    end
    
    % ===== Add current local equation system to global equation system =====
    
    % If the k-th local node has a direchlet boundary condition, the elements from
    % the k-th column should not be included in the global equation system as they
    % are already condidered in the right-side vector.
    col_idx(local_nodes_with_direchlet_conditions_bin_idx) = [];
    
    % global_node_to_equadion_system_mapping maps the global nodes to their
    % corresponding row and column within the global equation system.
    % E.g. the unknown potential of global node 6 maps to the 5th row/column of the 
    % global equation system as, for example, the 3rd global node has a direchlet
    % boundary condition. 
    %
    % loca_node_mapping maps the local node to the global nodes and row_idx and
    % col_idx contains the indices of the equation system to copy.
    global_row_idx = global_node_to_equation_system_index_mapping(local_node_mapping(row_idx));
    global_col_idx = global_node_to_equation_system_index_mapping(local_node_mapping(col_idx));
    
    A(global_row_idx, global_col_idx) = A(global_row_idx, global_col_idx) + ...
        A_loc(row_idx, col_idx);
    
    r(global_row_idx) = r(global_row_idx) + r_loc(row_idx);
    
end

V_unknown = A\r;
V_unknown

% ===== Same procedure as above, just with fully implemented algorithm =====

material_properties = eps * ones(Ne, 2);
volume_charges = 0 * ones(Ne, 1);
problem_type = 'E-Static';
plane_integration_order = 7;
curve_integration_order = 5;

[A, r] = Solver.assemble(node_coordinates, node_mapping, ...
                material_properties, direchlet_boundary_conditions, ...
                neumann_boundary_conditions, volume_charges, problem_type, ...
                plane_integration_order, curve_integration_order, 'Second order');

V_unknown2 = A\r         

diff = V_unknown2 - V_unknown