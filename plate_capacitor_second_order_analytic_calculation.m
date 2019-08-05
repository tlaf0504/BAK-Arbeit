clear all
close all
clc

xe_num = [1, 3/4, 0.5, 3/4, 1, 1];
ye_num = [1, 7/8, 3/4, 5/8, 0.5, 3/4];

N_unknowns = 17;
A_glob = zeros(N_unknowns, N_unknowns);
r_glob = zeros(N_unknowns);



node_1 = [0 0];
node_2 = [1 0];
node_3 = [1 0.5];
node_4 = [1 1];
node_5 = [0 0.5];
node_6 = [0 1];
node_7 = [0 0.75];
node_8 = [0 0.25];
node_9 = [0.4999999999999999 0];
node_10 = [1 0.25];
node_11 = [1 0.75];
node_12 = [0.5 1];
node_13 = [0.4999999999999999 0.5];
node_14 = [0.5 0.75];
node_15 = [0.75 0.875];
node_16 = [0.75 0.625];
node_17 = [0.25 0.625];
node_18 = [0.25 0.875];
node_19 = [0.5 0.25];
node_20 = [0.25 0.375];
node_21 = [0.25 0.125];
node_22 = [0.75 0.125];
node_23 = [0.75 0.375];

node_coordinates = [ ...
    node_1; node_2; node_3; node_4; node_5; node_6; node_7; node_8; node_9; ...
    node_10; node_11; node_12; node_13; node_14; node_15; node_16; node_17; ...
    node_18; node_19; node_20; node_21; node_22; node_23];

element_node_mapping = [ ...
    4 15 14 16 3 11;
    5 17 14 18 6 7;
    3 16 14 17 5 13;
    6 18 14 15 4 12;
    5 20 19 21 1 8;
    2 22 19 23 3 10;
    1 21 19 22 2 9;
    3 23 19 20 5 13];




% Analytic calculation for problem 'plate capacitor multimaterial'

syms zeta eta
syms x1 x2 x3 x4 x5 x6
syms y1 y2 y3 y4 y5 y6


N1 = (1 - zeta - eta) * (1 - 2*zeta - 2*eta);
N2 = 4*zeta * (1 - zeta - eta);
N3 = zeta * (2 * zeta - 1);
N4 = 4 * zeta * eta;
N5 = eta * (2 * eta - 1);
N6 = 4*eta * (1 - zeta - eta);

dN1_dZeta = diff(N1, zeta);
dN2_dZeta = diff(N2, zeta);
dN3_dZeta = diff(N3, zeta);
dN4_dZeta = diff(N4, zeta);
dN5_dZeta = diff(N5, zeta);
dN6_dZeta = diff(N6, zeta);

dN1_dEta = diff(N1, eta);
dN2_dEta = diff(N2, eta);
dN3_dEta = diff(N3, eta);
dN4_dEta = diff(N4, eta);
dN5_dEta = diff(N5, eta);
dN6_dEta = diff(N6, eta);

dN_dZeta = [dN1_dZeta; dN2_dZeta; dN3_dZeta; dN4_dZeta; dN5_dZeta; dN6_dZeta];
dN_dEta = [dN1_dEta; dN2_dEta; dN3_dEta; dN4_dEta; dN5_dEta; dN6_dEta];

xe = [x1; x2; x3; x4; x5; x6];
ye = [y1; y2; y3; y4; y5; y6];

J = [dN_dZeta.' * xe, dN_dZeta.' * ye;
    dN_dEta.' * xe, dN_dEta.' * ye];

det_J = det(J);

inv_J = inv(J);

% Integrant of plane-integral for element matrices: Integrant from electrostatic
% problem used. epsilon_x = epsilon_y = 1 for simplification.

F_kl = @(k,l, zeta, eta) ...
    ((dN_dZeta(k) * inv_J(1,1) + dN_dEta(k) * inv_J(1,2)) * ...
    (dN_dZeta(l) * inv_J(1,1) + dN_dEta(l) * inv_J(1,2)) + ...
    (dN_dZeta(k) * inv_J(2,1) + dN_dEta(k) * inv_J(2,2)) * ...
    (dN_dZeta(l) * inv_J(2,1) + dN_dEta(l) * inv_J(2,2))) * det_J;

A_loc = zeros(6,6);

[number_of_elements, ~] = size(element_node_mapping);

element_matrix_struct = struct();

for element_idx = 1 : number_of_elements
    nodes_of_current_element = element_node_mapping(element_idx, :);
    node_coordinates_of_current_element = node_coordinates(nodes_of_current_element, :);
    
    xe_num = node_coordinates_of_current_element(:, 1);
    ye_num = node_coordinates_of_current_element(:, 2);
    
    for k = 1 : 6
        for l = 1 : 6
            
            tmp = subs(F_kl(k, l), [x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6], ...
                [xe_num(1), xe_num(2), xe_num(3), xe_num(4), xe_num(5), xe_num(6),  ...
                ye_num(1), ye_num(2), ye_num(3), ye_num(4), ye_num(5), ye_num(6)]);
            
            val =  int(int(tmp, eta, 0, 1 - zeta), zeta, 0, 1);
            A_loc(k,l) = val;
        end
    end
    fprintf('====================\nElement matrix of element %d:\n', element_idx);
    disp(A_loc)
    fprintf('====================\n\n')
    
    element_matrix_struct.(sprintf('element_matrix_%d', element_idx)) = A_loc;
    
end

save('element_matrices', '-struct', 'element_matrix_struct')








