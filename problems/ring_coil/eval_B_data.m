clear al l
close all
clc

load('B_data.mat');

I = 1000;
mu_0 = 4*pi*10^(-7);
R = 0.25;

R_vec = [R, 0];

B_ideal = @(z) mu_0/2 * R^2/(R^2 + z^2)^(3/2) * I;

N_data_points = length(abs_magnetic_field);

B_ideal_data = zeros(N_data_points, 1);
z_data = zeros(N_data_points, 1);

for k = 1 : N_data_points
    coordinates_k = coordinates(k, :);
    
    z_data(k) = norm(coordinates_k - R_vec);
    
    B_ideal_data(k) = B_ideal(z_data(k));
    
end
scale = max(abs_magnetic_field) / max(B_ideal_data);

figure()
%plot(z_data, scale * B_ideal_data, 'o');
%hold on
plot(z_data,abs_magnetic_field, 'o')