classdef Postprocessor
    
    methods(Static)
        function plot_electrostatic_potential(problem_location)
            
            tmp = pwd;
            cd(problem_location)
            load('problem_setup.mat', 'problem_setup')
            load(fullfile(problem_location, 'internal', 'mesh_data.mat'), 'mesh_data');
            load(fullfile(problem_location, 'results', 'results.mat'), 'unknowns', ...
                'result_to_global_node_mapping');
            
            % Merge unknowns and dirichlet boundary conditions
            node_potentials = Postprocessor.get_node_potentials_struct(mesh_data, ...
                unknowns, ...
                result_to_global_node_mapping);
           
            % Plot mesh
            mesh_plot = Postprocessor.get_mesh_plot_data(mesh_data, problem_setup);

            Postprocessor.plot_node_potentials(node_potentials, mesh_plot);
            
            Postprocessor.calculate_energy(mesh_data, node_potentials.potentials, ...
                problem_setup);
            
            Postprocessor.plot_electric_field(mesh_data, node_potentials.potentials, ...
                problem_setup, mesh_plot)
            
            
            
            
            cd(tmp);
        end
         
        
        function calculate_energy(mesh_data, node_potentials, problem_setup)
            
            % Get prblem type class
            [problem_type, success] = Misc. ...
                get_problem_type_class_from_problem_type_number( ...
                problem_setup.problem_type);
            if ~success
                return
            end
            
            finite_element_class = ...
                Misc.get_finite_element_class_from_mesh_order(problem_setup.mesh_order);
            
            
            number_of_finite_elements = mesh_data.triangle_data.number_of_triangles;
            
            material_values = mesh_data.material_and_source_properties. ...
                material_values;
            
            vacuum_material = Misc.get_vacuum_material(problem_setup.problem_type);
            
            
            
            energy = 0;
            
            for k = 1 : number_of_finite_elements
                
                nodes_of_current_triangle = mesh_data.triangle_data.nodes(k, :);
                node_potentials_of_current_triangle = node_potentials(...
                    nodes_of_current_triangle);
                
                node_coordinates_of_current_triangle = mesh_data.node_data. ...
                    coordinates(nodes_of_current_triangle, :);
                
                xe = node_coordinates_of_current_triangle(:, 1);
                ye = node_coordinates_of_current_triangle(:, 2);
                
                
                
                energy_integrant = @(zeta, eta) ...
                    problem_type.get_energy_integrant(...
                    finite_element_class, xe, ye, ...
                    node_potentials_of_current_triangle, zeta, eta, ...
                    vacuum_material * material_values(k), ...
                    vacuum_material * material_values(k));
                
                
                energy = energy + 1/2 * ...
                    GaussianQuadrature.integrate_2D_normalized_triangle_region( ...
                    energy_integrant, 7);
                
            end
            
            Misc.print_message(sprintf('Calculated energy: %d', energy))
            
        end
        
        function node_potentials = get_node_potentials_struct(mesh_data, unknowns, ...
                result_to_global_node_mapping)
            
            
            
            N_nodes = mesh_data.node_data.number_of_nodes;
            
            potentials = zeros(N_nodes, 1);
            
            potentials(mesh_data.dirichlet_boundary_data.IDs) = ...
                mesh_data.dirichlet_boundary_data.values;
            
            potentials(result_to_global_node_mapping) = unknowns;
            
            node_potentials = struct(...
                'potentials', potentials, ...
                'coordinates', mesh_data.node_data.coordinates);
        end
        
        function plot_node_potentials(node_potentials, mesh_plot)
            
            N_nodes = length(node_potentials.potentials);
            
            max_V = max(node_potentials.potentials);
            min_V = min(node_potentials.potentials);
            
            
            x_min = min(node_potentials.coordinates(:, 1));
            x_max = max(node_potentials.coordinates(:, 1));
            
            y_min = min(node_potentials.coordinates(:, 2));
            y_max = max(node_potentials.coordinates(:, 2));
            
            number_of_grid_lines = 1000;
            
            dx = (x_max - x_min) / number_of_grid_lines;
            dy = (y_max - y_min) / number_of_grid_lines;
            
            x = y_min : dx : y_max;
            y = y_min : dy : y_max;
            
            [X, Y] = meshgrid(x, y);
            
            Z = griddata(node_potentials.coordinates(:, 1), node_potentials.coordinates(:, 2), ...
                node_potentials.potentials, X, Y);
             
            %C = 1/(min_V - max_V) * Z;
            [a,b] = size(Z);
            Z_tmp = zeros(a,b);
            
            figure()
            surf(X, Y, Z, 'EdgeColor', 'none')
            colormap jet
            colorbar
            xlabel('Distance in m')
            ylabel('Distance in m')
            zlabel('Potential')
            title('3D view of potential')
            
            
            h = figure();
            
            Misc.CloneFig(mesh_plot, h);
            hold on
            surf(X, Y, Z_tmp, Z, 'FaceAlpha',0.75, 'EdgeColor', 'none');
            colormap jet
            c = colorbar;

%             tick_labels_num = min_V + (max_V - min_V) * c.Ticks;
%             tick_labels_str = sprintf('%f;',tick_labels_num);
%             tick_labels = strsplit(tick_labels_str, ';');
%             tick_labels = tick_labels(1:end-1);
%             c.TickLabels = tick_labels;
            
            title('Potential with mesh')
            axis equal
            xlabel('Distance in m')
            ylabel('Distance in m')

        end
        
        function plot_electric_field(mesh_data, node_potentials, problem_setup, ...
                mesh_plot)
            
            % Get prblem type class
            [problem_type, success] = Misc. ...
                get_problem_type_class_from_problem_type_number( ...
                problem_setup.problem_type);
            if ~success
                return
            end
            
            finite_element_class = ...
                Misc.get_finite_element_class_from_mesh_order(problem_setup.mesh_order);
            
            
            number_of_finite_elements = mesh_data.triangle_data.number_of_triangles;
            
            % Preallocate memory
            electric_field = zeros(number_of_finite_elements, 2);
            coordinates = zeros(number_of_finite_elements, 2);
            
            % Calculate field at center point of triangle
            zeta = 1/3;
            eta = 1/3;
            
            dN_dZeta = finite_element_class. ...
                get_zeta_shape_function_derivative_matrix(zeta, eta);
            
            dN_dEta = finite_element_class. ...
                get_eta_shape_function_derivative_matrix(zeta, eta);
            
            shape_functions = finite_element_class.get_shape_function_matrix(zeta, eta);
            
            for k = 1 : number_of_finite_elements
                
                nodes_of_current_triangle = mesh_data.triangle_data.nodes(k, :);
                node_potentials_of_current_triangle = node_potentials(...
                    nodes_of_current_triangle);
                
                node_coordinates_of_current_triangle = mesh_data.node_data. ...
                    coordinates(nodes_of_current_triangle, :);
                
                xe = node_coordinates_of_current_triangle(:, 1);
                ye = node_coordinates_of_current_triangle(:, 2);
                
                % x-components
                electric_field(k, 1) = ...
                    -(node_potentials_of_current_triangle' * dN_dZeta) * ...
                    (ye' * dN_dEta) + ...
                    (node_potentials_of_current_triangle' * dN_dEta) *...
                    (ye' * dN_dZeta);
                
                % y-component
                electric_field(k, 2) = ...
                    (node_potentials_of_current_triangle' * dN_dZeta) * ...
                    (xe' * dN_dEta) - ...
                    (node_potentials_of_current_triangle' * dN_dEta) *...
                    (xe' * dN_dZeta);
                
                coordinates(k, 1) = xe' * shape_functions;
                coordinates(k, 2) = ye' * shape_functions;

            end
            
            h = figure();
            Misc.CloneFig(mesh_plot, h)
            hold on
            quiver(coordinates(:, 1), coordinates(:, 2), electric_field(:,1), electric_field(:, 2))
            title('Electric field')
            
            h = figure();
            Misc.CloneFig(mesh_plot, h)
            hold on
            
            
            abs_electric_field = sqrt(electric_field(:, 1) .^2 + electric_field(:, 2) .^2);
            
            max_V = max(abs_electric_field);
            min_V = min(abs_electric_field);
            
            % Generate colors for plot:
            % Minimum value of the node potentials should have blue colour. 
            % Maximum value should have red colour.
            % 
            % Colour generation is done using HSV colour-space.
            % The 'Hue' value valculated in a way that nodes with the minimum value
            % equal to a hue of 240° (blue color, in Matlab: 8/9 as the value 1
            % corresponds to a hue of 360°), and the nodes with a maximum value equal
            % to a hue of 0 (red color).
            %
            % The transormation is done by a simple linear function determined by the
            % two points:
            % @ max_V: hue = 0
            % @ min_V: hue = 8/9
            %
            hsv_hue = 2/3 * 1/(min_V - max_V) * (abs_electric_field - max_V);
            
            hsv_data = [hsv_hue, ones(length(abs_electric_field), 2)];
            
            
            scatter(coordinates(:, 1), ...
                coordinates(:, 2), ...
                20, hsv2rgb(hsv_data), 'filled');
            title('Absolute value of electric field')
            
        end
        
        function mesh_plot = get_mesh_plot_data(mesh_data, problem_setup)
            
            number_of_plot_points_per_edge = 10;
            
            
            number_of_finite_elements = mesh_data.triangle_data.number_of_triangles;
            
            finite_element_class = ...
                Misc.get_finite_element_class_from_mesh_order(problem_setup.mesh_order);
            
            
            % Zeta and eta coordinates for the first edge (the one on the zeta axis)
            zeta_1 = 0 : 1/(number_of_plot_points_per_edge - 1) : 1;
            eta_1 = zeros(1, number_of_plot_points_per_edge);
            
            % Zeta and eta coordinates for the second edge (the one on the eta axis)
            zeta_2 = zeros(1, number_of_plot_points_per_edge);
            eta_2 = 0 : 1/(number_of_plot_points_per_edge - 1) : 1;
            
            % Zeta and eta coordinates for the thirs edge
            zeta_3 = 0 : 1/(number_of_plot_points_per_edge - 1) : 1;
            eta_3 = 1 - zeta_3;
            
            shape_functions_edge_1 = finite_element_class. ...
                get_shape_function_matrix(zeta_1, eta_1);
            
            shape_functions_edge_2 = finite_element_class. ...
                get_shape_function_matrix(zeta_2, eta_2);
            
            shape_functions_edge_3 = finite_element_class. ...
                get_shape_function_matrix(zeta_3, eta_3);
            
            mesh_plot = figure();
            
            for k = 1 : number_of_finite_elements
                
                nodes_of_current_triangle = mesh_data.triangle_data.nodes(k, :);
                
                node_coordinates_of_current_triangle = mesh_data.node_data. ...
                    coordinates(nodes_of_current_triangle, :);
                
                xe = node_coordinates_of_current_triangle(:, 1);
                ye = node_coordinates_of_current_triangle(:, 2);
                
                 
                
                x_edge_1 = xe' * shape_functions_edge_1;
                y_edge_1 = ye' * shape_functions_edge_1;
                
                x_edge_2 = xe' * shape_functions_edge_2;
                y_edge_2 = ye' * shape_functions_edge_2;
                
                x_edge_3 = xe' * shape_functions_edge_3;
                y_edge_3 = ye' * shape_functions_edge_3;
                
                plot(x_edge_1, y_edge_1, 'color', 'black')
                hold on
                plot(x_edge_2, y_edge_2, 'color', 'black')
                hold on
                plot(x_edge_3, y_edge_3, 'color', 'black')
                hold on         
            end
            
            hold off
            title('Mesh')
            axis equal
            xlabel('Distance in m')
            ylabel('Distance in m')
    
        end
        
        
    end
end

