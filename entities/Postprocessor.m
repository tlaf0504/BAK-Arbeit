classdef Postprocessor
    
    methods(Static)
        function plot_electrostatic_potential(problem_location)
            
            tmp = pwd;
            cd(problem_location)
            load('problem_setup.mat', 'problem_setup')
            load(fullfile(problem_location, 'internal', 'mesh_data.mat'), 'mesh_data');
            load(fullfile(problem_location, 'results', 'results.mat'), 'unknowns', ...
                'result_to_global_node_mapping');
            
            
            node_potentials = Postprocessor.get_node_potentials_struct(mesh_data, ...
                unknowns, ...
                result_to_global_node_mapping);
            
            Postprocessor.plot_node_potentials(node_potentials);
            
            
            Postprocessor.calculate_energy(mesh_data, node_potentials.potentials, ...
                problem_setup);
            
            Postprocessor.plot_electric_field(mesh_data, node_potentials.potentials, ...
                problem_setup)
            
            
            cd(tmp);
        end
        
        function postprocessor_test(problem_location)
            tmp = pwd;
            cd(problem_location)
            load('problem_setup.mat', 'problem_setup')
            load(fullfile(problem_location, 'internal', 'mesh_data.mat'), 'mesh_data');
            load(fullfile(problem_location, 'results', 'results.mat'), 'unknowns', ...
                'result_to_global_node_mapping');
            
            
            
            
            a = triangulation(double(mesh_data.triangle_data.nodes), ...
               mesh_data.node_data.coordinates);
            
%             a = delaunay(mesh_data.node_data.coordinates(:, 1), ...
%                 mesh_data.node_data.coordinates(:, 2));
%             triplot(a)
            
            %triplot(a)
            
            keyboard
            
            
            
            
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
        
        function plot_node_potentials(node_potentials)
            
            N_nodes = length(node_potentials.potentials);
            
            max_V = max(node_potentials.potentials);
            min_V = min(node_potentials.potentials);
            
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
            hsv_hue = 2/3 * 1/(min_V - max_V) * (node_potentials.potentials - max_V);
            
            hsv_data = [hsv_hue, ones(N_nodes, 2)];
            
            figure()
            scatter3(node_potentials.coordinates(:, 1), ...
                node_potentials.coordinates(:, 2), ...
                node_potentials.potentials, 20, hsv2rgb(hsv_data), 'filled');
            title('Potential V in nodes')
            
        end
        
        function plot_electric_field(mesh_data, node_potentials, problem_setup)
            
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
            
            figure()
            quiver(coordinates(:, 1), coordinates(:, 2), electric_field(:,1), electric_field(:, 2))
            keyboard
        end
        
    end
end

