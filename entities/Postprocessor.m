classdef Postprocessor
    
    methods(Static)
        function plot_electrostatic_potential(problem_location)
            tmp = pwd;
            cd(problem_location)
            load('problem_setup.mat', 'problem_setup')
            load(fullfile(problem_location, 'internal', 'mesh_data.mat'), 'mesh_data');
            load(fullfile(problem_location, 'results', 'results.mat'), 'unknowns', ...
                'result_to_global_node_mapping');
            
            N_nodes = length(mesh_data.node_data{1});
            
            %                         x                         y                       value
            node_plot_data = [mesh_data.node_data{2}(:, 1), mesh_data.node_data{2}(:, 2), zeros(N_nodes, 1)];
            
            node_plot_data(mesh_data.dirichlet_boundary_data{1}, 3) = ...
                mesh_data.dirichlet_boundary_data{2};
            
            node_plot_data(result_to_global_node_mapping, 3) = unknowns;
            
            max_V = max(node_plot_data(:, 3));
            min_V = min(node_plot_data(:, 3));
            
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
            hsv_hue = 2/3 * 1/(min_V - max_V) * (node_plot_data(:, 3) - max_V);
            
            hsv_data = [hsv_hue, ones(N_nodes, 2)];
            
            figure()
            scatter3(node_plot_data(:, 1), node_plot_data(:, 2), node_plot_data(:, 3), 20, hsv2rgb(hsv_data), 'filled');
            title('Potential V in nodes')
            cd(tmp);
        end
    end
end

