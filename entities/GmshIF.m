classdef GmshIF
    
    properties(Constant)
        
        supported_mesh_file_versions = [1,2,4];
        
        % Format string for specific mesh file version.
        % See http://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options
        % command line option '-format'
        mesh_file_version_formats = containers.Map(...
            GmshIF.supported_mesh_file_versions, {'msh1', 'msh22', 'msh4'});
    end
    
    methods(Static)
        function success = mesh(problem_location)
            
            % Plot section separator to console
            Misc.print_message(sprintf('%s\n%s\n%s', ...
                Misc.console_section_separator, ...
                'Mesher', Misc.console_section_separator));
            
            
            
            success = 1;
            tmp = pwd;
            cd(problem_location);
            
            load('problem_setup.mat', 'problem_setup');
            
            if problem_setup.state < Setup.setup_state_initializing_finished
                msg = sprintf(['Error. Problem initilization were not ', ...
                    'finished. Maybe the setup procedure was interrupted.', ...
                    'Please rerun the complete setup procedure.']);
                Misc.print_error_message(msg)
                success = 0;
                return
            end
            
            Misc.print_message(sprintf('Meshing file %s...\n', ...
                problem_setup.geometry_file));
            tic
            
            % Check if mesh order is valid
            if ~ismember(problem_setup.mesh_order, Misc.possible_mesh_orders)
                msg = sprintf(['Error. Specified mesh_order %d is not one ', ...
                    'of the following allowed mesh orders:\n%s'], ...
                    problem_setup.mesh_order, ...
                    sprintf('%s ', Misc.possible_mesh_orders));
                Misc.print_error_message(msg);
                success = 0;
                return
            end
            
            % Mesh geometry file
            % Use mesh file version 2
            [success, output_filename] = GmshIF.gmshMesh2D(problem_setup, 2);
            if ~success
                return
            end
            
            problem_setup.mesh_file = output_filename;
            problem_setup.state = Setup.setup_state_meshing_finished;
            Setup.update_problem_setup_file(problem_setup);
            
            time = toc;
            Misc.print_message(sprintf('Done\nElapsed time is %d seconds.\n', time));
            Misc.print_message(Misc.console_section_separator);
            Misc.print_message('\n');
            
            cd(tmp);
            
        end
        
        function [success, output_filename] = gmshMesh2D(problem_setup, ...
                mesh_file_version)
            
            success = 1;
            output_filename = [];
            global gmsh_exec;
            
            if ~ismember(mesh_file_version, GmshIF.supported_mesh_file_versions)
                Misc.print_error_message(['Mesh file version %d not ', ...
                    'supported. Followind versions are supported:\n%s'], ...
                    mesh_file_version, sprintf('%d ', ...
                    GmshIF.supported_mesh_file_versions));
                success = 0;
                return
            end
            
            output_filename = sprintf('%s.msh', problem_setup.problem_name);
            mesh_file_format = GmshIF.mesh_file_version_formats(mesh_file_version);
            
            command = sprintf(['%s -2 -o %s -format %s -order %d %s'], ...
                gmsh_exec, ...
                output_filename, ...
                mesh_file_format, ...
                problem_setup.mesh_order, ...
                problem_setup.geometry_file);
            
            [state, msg] = system(command);
            if state ~= 0
                msg = sprintf(['Error when executing command "%s".\nSystem ', ...
                    'returned following message:\n\n%s'], command, msg);
                Misc.print_error_message(msg);
                success = 0;
                return
            end
            
            GmshIF.plot_gmsh_log(msg);
        end
        
        function plot_gmsh_log(log)
            fprintf(['Gmsh output:\n', ...
                '------------------------------------------------------------\n', ...
                '%s\n', ...
                '------------------------------------------------------------\n'], ...
                log);
        end
    end
end