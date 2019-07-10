classdef GmshIF
    
    properties(Constant)
        
        % Gmsh defines multiple mesh-file versions. 
        % Check http://gmsh.info/doc/texinfo/gmsh.html#File-formats for supported
        % mesh file versions.
        % 
        % Following versions are supported by the mesher:
        % Version 1 <-> 1
        % Version 2.2 <-> 2
        % Version 4 <-> 4
        % 
        % Attention: Not all mesh file versiond spcified here are supported by the
        % parser.
        supported_mesh_file_versions = [1,2,4];
        
        % Format string for specific mesh file version.
        % See http://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options
        % command line option '-format'
        mesh_file_version_formats = containers.Map(...
            GmshIF.supported_mesh_file_versions, {'msh1', 'msh22', 'msh4'});
    end
    
    methods(Static)
        function success = mesh(problem_location)
            
        % function success = mesh(problem_location)
        %
        % DESCRIPTION: Main function for meshing the problem.
        %
        %     The meshing procedure calls Gmsh via commandline. The final mesh file
        %     is stored to the problem directory using the following filename:
        %         <problem_setup.problem_name>.msh
        %
        %     Currently the parser only supports mesh file version 2.2
        %
        % INPUT:
        %     problem_location.....Directory of the problem
        %
        % OUTPUT:
        %     success.....If 1 is returned, the procedure was successful. If 0 was
        %         returned an error occurred.
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            % Plot section separator to console
            Misc.print_message(sprintf('%s\n%s\n%s', ...
                Misc.console_section_separator, ...
                'Mesher', Misc.console_section_separator));
            
            
            % Define return parameter
            success = 1;
            
             % Change working directory to problem directory
            tmp = pwd; % Backup current working directory
            cd(problem_location);
            
            load('problem_setup.mat', 'problem_setup');
            
            % Check input parameters
            success = GmshIF.check_input_parameters(problem_setup);
            if ~success
                return
            end
            
            
            Misc.print_message(sprintf('Meshing file %s...\n', ...
                problem_setup.geometry_file));
            tic

            
            % Mesh geometry file
            [success, output_filename] = GmshIF.gmshMesh2D(problem_setup);
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
        
        function [success, output_filename] = gmshMesh2D(problem_setup)
            
        % function [success, output_filename] = gmshMesh2D(problem_setup, ...
        %        mesh_file_version)
        %
        % DESCRIPTION: Function for creating a 2D mesh with Gmsh
        %
        %     The meshing procedure calls Gmsh via commandline. The final mesh file
        %     is stored to the problem directory using the following filename:
        %         <problem_setup.problem_name>.msh
        %
        % INPUT:
        %     problem_setup.....Structure containing all neccesary information about
        %         problem.
        %
        % OUTPUT:
        %     success.....If 1 is returned, the procedure was successful. If 0 is
        %         returned, an error occurred.
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            success = 1;
            global gmsh_exec;
            
            mesh_file_version = problem_setup.mesh_file_version;
            
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
            
            GmshIF.print_gmsh_log(msg);
        end
        
        function print_gmsh_log(log)
            
        % function plot_gmsh_log(log)
        %
        % DESCRIPTION: Function for printing the log returned from Gmsh to the matlab
        %     console.
        %
        % INPUT:
        %     log....Log returned from Gmsh
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            fprintf(['Gmsh output:\n', ...
                '------------------------------------------------------------\n', ...
                '%s\n', ...
                '------------------------------------------------------------\n'], ...
                log);
        end
        
        function success = check_input_parameters(problem_setup)
            
        % function success = check_input_parameters(problem_setup)
        %
        % DESCRIPTION: Function for checking the input parameters of the procedure.
        %
        %
        %     Following parameters are evaluated:
        %         -) problem_setup.state: Evaluating of the initialization procedure
        %             was completed.
        %         -) problem_setup.mesh_order: Evaluating if the user specified a
        %             valid mesh order
        %         -) prolem_setup.mesh_file_version: Evaluating if the user specified
        %             a valid mesh file verison.
        %
        %     The meshing procedure calls Gmsh via commandline. The final mesh file
        %     is stored to the problem directory using the following filename:
        %         <problem_setup.problem_name>.msh
        %
        % INPUT:
        %     problem_setup.....Structure containing all neccesary information about
        %         problem.
        %
        % OUTPUT:
        %     success.....If 1 is returned, the procedure was successful. If 0 is
        %         returned, an error occurred.
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            % Check if inititalization procedure was finished successfully
            success = GmshIF.check_setup_state(problem_setup.state);
            if ~success
                return
            end
            
            % Check if the specified mesh order is valid
            success = GmshIF.check_mesh_order(problem_setup.mesh_order);
            if ~success
                return
            end
            
            % Check if the specified mesh file version is valid
            success = GmshIF.check_mesh_file_version(problem_setup.mesh_file_version);
            if ~success
                return
            end
        end
        
        function success = check_setup_state(setup_state)
            
        % function success = check_setup_state(setup_state)
        %
        % DESCRIPTION: Function for evaluating the state of the solution process.
        %
        %     The meshing procedure requires a completed initialization procedure.
        %     This function uses the problem_setup.state parmeter to determine this.
        %
        % INPUT:
        %     setup_state.....Content of problem_setup.state
        %
        % OUTPUT:
        %     success.....If 1 is returned, the initialization procedure was
        %         completed. If 0 is returned it was interrupted.
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            success = 1;
            
            if setup_state < Setup.setup_state_initializing_finished
                msg = sprintf(['Error. Problem initilization were not ', ...
                    'finished. Maybe the setup procedure was interrupted.', ...
                    'Please rerun the complete setup procedure.']);
                Misc.print_error_message(msg)
                success = 0;
                return
            end
            
        end
        
        function success = check_mesh_file_version(mesh_file_version)
            
        % function success = check_mesh_file_version(mesh_file_version)
        %
        % DESCRIPTION: Function for evaluating if the user specified a valid mesh
        %     file version.
        %
        %     GmshIF.supported_mesh_file_versions specifies the supported mesh file
        %     versions.
        %
        % INPUT:
        %     mesh_file_version.....Integer number specifying the version of the
        %         output mesh file.
        %
        % OUTPUT:
        %     success.....If 1 is returned, the input is correct. If 0 is returned,
        %         it is incorrect.
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            success = 1;
            
            if ~ismember(mesh_file_version, GmshIF.supported_mesh_file_versions)
                msg = sprintf(['Specified mesh file version %d is not supported. ', ...
                    'Supported values for mesh_file_version are: %s'], ...
                    mesh_file_version, ...
                    sprintf('%d ', GmshIF.supported_mesh_file_versions));
                
                Misc.print_error_message(msg);
                success = 0;
                return
            end
        end
        
        function success = check_mesh_order(mesh_order)
            
        % function success = check_mesh_order(mesh_order)
        %
        % DESCRIPTION: Function for evaluating if the user specified a valid mesh
        %     order.
        %
        %     Misc.supported_mesh_orders specifies the supported mesh order.
        %
        % INPUT:
        %     mesh_order.....Integer number specifying the order of the FEM mesh.
        %
        % OUTPUT:
        %     success.....If 1 is returned, the input is correct. If 0 is returned,
        %         it is incorrect.
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            
            success = 1;
            
            if ~ismember(mesh_order, Misc.supported_mesh_orders)
                msg = sprintf(['Error. Specified mesh_order %d is not one ', ...
                    'of the following allowed mesh orders:\n%s'], ...
                    problem_setup.mesh_order, ...
                    sprintf('%s ', Misc.possible_mesh_orders));
                Misc.print_error_message(msg);
                success = 0;
                return
            end    
        end
    end
end