
% Setup.m
%
% DESCRIPTION: File containing the class "Setup" which implements the setup procedure
%     for a new or existing problem.
%
%
%
%
% Author: Tobias Lafer
% e-mail: lafer@student.tugraz.at
%
% Created on: 26.06.2019
% Last modified 26.06.2019 by Tobias Lafer


classdef Setup 
% classdef Setup
%
% DESCRIPTION: Class implementing the setup procedure for a new or existing problem
%
%    All methods are implemented statically, so no instantiation is required.
%    User should only use one of the two following methods:
%        -) Setup.new_setup(): Sets up a new problem
%        -) Setup.existing_setup(): Loads an existing problem from file 
%               "problem_setup.mat" in problem directory.
%
%    DEPENDENCIES:
%        -) Class 'Misc' for message printing functions and file existence checking
%           functions
%       
%
% Author: Tobias Lafer
% e-mail: lafer@student.tugraz.at
%
% Created on: 26.06.2019
% Last modified 26.06.2019 by Tobias Lafer
    
    properties(Constant)
        % ===== Setup states
        % The procedure of solving a problem can be divided into different steps.
        % Each step gets asigned to a 'state' which is stored in the problem_setup
        % structure helping to determine which step to apply when loading an existing
        % setup.
        % 
        % The procedure is mainly divided into 5 steps. Each step has its own class
        % containing all step-specific methods and properties. 
        %
        % The five steps are:
        % -) Preprocessing and Initializing, implemented in this file
        % -) Meshing, done by the external software package 'Gmsh'
        %    (http://gmsh.info/). Interface to this software is implemented in class
        %    'GmshIF'
        % -) Parsing mesh- and settings file, implemented in class 'Parser'
        % -) Solving the problem, implemented in class 'Solver'
        % -) Postprocessing (Figure generation, energy calculation etc.), implemented
        %    in class 'Postprocessor'
        %
        % In case of a new setup, the problem_setup structure is created with 
        % problem_setup.state = setup_state_initializing, indicating that the 
        % initialozation process was not finished yet. After initialization,
        % problem_setup.state gets assigned to setup_state_initializing_finished
        % indicating that the problem is ready for meshing and so on and so forth.
        %
        
        setup_state_initializing = 0; % Initialization not finished
        setup_state_initializing_finished = 1; % Initialization finished
        setup_state_meshing_finished = 2; % Meshing finished
        setup_state_parsing_finished = 3; % Parsing finished
        setup_state_solving_finished = 4; % Solving finished
        setup_state_postprocessing_finished = 5; % Postprocessing finished
    end
    methods(Static)
        
        
        function success = new_setup(problem_location, ...
                geometry_filename, settings_filename, mesh_order, problem_type, ...
                mesh_file_version, load_pregenerated_mesh, mesh_filename)
            
        % function success = new_setup(problem_location, ...
        %         geometry_filename, settings_filename, mesh_order, problem_type, ...
        %         mesh_file_version)
        %
        % DESCRIPTION: Function applying setup procedure for a new problem.
        % 
        %    This function should be called in case of a new of modified geometry
        %    file.
        %
        %    If the setup procedure was successful, a file 'problem_setup.mat' is
        %    created in the problem directory containing the 'problem_setup'
        %    structure.
        %
        % INPUT:
        %     problem_location.....Path to the directory containing the .geo file.
        %         An own directory for each problem is recommended because the
        %         program will create subdirectories and additional files.
        %     geometry_filename.....Name of the geometry file WITH file extension.
        %         e.g. 'capacitor.geo'
        %     settings_filenam.....Name of the settings file WITH extension.
        %         e.g. 'settings.set'
        %     mesh_order.....Integer value containing the order of the triangluar
        %         mesh created by the program. Following values are supported:
        %             -) 1: First order (linear) mesh
        %             -) 2: Second order (quadratic) mesh
        %             -) 3: Third order (cubic) mesh
        %     problem_type.....String representing the type of the problem to be
        %         solved. See Misc.supported_problem_types for further information
        %         about the supported problem types.
        %     mesh_file_version.....Format of the generate mesh file. See 
        %         GmshIF.supported_mesh_file_versions for further information.
        %     load_pregenerated_mesh.....Set this flag to 1 to load a pregenerated
        %         mesh. Settings this flag to 0 enables automatic mesh generation by 
        %         the software.
        %     mesh_filename.....Optional argument. Specifies the mesh file to be 
        %         loaded when parameter 'load_pregenerated_mesh' is set to 1.
        %
        %
        % OUTPUT:
        %     success.....Logical value indicating the success or failure of the
        %         setup process.
        %
        % CREATED FILES:
        %   problem_setup.mat.....MAT-File containing the 'problem_setup' structure,
        %       which stores all neccesary information about the problem. 
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 26.06.2019 by Tobias Lafer
        %
        
            
            success = 1;
            
            % Backup current working directory and go to problem directory
            tmp = pwd; 
            cd(problem_location)
            
            % Plot section separator to console
            Misc.print_message(sprintf('%s\n%s\n%s', ...
                Misc.console_section_separator, ...
                'Setup', Misc.console_section_separator));
            
            Misc.print_message('Setting up new problem:')
            
            if load_pregenerated_mesh
                msg = sprintf(['\tProblem location: %s\n\tGeometry-file: %s', ...
                    '\n\tSettings-file: %s\n\tUsing pregenerated mesh-file: %s'], ...
                    problem_location, geometry_filename, ...
                    settings_filename, mesh_filename);
            else
                
                msg = sprintf(['\tProblem location: %s\n\tGeometry-file: %s', ...
                    '\n\tSettings-file: %s\n\tMeshig is done by software'], ...
                    problem_location, geometry_filename, ...
                    settings_filename);
            end
            
           
            Misc.print_message(msg);
            
            
            % Check existence of geometry and settings file
            if  ~Setup.check_geometry_and_settings_file_existence(problem_location, ...
                    geometry_filename, settings_filename)
                success = 0;
                return
            end
            
            % Check if pregenerated mesh file exists
            if load_pregenerated_mesh
                if ~Misc.check_file_existence(problem_location, mesh_filename)
                    msg = sprintf(['Mesh-file "%s" in folder "%s" not ', ...
                        'found'], mesh_filename, problem_location);
                    Misc.print_error_message(msg);
                    success = 0;
                    return
                end
            end
            
            % Check if user specified a valid problem type
            success = Setup.evaluate_problem_type(problem_type);
            if ~success
                return
            end

            % Setting up folders for results, plots etc.
            Setup.setup_solution_hierachy();
            
            % Create new structure for storing problem information
            problem_setup = Setup.new_problem_setup(problem_location, ...
                geometry_filename, settings_filename, mesh_order, problem_type, ...
                mesh_file_version, load_pregenerated_mesh, mesh_filename);
            
            if load_pregenerated_mesh
                problem_setup.state = Setup.setup_state_meshing_finished;
            else
                problem_setup.state = Setup.setup_state_initializing_finished;
            end
            
            Misc.print_message(Misc.console_section_separator);
            Misc.print_message('\n');
            
            % Create .mat file
            save('problem_setup.mat', 'problem_setup')
            
            % Go back to previous working directory
            cd(tmp)
            
        end
        
        
        
        
        
        
        function [success, next_computation_step] = load_setup(problem_location)
            
            success = 1;
            next_computation_step = 0;
            
            % Backup current working directory and go to problem directory
            tmp = pwd; 
            cd(problem_location)
            
            % Plot section separator to console
            Misc.print_message(sprintf('%s\n%s\n%s', ...
                Misc.console_section_separator, ...
                'Setup', Misc.console_section_separator));
            
            Misc.print_message('Loading existing problem:\n')
            msg = sprintf('\tProblem location: "%s"', problem_location);
            Misc.print_message(msg);
            
            load('problem_setup.mat', 'problem_setup');
            
            

            % Go through each computation step until the current state of the proble
            % is reached. (The value stored in problem_setup.state).
            % In each step, the required files etc. are evaluated.
            
            checking_setup_completed = false;
            
            % The first possible step is the one after the initialization is
            % completed. (Otherwise no setup-file ("problem_setup.mat") would exist).
            current_setup_state = Setup.setup_state_initializing_finished;
            
            while ~checking_setup_completed
                switch current_setup_state
                    
                    
                    case Setup.setup_state_initializing_finished
                        [success, next_loading_step, loading_completed, ...
                            next_computation_step] = ...
                            Setup.loading_setup_check_initialization_step(...
                            problem_setup, current_setup_state);
                        if ~success
                            return
                        end
                        
                        if loading_completed
                            break;
                        else
                            current_setup_state = next_loading_step;
                        end

      
                    case Setup.setup_state_meshing_finished
                        [success, next_loading_step, loading_completed, ...
                            next_computation_step] = ...
                            Setup.loading_setup_check_meshing_step( ...
                            problem_setup, current_setup_state);
                        
                        if ~success
                            return
                        end
                        
                        if loading_completed
                            break;
                        else
                            current_setup_state = next_loading_step;
                        end
                        
                        
                    case Setup.setup_state_parsing_finished
                        [success, next_loading_step, loading_completed, next_computation_step] = ...
                            Setup.loading_setup_check_parsing_step( ...
                            problem_setup, current_setup_state);
                        
                        if ~success
                            return
                        end
                        
                        if loading_completed
                            break;
                        else
                            current_setup_state = next_loading_step;
                        end
                        
                        
                        
                        
                   case Setup.setup_state_solving_finished
                       [success, next_loading_step, loading_completed, next_computation_step] = ...
                           Setup.loading_setup_check_solving_step( ...
                           problem_setup, current_setup_state);
                       
                       if ~success
                            return
                        end
                        
                        if loading_completed
                            break;
                        else
                            current_setup_state = next_loading_step;
                        end
                      
                        
                    case Setup.setup_state_postprocessing_finished
                        [success, next_loading_step, loading_completed, next_computation_step] = ...
                            Setup.loading_setup_check_postprocessing_step( ...
                            problem_setup, current_setup_state);
                        
                        if ~success
                            return
                        end
                        
                        if loading_completed
                            break;
                        else
                            current_setup_state = next_loading_step;
                        end
                        
                        
                end
            end
            
            if ~success
                Misc.print_error_message(...
                    ['Loading setup failed. Please rerun the complete ', ...
                    'setup.']);
            else
                Misc.print_message('\nLoading completed.')
            end
            
             cd(tmp)
             Misc.print_message(Misc.console_section_separator);
             Misc.print_message('\n');

        end
        
        function success = check_geometry_and_settings_file_existence( ...
                problem_location, geometry_filename, settings_filename)
            
        % function success = check_geometry_and_settings_file_existence( ...
        %        problem_location, geometry_filename, settings_filename)
        %
        % DESCRIPTION: Checks the existence of the getometry and settings files in
        %     the problem directory
        %
        % INPUT:
        %     problem_location.....Path to the directory containing the .geo file.
        %         An own directory for each problem is recommended because the
        %         program will create subdirectories and additional files.
        %     geometry_filename.....Name of the geometry file WITH file extension.
        %         e.g. 'capacitor.geo'
        %     settings_filenam.....Name of the settings file WITH extension.
        %         e.g. 'settings.set'
        %
        % OUTPUT:
        %     success.....Logical value indicating if the two files exist or not.
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 26.06.2019 by Tobias Lafer
        %
            
            success = 1;
            % Check for specifed geometry file in problem folder
            if ~Misc.check_file_existence(problem_location, geometry_filename)
                msg = sprintf(['Geometry file "%s" in folder "%s" not ', ...
                    'found'], geometry_filename,  problem_location);
                Misc.print_error_message(msg);
                success = 0;
                return
            end
            
            % Check for settings geometry file in problem folder
            if ~Misc.check_file_existence(problem_location, settings_filename)
                msg = sprintf(['Settings file "%s" in folder "%s" not ', ...
                    'found'], settings_filename, problem_location);
                Misc.print_error_message(msg);
                success = 0;
                return
            end
            
        end
        
        function setup_solution_hierachy()
            
            foldername = 'internal';
            Setup.create_or_cleanup_folder(foldername);
            
            foldername = 'results';
            Setup.create_or_cleanup_folder(foldername);
            
            foldername = 'plots';
            Setup.create_or_cleanup_folder(foldername);
        end
        
        
        function create_or_cleanup_folder(folder)
            
        % function create_or_cleanup_folder(folder)
        %
        % DESCRIPTION: Function for either creating or cleaning up (deleting all
        %     contained files) from a folder.
        %
        % INPUT:
        %     folder.....Folder to be created or cleaned up 
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 03.07.2019 by Tobias Lafer
        %
            if exist(folder, 'dir')
                rmdir(folder, 's');
            end
            
            mkdir(folder)
        end
        
        
        function problem_setup = new_problem_setup(location, ...
                geometry_filename, settings_filename, mesh_order, problem_type, ...
                mesh_file_version, load_pregenerated_mesh, mesh_filename)
            
        % function problem_setup = new_problem_setup(location, ...
        %        geometry_filename, settings_filename, mesh_order, problem_type, ...
        %         mesh_file_version)
        %
        % DESCRIPTION: Function for creating a new strucutre containing all
        %     information about the problem.
        %
        %     The name of the structure is 'problem_setup' and is kept over the
        %     whole program.
        %     
        %     This structure is stored to a file named 'problem_setup.mat' in the 
        %     'internals' folder of the problem directory.
        %
        %     All members are created within this function and are updated along the
        %     program.
        %
        %     INPUT:
        %         location.....Path to the problem-directory
        %         geometry_filenams.....Name of the geometry file (including
        %             extension)
        %         settings_filename.....Name of the settings file (including 
        %             extension)
        %         mesh_order.....Integer number representing the order of the 
        %             triangluar mesh. (The order of the shape function polynomials)
        %             Must be either 1, 2 or 3.
        %             Attention: this parameter might change in future releases when
        %             more mesh types (square mesh etc.) or higher orders are
        %             supported.
        %        problem_type.....String specifying the type of the problem to be
        %            solved.
        %            Every problem type is assigned to an integer which is stored in
        %            the structure and used in the further problem for problem-type
        %            dependent decisions. (e.g. how to calculate the element
        %            matrices)
        %            See Misc.supported_problem_types for supported problem types.
        %        mesh_file_version.....Format of the generate mesh file. See 
        %            GmshIF.supported_mesh_file_versions for further information.
        %        load_pregenerated_mesh.....Set this flag to 1 to load a pregenerated
        %            mesh. Settings this flag to 0 enables automatic mesh generation by 
        %            the software.
        %        mesh_filename.....Optional argument. Specifies the mesh file to be 
        %            loaded when parameter 'load_pregenerated_mesh' is set to 1.
        %
        % OUTPUT:
        %     problem_setup.....Structuce containing all necessary information about
        %     the problem.
        %     It contains following members:
        % 
        %        -) problem_name.....Name of the problem, derived from the geometry
        %            filename. E.g. 'cylinder_cyp.geo' leads to the problem name 
        %            'cylinder_cap'
        %        -) problem_location.....Directory of the problem
        %        -) problem_type.....Integer number specifying the problem type.
        %            Integer is derived from the input parameter problem_type using
        %            Misc.supported_problem_types
        %        -) geometry_file.....Name of the geometry file. E.g. 
        %            'cylinder_cap.geo'
        %        -) settings_file.....Name of the settings file. E.g. 'settings.set'
        %        -) geometry_file_hash.....MD5 hash of geometry file.
        %            Used to check if geometry file has changed when loading an
        %            existing solution.
        %        -) settings_file_hash.....Same as for geometry file but for
        %            settings file.
        %        -) mesh_file_hash.....Same as for the geometry file but for the mesh
        %            file.
        %        -) mesh_data_file_hash.....Same as for geometry file but for mesh
        %            data file (.mat file containing the parsed data from the .msh 
        %            file)
        %        -) result_file_hash.....Same as for the geometry file but for the
        %            result file.
        %        -) mesh_order.....Order of the triangular mesh
        %        -) state.....State of the solution process. See property group
        %            "Setup states" in properties(Constant) section of this class.
        %        -) mesh_file.....Name of the mesh file generated later in this
        %            program using Gmsh.
        %        -) mesh_data_file.....Name of the .mat file containing the parsed
        %            information from the generated mesh file.
        %            This file is stored in the 'internals' subfolder of the problem
        %            directory.
        %        -) result_data_file.....Name of the .mat file containing the results
        %            for the problem. This file is stored in the 'results' subfolder
        %            of the problem directory.
        %        -) mesh_file_version.....Format of the generate mesh file. See 
        %            GmshIF.supported_mesh_file_versions for further information.
        %        -) load_pregenerated_mesh.....Flag indicating wheter a pregernated
        %            mesh should be used (1) or if the mesh should be generated by 
        %            the software (0)
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 04.08.2019 by Tobias Lafer
        %
        
            problem_setup = struct();
        
            [~, name, ~] = fileparts(geometry_filename);
            
            problem_setup.problem_name = name;
            problem_setup.problem_location = location;
            problem_setup.problem_type = Misc.supported_problem_types(problem_type);
            problem_setup.geometry_file = geometry_filename;
            problem_setup.settings_file = settings_filename;
            problem_setup.geometry_file_hash = DataHash(geometry_filename);
            problem_setup.settings_file_hash = DataHash(settings_filename);
            problem_setup.mesh_order = mesh_order;
            problem_setup.state = Setup.setup_state_initializing;
            problem_setup.mesh_file_version = mesh_file_version;
            problem_setup.load_pregenerated_mesh = load_pregenerated_mesh;
            
            if load_pregenerated_mesh
                problem_setup.mesh_file = mesh_filename;
                problem_setup.mesh_file_hash = DataHash(mesh_filename);
            else
                problem_setup.mesh_file = '';
                problem_setup.mesh_file_hash = '';
            end
            
            % Placeholders
            problem_setup.mesh_data_file = '';
            problem_setup.result_file = '';
       
            problem_setup.mesh_data_file_hash = '';
            problem_setup.result_file_hash = '';
        end
        
        
        function update_problem_setup_file(problem_setup)
            
        % function update_problem_setup_file(problem_setup)
        %
        % DESCRIPTION: Function for updating the problem_setup.mat file.
        %
        %
        % INPUT:
        %     problem_setup.....Structuce containing all necessary information about
        %     the problem.
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            save(fullfile(problem_setup.problem_location, 'problem_setup.mat'), ...
                'problem_setup');
        end
        
        function success = evaluate_problem_type(problem_type)
            
        % function success = evaluate_problem_type(problem_type)
        %
        % DESCRIPTION: Function evaluating the problem_type string specified by the
        %     user.
        %
        %     Misc.supported_problem_types defines the problem types supported and
        %     their corresponding ID. This function checks if the parameter
        %     problem_type in the main script war correctly set with one of the
        %     suppported problems.
        %
        %
        % INPUT:
        %     problem_type.....String specifying the type of the problem. See 
        %         Misc.supported_problem_types for all supported problem.
        % OUTPUT:
        %     success.....If function returns 1, the content of problem_type is
        %     valid. If 0 is returned, it is invalid.
        %
        %
        % Author: Tobias Lafer
        % e-mail: lafer@student.tugraz.at
        %
        % Created on: 26.06.2019
        % Last modified 05.07.2019 by Tobias Lafer
        %
            success = 1;
            
            supported_problems_strings = keys(Misc.supported_problem_types);
            tmp = '';
            
            if ~ismember(problem_type, supported_problems_strings)
                
                for k = 1 : length(supported_problems_strings)
                    tmp = [tmp, '"', supported_problems_strings{k}, '" '];
                end
                
                msg = sprintf(['Problem type "%s" is not supported.', ...
                    ' Supported problems:\n\t%s'], problem_type, tmp);
                Misc.print_error_message(msg);
                
                success = 0;
                return
            end
        end
        
        function success = check_file_hash(file, old_hash)
            
            success = 1;
            hash = DataHash(file);
            
            if hash ~= old_hash
                msg = sprintf(['File "%s" has changed since the ', ...
                    'last run. Please rerun the complete setup. \n', ...
                    'Old hash was %s, new hash is %s'], ...
                    file, hash, old_hash);
                Misc.print_message(msg);
                success = 0;
                return
                
            end
        end
        
        function [success, next_loading_step, loading_completed, next_computation_step] = ...
                loading_setup_check_initialization_step(...
                problem_setup, current_setup_state)
            
            success = 1;
            loading_completed = 0;
            next_loading_step = 0;
            next_computation_step = 0;
           
            
            problem_location = problem_setup.problem_location;
            
            % Check existence of geometry and settings file
            if  ~Setup.check_geometry_and_settings_file_existence(problem_location, ...
                    problem_setup.geometry_file, problem_setup.settings_file)
                success = 0;
                return
            end
            
            
            % Check if geometry file has not changed since the last run
            if ~Setup.check_file_hash(problem_setup.geometry_file, ...
                    problem_setup.geometry_file_hash)
                success = 0;
                return
            end
            
            % Check if settings file has not changed since the last run
            if ~Setup.check_file_hash(problem_setup.settings_file, ...
                    problem_setup.settings_file_hash)
                success = 0;
                return
            end
            
            msg = sprintf('\tUsing geometry file "%s"', problem_setup.geometry_file);
            Misc.print_message(msg);
            msg = sprintf('\tUsing settings file "%s"', problem_setup.settings_file);
            Misc.print_message(msg);
            
            if problem_setup.load_pregenerated_mesh
                
                if ~Misc.check_file_existence(problem_location, ...
                        problem_setup.mesh_file)
                    
                    msg = sprintf(['Mesh-file "%s" in folder "%s" not ', ...
                        'found'], problem_setup.mesh_file, problem_location);
                    Misc.print_error_message(msg);
                    success = 0;
                    return
                end
                
                % Check if pregenerated mesh-file has not changed
                if ~Setup.check_file_hash(problem_setup.mesh_file, ...
                        problem_setup.mesh_file_hash)
                    success = 0;
                    return
                end
                
                
                msg = sprintf('\tUsing mesh file "%s"', problem_setup.mesh_file);
                Misc.print_message(msg);
            end
            
            
            % If the previous run was terminated after the
            % initialization, continue the computation procedure with
            % the meshing of parsing procdure. (Depending on wheter a
            % pregenerated mesh should be used or not)
            if current_setup_state == problem_setup.state
                if problem_setup.load_pregenerated_mesh
                    next_computation_step = Setup.setup_state_meshing_finished;
                    Misc.print_message(['\nUsing pregenerated mesh-file. ', ...
                        'Continuing computation with meshing ', ...
                        ' procedure.']);
                else
                    next_computation_step = Setup.setup_state_initializing_finished;
                    Misc.print_message(['\nContinuing computation with mesh-parsing ', ...
                        ' procedure.']);
                end
                
                loading_completed = true;
                
                % Decide wheter to continue the loading procedure
            else
                
                % If a pregenerated mesh is used, evaluation of the
                % meshing procedure can be skipped.
                if problem_setup.load_pregenerated_mesh
                    next_loading_step = Setup.setup_state_parsing_finished;
                else
                    next_loading_step = Setup.setup_state_meshing_finished;
                end
            end
        end
        
        function [success, next_loading_step, loading_completed, next_computation_step] = ...
                loading_setup_check_meshing_step( ...
                problem_setup, current_setup_state)
            
            success = 1;
            loading_completed = 0;
            next_loading_step = 0;
            next_computation_step = 0;
           
            problem_location = problem_setup.problem_location;
            
            
            if ~Misc.check_file_existence(problem_location, ...
                    problem_setup.mesh_file)
                success = 0;
                
                msg = sprintf(['Mesh-file "%s" in folder "%s" not ', ...
                    'found'], problem_setup.mesh_file, problem_location);
                Misc.print_error_message(msg);
                return
            end
            
            % Check if hashes are the same
            if ~Setup.check_file_hash(problem_setup.mesh_file, ...
                    problem_setup.mesh_file_hash)
                success = 0;
                return
            end
            
            msg = sprintf('\tUsing mesh file "%s"', problem_setup.mesh_file);
            Misc.print_message(msg);
            
            % If the previous run was terminated after the
            % meshing step, continue the computation procedure with
            % parsing
            if current_setup_state == problem_setup.state
                
                next_computation_step = Setup.setup_state_meshing_finished;
                loading_completed = true;
                
                Misc.print_message(['\nContinuing computation with mesh-parsing ', ...
                    ' procedure.']);
                
                % Decide wheter to continue the loading procedure
            else
                next_loading_step = Setup.setup_state_parsing_finished;
            end
            
        end
        
        
        function [success, next_loading_step, loading_completed, next_computation_step] = ...
                loading_setup_check_parsing_step( ...
                problem_setup, current_setup_state)
            
            success = 1;
            loading_completed = 0;
            next_loading_step = 0;
            next_computation_step = 0;
           
            problem_location = problem_setup.problem_location;
            
            result_file = ...
                fullfile('internal', problem_setup.mesh_data_file);
            
            % Check if mesh data file exists
            if ~Misc.check_file_existence(problem_location, ...
                    result_file )
                
                msg = sprintf(['Mesh data file "%s" in folder "%s" not ', ...
                    'found'], problem_setup.mesh_data_file, ...
                    fullfile(problem_location, 'internal'));
                Misc.print_error_message(msg);
                
                success = 0;
                
            end
            
            % Check if hashes are the same
            if ~Setup.check_file_hash(result_file, ...
                    problem_setup.mesh_data_file_hash)
                success = 0;
                return
            end
            
            msg = sprintf('\tUsing mesh data file "%s"', fullfile('internal', ...
                problem_setup.mesh_data_file));
            Misc.print_message(msg);
            
            % If the previous run was terminated after the
            % meshing step, continue the computation procedure with
            % parsing
            if current_setup_state == problem_setup.state
                
                next_computation_step = Setup.setup_state_parsing_finished;
                loading_completed = true;
                Misc.print_message(['\nContinuing computation with assembling ', ...
                    'and solving procedure.']);
                
            else
                next_loading_step = Setup.setup_state_solving_finished;
            end
            
            
        end
        
        function [success, next_loading_step, loading_completed, next_computation_step] = ...
                loading_setup_check_solving_step( ...
                problem_setup, current_setup_state)
            
            
            success = 1;
            loading_completed = 0;
            next_loading_step = 0;
            next_computation_step = 0;
           
            problem_location = problem_setup.problem_location;
            
            
            result_file = ...
                fullfile('results', problem_setup.result_file);
            
            % Check if mesh data file exists
            if ~Misc.check_file_existence(problem_location, ...
                    result_file )
                success = 0;
                msg = sprintf(['result file "%s" in folder "%s" not ', ...
                    'found'], problem_setup.result_file, ...
                    fullfile(problem_location, 'results'));
                Misc.print_error_message(msg);
                return
            end
            
            % Check if hashes are the same
            if ~Setup.check_file_hash(result_file, ...
                    problem_setup.result_file_hash)
                success = 0;
                return
            end
            
            msg = sprintf('\tUsing result file "%s"', fullfile('results', ...
                problem_setup.result_file));
            Misc.print_message(msg);
            
            % If the previous run was terminated after the
            % meshing step, continue the computation procedure with
            % parsing
            if current_setup_state == problem_setup.state
                
                next_computation_step = Setup.setup_state_solving_finished;
                loading_completed = true;
                
                Misc.print_message('\nContinuing computation with postprocessing.');
                
                % Decide wheter to continue the loading procedure
            else
                next_loading_step = Setup.setup_state_postprocessing_finished;
            end
            
            
            
        end
        
        function [success, next_loading_step, loading_completed, next_computation_step] = ...
                loading_setup_check_postprocessing_step( ...
                problem_setup, current_setup_state)
            
            success = 1;
            loading_completed = 1;
            next_loading_step = 0;
            next_computation_step = Setup.setup_state_postprocessing_finished;
        end
        
        
    end
    
   
end

