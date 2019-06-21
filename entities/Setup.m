classdef Setup
    
    properties(Constant)
        setup_state_initializing = 0;
        setup_state_initializing_finished = 1;
        setup_state_meshing_finished = 2;
        setup_state_parsing_finished = 3;
        setup_state_solving_finished = 4;
        setup_state_postprocessing_finished = 5;
    end
    methods(Static)
        
        
        function success = new_setup(problem_location, ...
                geometry_filename, settings_filename, mesh_order)
            success = 1;
            
            tmp = pwd;
            cd(problem_location)
            
            if  ~Setup.check_geometry_and_settings_file_existence(problem_location, ...
                    geometry_filename, settings_filename)
                success = 0;
                return
            end

            % Setting up folders for results, plots etc.
            Setup.setup_solution_hierachy();
            
            % Create new structure for storing problem information
            problem_setup = Setup.new_problem_setup(problem_location, ...
                geometry_filename, settings_filename, mesh_order);
            
            problem_setup.state = Setup.setup_state_initializing_finished;
            
            save('problem_setup.mat', 'problem_setup')
            
            cd(tmp)
            
        end
        
        function [problem_setup, existing_solution] = check_for_existing_solution()
            
            % If no setup file is found, setup the folder hierachy for the solution.
            % Old folders and their contents are deleted in case of a new problem.
            
            
            if Misc.check_file_existence(pwd, 'problem_setup.mat')
                
                Misc.print_message(['Problem setup file "problem_setup.mat" found. \nLoading ', ...
                    'configuarion.'])
                
                S = load('problem_setup.mat', 'problem_setup');
                problem_setup = S.problem_setup;
                clear S;
                
                if check_for_newer_problem_version(problem_setup)
                    % Clean up existing direcories and setup structure
                    Misc.setup_solution_hierachy();
                    problem_setup = Setup.create_problem_setup_structure();
                    
                    existing_solution = 0;
                else
                    existing_solution = 1;
                end
                
                
            else
                
                Misc.print_message(['No previous solution setups found.'])
                
                existing_solution = 0;
                problem_setup = [];
                
            end
        end
        
        function success = check_geometry_and_settings_file_existence(problem_location, ...
                geometry_filename, settings_filename)
            
            success = 1;
            % Check for specifed geometry file in problem folder
            if ~Misc.check_file_existence(problem_location, geometry_filename)
                success = 0;
                msg = sprintf(['Error. Geometry file "%s" in folder "%s" not ', ...
                    'found'], problem_location, geometry_filename);
                Misc.print_error_message(msg);
                success = 0;
                return
            end
            
            % Check for specifed geometry file in problem folder
            if ~Misc.check_file_existence(problem_location, settings_filename)
                msg = sprintf(['Error. Settings file "%s" in folder "%s" not ', ...
                    'found'], problem_location, geometry_filename);
                Misc.print_error_message(msg);
                success = 0;
                return
            end
            
        end
        
        
        function problem_setup = create_problem_setup_structure()
            problem_setup = struct();
            
            problem_setup.problem_name = '';
            problem_setup.problem_location = '';
            problem_setup.geometry_file = '';
            problem_setup.settings_file = '';
            problem_setup.mesh_file = '';
            problem_setup.problem_type = '';
            problem_setup.geometry_modification_date = '';
        end
        
        
        function setup_solution_hierachy()
            
            foldername = 'internals';
            Setup.check_and_cleanup_folder(foldername);
            
            foldername = 'results';
            Setup.check_and_cleanup_folder(foldername);
            
            foldername = 'plots';
            Setup.check_and_cleanup_folder(foldername);
        end
        
        
        function check_and_cleanup_folder(folder)
            if exist(folder, 'dir')
                rmdir(folder, 's');
            end
            
            mkdir(folder)
        end
        
        
        function problem_modified = check_for_newer_problem_version(problem_setup)
            
            file_info = dir([problem_setup.location, problem_setup.geometry_file]);
            
            if (problem_setup.geometry_modification_date < file_info.datenum)
                Misc.print_message(['Geometry file was modified since last execution. ', ...
                    'Cleaning up hierachies.'])
                
                problem_modified = 1;
            else
                problem_modified = 0;
            end
        end
        
        
        function problem_setup = clean_problem_setup_structure(old_problem_setup)
            problem_setup = Setup.create_problem_setup_structure();
            
            problem_setup.problem_name = old_problem_setup.problem_name;
            problem_setup.problem_location = old_problem_setup.problem_location;
            problem_setup.geometry_file = old_problem_setup.geometry_file;
            problem_setup.geometry_modification_date = old_problem_setup. ...
                geometry_modification_date;
            
        end
        
        
        function problem_setup = new_problem_setup(location, ...
                geometry_filename, settings_filename, mesh_order)
            
            problem_setup = Setup.create_problem_setup_structure();
            
            [~, name, ~] = fileparts(geometry_filename);
            file_info_geo = dir(fullfile(location, geometry_filename));
            file_info_set = dir(fullfile(location, settings_filename));
            
            
            problem_setup.problem_name = name;
            problem_setup.problem_location = location;
            problem_setup.geometry_file = geometry_filename;
            problem_setup.settings_file = settings_filename;
            problem_setup.geometry_file_modification_date = file_info_geo.datenum;
            problem_setup.settings_file_modification_data = file_info_set.datenum;
            problem_setup.mesh_order = mesh_order;
            problem_setup.state = Setup.setup_state_initializing;
            
            % Placeholders
            problem_setup.mesh_file = '';
            problem_setup.mesh_data_file = '';
            problem_setup.result_data_file = '';
            
        end
        
        function update_problem_setup_file(problem_setup)
            save(fullfile(problem_setup.problem_location, 'problem_setup.mat'), ...
                'problem_setup');
        end
        
        
    end
end

