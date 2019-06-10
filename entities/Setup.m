classdef Setup
    methods(Static)
        
        
        function [problem_setup, success] = setup(problem_location, geometry_filename)
            success = 1;
            cd(problem_location)
            
            [problem_setup, existing_solution] = ...
                Setup.check_for_existing_solution();
            
            if existing_solution
                return
            else
                if ~Misc.check_file_existence(problem_location, geometry_filename)
                    success = 0;
                    return
                end
                
                % Setting up folders for results, plots etc.
                Setup.setup_solution_hierachy();
                
                % Create new structure for storing problem information
                problem_setup = Setup.new_problem_setup(problem_location, ...
                    geometry_filename);
            end
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
                
                Misc.print_message(['No previous solution setups found.\nCreating new setup ', ...
                    'file "problem_setup.mat"'])
                
                existing_solution = 0;
                problem_setup = [];
                
            end
        end
        
        
        function problem_setup = create_problem_setup_structure()
            problem_setup = struct();
            
            problem_setup.problem_name = '';
            problem_setup.problem_path = '';
            problem_setup.geometry_file = '';
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
            problem_setup.problem_path = old_problem_setup.problem_path;
            problem_setup.geometry_file = old_problem_setup.geometry_file;
            problem_setup.geometry_modification_date = old_problem_setup. ...
                geometry_modification_date;
            
        end
        
        
        function problem_setup = new_problem_setup(location, ...
                geometry_filename)
            
            problem_setup = Setup.create_problem_setup_structure();
            
            [~, name, ~] = fileparts(geometry_filename);
            file_info = dir(geometry_filename);
            
            problem_setup.problem_name = name;
            problem_setup.problem_path = location;
            problem_setup.geometry_file = geometry_filename;
            problem_setup.mesh_file = [name, '.msh'];
            problem_setup.geometry_modification_date = file_info.datenum;   
        end
        
        
    end
end

