clear all
close all
clc

addpath(fullfile(pwd, 'entities'))
addpath(fullfile(pwd, 'logo'))

global ui;


% ============================= User area ===========================================
% Gmsh location
gmsh_location = '/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/gmsh-4.3.0-Linux64/bin';

% Starts user interface.
% When the UI is used, the parameters below have no effect.
ui = 1;


% Problem location
if strcmp(getenv('username'),'baumgartner') % workaround that both can execute the code properly
    problem_location = 'D:\LV\Studienarbeiten\Lafer\problems\cylinder_cap';
else
    problem_location = '/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/plate_capacitor_multi_material_2';
    %problem_location = '/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/plate_cap_multi_material';
    %problem_location = '/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/cylinder_cap';
    %problem_location = '/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/fem_test';
    %problem_location = '/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/simple_cap';
    %problem_location = '/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/cap_multi_material';
end

% Integer number representing the type of the problem.
% See Misc.supported_problem_types for supported problem types.
problem_type = 'Electrostatic';


% Set this flag to 1 if you want to load a existing setup. If this flag is set,
% the parameters below have no effect.
load_existing_setup = 1;

% Name of the geometry file
%geometry_file = 'cylinder_cap_simple.geo';
%geometry_file = 'cylinder_cap_2_materials.geo';
%geometry_file = 'fem_test.geo';
geometry_file = 'cap.geo';
%geometry_file = 'cap_multi_material.geo';

% Name of the settings file
%settings_file = 'cylinder_cap_2_materials.set';
%settings_file = 'cylinder_cap_simple.set';
%settings_file = 'fem_test.set';
settings_file = 'cap.set';
%ettings_file = 'cap_multi_material.set';

% Element order.
% 1: Linear, 2: Quadratic, 3: Cubic
mesh_order = 2 ;

% Mesh file version
% See GmshIF.supported_mesh_file_versions for further information
mesh_file_version = 2;

% Load pregenerated mesh file
load_mesh = false;
mesh_file = 'cap.msh';

% ========================= Do not mofidy code blow here ============================

% Matlab uses double by default.
mesh_order = int32(mesh_order);
mesh_file_verion = int32(mesh_file_version);

global gmsh_exec;

if isunix
    gmsh_exec = fullfile(gmsh_location, 'gmsh');
elseif ispc
    gmsh_exec = fullfile(gmsh_location, 'gmsh.exe');
elseif ismac
    error("I don't like apples.")
end


if ui
    main_window()
else
    computation_step = 0;


computation_completed = false;

while(~computation_completed)
    
    
    switch computation_step
        
        % Step for setting up a new problem
        case Setup.setup_state_initializing
            
            if load_existing_setup
                [success, computation_step] = Setup.load_setup(problem_location);
                if ~success
                    return
                end
            else
                success = Setup.new_setup(problem_location, geometry_file, settings_file, ...
                    mesh_order, problem_type, mesh_file_version, load_mesh, mesh_file);
                if ~success
                    return
                end
                
                % Skip meshing if a pregenerated mesh should be used
                if load_mesh
                    computation_step = Setup.setup_state_meshing_finished;
                else
                    computation_step = Setup.setup_state_initializing_finished;
                end
            end
            
            
        
        % Mesh generation step
        case Setup.setup_state_initializing_finished
            
            success = GmshIF.mesh(problem_location);
            if ~success
                return
            end
            computation_step = Setup.setup_state_meshing_finished;
        
        % Parsing step    
        case Setup.setup_state_meshing_finished
            success = Parser.parse(problem_location);
            if ~success
                return
            end
            computation_step = Setup.setup_state_parsing_finished;
            
       % Assembling and solving step
        case Setup.setup_state_parsing_finished
            success = Solver.solve(problem_location);
            if ~success
                return
            end
            computation_step = Setup.setup_state_solving_finished;
            
        % Postprocessing step    
        case Setup.setup_state_solving_finished
            success = Postprocessor.run(problem_location);
            if ~success
                return
            end
            
            computation_completed = true;
    end
end
end

