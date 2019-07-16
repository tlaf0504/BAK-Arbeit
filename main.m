clear all
close all
clc

addpath(fullfile(pwd, 'entities'))

global ui;


% ============================= User area ===========================================
% Gmsh location
gmsh_location = '/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/gmsh-4.3.0-Linux64/bin';

% Starts user interface.
% When the UI is used, the parameters below have no effect.
ui = 0;


% Problem location
if strcmp(getenv('username'),'baumgartner') % workaround that both can execute the code properly
    problem_location = 'D:\LV\Studienarbeiten\Lafer\problems\cylinder_cap';
else
    
    problem_location = '/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/cylinder_cap';
    %problem_location = '/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/fem_test';
    %problem_location = '/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/simple_cap';
end

% Integer number representing the type of the problem.
% See Misc.supported_problem_types for supported problem types.
problem_type = 'Static Current';


% Set this flag to 1 if you want to load a existing setup. If this flag is set,
% the parameters below have no effect.
load_existing_setup = 1;

% Name of the geometry file
geometry_file = 'cylinder_cap.geo';
%geometry_file = 'fem_test.geo';
%geometry_file = 'cap.geo';

% Name of the settings file
settings_file = 'cylinder_cap.set';
%settings_file = 'fem_test.set';
%settings_file = 'cap.set';

% Element order.
% 1: Linear, 2: Quadratic, 3: Cubic
mesh_order = 2;

% Mesh file version
% See GmshIF.supported_mesh_file_versions for further information
mesh_file_version = 2;

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
    success = Setup.new_setup(problem_location, geometry_file, settings_file, ...
        mesh_order, problem_type, mesh_file_version);
    if ~success
        return
    end
    
    success = GmshIF.mesh(problem_location);
    if ~success
        return
    end
    
    success = Parser.parse(problem_location);
    if ~success
        return
    end
    
    success = Solver.solve(problem_location);
    if ~success
        return
    end
    
    Postprocessor.plot_electrostatic_potential(problem_location);
    
    
  
end

