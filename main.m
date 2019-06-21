clear all
close all
clc

addpath(fullfile(pwd, 'entities'))

global ui;


% ============================= User area ===========================================
% Gmsh location
% Please enter the path to the MAIN FOLDER of the Gmsh SDK package. You can download 
% it from http://gmsh.info/#Download.
%
% Extract the archive to a folder of your choice and enter the path of the main
% folder in the parameter below.
%
% e.g. gmsh_location = '/home/<user>/gmsh/gmsh-4.3.0-Linux64'
%
%
% Please ensure you download the SDK, as this software uses the powerful Gmsh 
% built-in API, which is only delivered with the SDK version.
gmsh_location = '/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/gmsh-4.3.0-Linux64';

% Starts user interface.
% When the UI is used, the parameters below have no effect.
ui = 0;


% Problem location
problem_location = '/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/cylinder_cap';

% Set this flag to 1 if you want to load a existing setup. If this flag is set,
% the parameters below have no effect.
load_existing_setup = 1;

% Name of the geometry file
geometry_file = 'cylinder_cap.geo';

% Name of the settings file
settings_file = 'cylinder_cap.set';

% Element order. 
% 1: Linear, 2: Quadratic, 3: cubic
mesh_order = 1;

% ========================= Do not mofidy code blow here ============================
if ~GmshAPI.check_installation(gmsh_location)
    return
end

addpath(fullfile(gmsh_location, 'include'))
addpath(fullfile(gmsh_location, 'lib'))
if not(libisloaded('libgmsh'))
    loadlibrary('libgmsh.so', 'gmshc.h')
end

if ui
    main_window()
else
    success = Setup.new_setup(problem_location, geometry_file, settings_file, mesh_order);
    if ~success
        return
    end
    
    GmshAPI.mesh(problem_location);
    Parser.parse(problem_location);
    
    
end

% tic
% Parser.parse('/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/cylinder_cap', 'cylinder_cap.geo', 'cylinder_cap.set')
% toc


        