classdef GmshAPI
    
    properties(Constant)
        
        supported_mesh_file_versions = [1,2,4];
        
        % Exporting a specific mesh file format is depended of the extension while
        % exporting a file. 
        % See http://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options,
        % command line option "-format"
        %
        % Supported mesh file versions of this software: Version 2.2
        % (More support maybe added in future)
        %
        % For mesh file verion 1, the extension '.msh1' is used.
        % For mesh file version 2(.2), the extension '.msh22' is used
        % For mesh file version 4, the estension '.msh41' is used. 
        %
        mesh_file_extension_map = containers.Map( ...
            GmshAPI.supported_mesh_file_versions, ... % keys
            {'.msh1', '.msh22', '.msh41'} ... % values
            );
    end
    
    methods(Static)
        function mesh(problem_location)
            tmp = pwd;
            cd(problem_location);
            
            load('problem_setup.mat', 'problem_setup');
            
            if problem_setup.state < Setup.setup_state_initializing_finished
                msg = sprintf(['Error. Problem initilization were not ', ...
                    'finished. Maybe the setup procedure was interrupted.', ...
                    'Please rerun the complete setup procedure.']);
                Misc.print_error_message(msg) 
            end
            
            % Check if mesh order is valid
            if ~ismember(problem_setup.mesh_order, Misc.possible_mesh_orders)
                msg = sprintf(['Error. Specified mesh_order %d is not one ', ...
                    'of the following allowed mesh orders:\n%s'], ...
                    problem_setup.mesh_order, ...
                    sprintf('%s ', Misc.possible_mesh_orders));
                Misc.print_error_message(msg);
            end
            
            
            % Initialize Gmsh API
            if ~GmshAPI.init()
                Misc.print_error_message(['Error while initializing Gmsh API. ', ...
                    'Check Gmsh installaiton.'])
                return
            end

            GmshAPI.load_geometry(problem_setup);
            GmshAPI.generate_2D_mesh(problem_setup.mesh_order);
            problem_setup = GmshAPI.export_mesh_file(problem_setup, 2);
            GmshAPI.finalize();
            
            problem_setup.state = Setup.setup_state_meshing_finished;
            Setup.update_problem_setup_file(problem_setup);     
            cd(tmp)
            
            
            
        end
        
        function load_geometry(problem_setup)
            file = fullfile(problem_setup.problem_location, ...
                problem_setup.geometry_file);
            ierr = int32(0);
            
            calllib('libgmsh', 'gmshOpen', file, ierr);
        end
        
        function success = check_installation(gmsh_location)
            
            success = 1;
            % Check gmsh installation
            gmsh_header_location = fullfile(gmsh_location, 'include');
            if ~Misc.check_file_existence(gmsh_header_location, 'gmshc.h')
                msg = sprintf(['Error. Expected Gmsh API header file ', ...
                    '"gmshc.h" in directory "%s" not found. Check path ', ...
                    'and installation of gmsh.'], gmsh_header_location);
                Misc.print_error_message(msg);
                success = 0;
            end
            
            gmsh_library_location = fullfile(gmsh_location, 'lib');
            if ~Misc.check_file_existence(gmsh_library_location, 'libgmsh.so')
                msg = sprintf(['Error. Expected Gmsh library file ', ...
                    '"libgmsh.so" in directory "%s" not found. Check path ', ...
                    'and installation of gmsh.'], gmsh_library_location);
                Misc.print_error_message(msg);
                success = 0;
            end
        end
        
        function success = init()
            ret = int32(0);
            success = 1;
            
            %Initialize Gmsh API
            calllib('libgmsh', 'gmshInitialize', 0, {}, 1, ret);
            if ret ~= 0
                success = 0;
                return
            end
        end
        
        function finalize()
            ierr = int32(0);
            calllib('libgmsh', 'gmshFinalize', ierr)
        end
        
        function set_mesh_order(order)
            ierr = int32(0);
            order = int32(order);
            calllib('libgmsh', 'gmshModelMeshSetOrder', order, ierr);
        end
        
        function generate_mesh(dimension)
            ierr = int32(0);
            dimension = int32(dimension);
            calllib('libgmsh', 'gmshModelMeshGenerate', dimension, ierr);
        end
        
        function generate_2D_mesh(order)
            GmshAPI.set_mesh_order(order);
            GmshAPI.generate_mesh(2); 
        end
        
        function refine_mesh()
            calllib('libgmsh', 'refine')
        end
        
        function problem_setup = export_mesh_file(problem_setup, version)
            
            if ~ismember(version, GmshAPI.supported_mesh_file_versions)
                msg = sprintf(['Error. Mesh file version %d not supported. ', ...
                    'Supported versins are: %s'], version, ...
                    sprintf('%d ', GmshAPI.supported_mesh_file_versions));
                Misc.print_error_message(msg);
            end
            
            problem_location = problem_setup.problem_location;
            mesh_file_version_extension = GmshAPI.mesh_file_extension_map(version);
            
            mesh_file_name = [problem_setup.problem_name, mesh_file_version_extension];
            problem_setup.mesh_file = mesh_file_name;
            
            GmshAPI.write_file(fullfile(problem_location, mesh_file_name));
        end
        
        function write_file(name)
            ierr = int32(0);
            calllib('libgmsh', 'gmshWrite', name, ierr);
        end
    end
end