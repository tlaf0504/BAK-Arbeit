classdef Misc
    
    properties(Constant)
        supported_mesh_orders = [1,2,3];
        console_section_separator = repmat('=', 1, 80); 
        console_subsection_separator = repmat('-', 1, 80);
        
        % Supported problem types.
        %     -) 'Electrostatic': Electrostatic problem. ID = 1
        %     -) 'Static Current': Static currentflow problem. ID = 2
        supported_problem_types = containers.Map(...
            {'Electrostatic', 'Static Current', 'RotMagnetostatic', 'PlaneMagnetostatic'}, ...        
            [1, 2, 3, 4] ...
            );
        
    end
    
    methods(Static)
        function print_error_message(message)
            message = ['Error:\n\t', message, '\n\r'];
            fprintf([message, '\n\r']);
            
            global ui;
            if ui
                % Log area in GUI handles does not support simple
                % fprintf-functionalities. It rather handles each line as a element
                % of a cell array.
                % The following code converts each line of the message into an
                % element of a cell array.
                tmp = replace(message, '\r', '');
                tmp = replace(tmp, '\t', '');
                tmp = split(tmp, '\n');
                global msg_dest;
                msg_dest.Value = [msg_dest.Value; tmp];
                
                msgbox(tmp,'Error', 'error', 'modal')
                uiwait()
            end
        end
        
        
        function print_warning_message(message)
            message = ['Warning:\n\t', message, '\n\r'];
            fprintf([message, '\n\r']);
            global ui;
            if ui
                tmp = replace(message, '\r', '');
                tmp = replace(tmp, '\t', '');
                tmp = split(tmp, '\n');
                global msg_dest;
                msg_dest.Value = [msg_dest.Value; tmp];
                
                msgbox(tmp,'Warning','warn', 'modal')
                uiwait()
            end
        end
        
        function print_message(message)
            message = [message, '\n'];
            fprintf(message);
            global ui;
            if ui
                tmp = replace(message, '\r', '');
                tmp = replace(tmp, '\t', '');
                tmp = split(tmp, '\n');
                global msg_dest;
                msg_dest.Value = [msg_dest.Value; tmp];
            end
        end
        
        
        function file_exists = check_file_existence(location, filename)
            
            file_exists = 1;
            path = fullfile(location, filename);
            
            releasedate=datevec(version('-date'),'');
            if releasedate(1)>=2018
                if ~isfile(path)
                    file_exists = 0;
                end
            else
                if ~exist(path,'file')
                    file_exists = 0;
                end
            end
            
        end
        
        function [problem_type, success] = ...
                get_problem_type_class_from_problem_type_number( ...
                problem_type_number)
            
            success = 1;
            problem_type = nan;
            
            switch problem_type_number
                case 1
                    problem_type = ElectrostaticProblem;
                case 2
                    problem_type = StaticCurrentProblem;
                case 3
                    problem_type = RotationSymmetricMagnetostaticProblem;
                case 4
                    problem_type = PlaneMagnetostaticProblem;
                otherwise
                    msg = sprintf(['Wrong problem type input. Following problem ', ...
                        'types are supported: "Electrostatic", "StaticCurrent", "RotationSymmetricMagnetostaticProblem".']);
                    Misc.print_error_message(msg);
                    
                    success = 0;
                    return
            end
        end
        
        
        function finite_element_type = ...
                get_finite_element_class_from_mesh_order(mesh_order)

            if mesh_order == 1
                finite_element_type = FirstOrderTriangleElement;
                
            elseif mesh_order == 2
                finite_element_type = SecondOrderTriangleElement;
                
            elseif mesh_order == 3
                finite_element_type = ThirdOrderTriangleElement;
            end
        end
        
        function vacuum_material = get_vacuum_material(problem_type)
            if problem_type == 1 % Electrostatic problem
                % Vacuum permittivity
                vacuum_material = 8.854187812813e-12;
            elseif problem_type == 2 % Static current problem
                % Contuctivity is given as an absolute value for static current
                % problems. So 'vacuum_material' is simply set to 1.
            elseif problem_type == 3 || problem_type == 4
                % Plane or rotationsymmetric magnetostatic problem
                vacuum_material = 4 * pi * 10^-7;
            end
        end
        
        function CloneFig(inFigNum,OutFigNum)
            % This function was taken from
            % https://www.mathworks.com/matlabcentral/fileexchange/26587-clone-figure
            % (18.07.2019)
            % Big thanks to Matt Fetterman for his work.
            
            % this program copies a figure to another figure
            % example: CloneFig(1,4) would copy Fig. 1 to Fig. 4
            % Matt Fetterman, 2009
            % pretty much taken from Matlab Technical solutions:
            % http://www.mathworks.com/support/solutions/en/data/1-1UTBOL/?solution=1-1UTBOL
            hf1=figure(inFigNum);
            hf2=figure(OutFigNum);
            clf;
            Misc.compCopy(hf1,hf2);
            
        end
        
        function compCopy(op, np)
            % This function was taken from
            % https://www.mathworks.com/matlabcentral/fileexchange/26587-clone-figure
            % (18.07.2019)
            % Big thanks to Matt Fetterman for his work.
            
            
            %COMPCOPY copies a figure object represented by "op" and its %
            %descendants to another figure "np" preserving the same hierarchy.
            ch = get(op, 'children');
            if ~isempty(ch)
                nh = copyobj(ch,np);
                for k = 1:length(ch)
                    Misc.compCopy(ch(k),nh(k));
                end
            end
            return;
        end
        
    end
    
end