classdef Misc
    
    properties(Constant)
        supported_mesh_orders = [1,2,3];
        console_section_separator = repmat('=', 1, 80); 
        console_subsection_separator = repmat('-', 1, 80);
        
        % Supported problem types.
        %     -) 'Electrostatic': Electrostatic problem. ID = 1
        %     -) 'Static Current': Static currentflow problem. ID = 2
        supported_problem_types = containers.Map(...
            {'Electrostatic', 'Static Current'}, ...        
            [1, 2] ...
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
    end
    
end