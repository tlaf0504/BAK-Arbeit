% <one blank line> (Without comment key)
% coding_standard.m
%
% DESCRIPTION: Brief description of file content
%
%    Detailed description: Leading blank line and  intended by 4 spaces.
%    Multiline description also intended by 4 spaces
%
%    Line breaks should apply after 85 characters. (Default setup of Matlab)
%
%
%
%
% Author: <authors name>
% e-mail: jhon.doe@blablabla.com
%
% Created on <date of creation>
% Last modified <date of last modification> by <name and e-mail of editor>
%



classdef ExampleClass
    
% classdef ExampleClass
%
% DESCRIPTION: Brief description
%
%    Detailed description: Leading blank line and  intended by 4 spaces.
%    Multiline description also intended by 4 spaces.
%    At least two two blank lines after description.
%
%
%
% Author: <authors name>
% e-mail: jhon.doe@blablabla.com
%
% Created on <date of creation>
% Last modified <date of last modification> by <name and e-mail of editor>
%

    properties(Constant)
        property1 = 1; % Single line property description
        
        % Multiline property descriptions start one line above, leaving at least one
        % blank line to the previous propery.
        property2 = 2;
        
        % ===== Property group
        % Description of property group shall start with a space followed by five '='
        % characters and again followed by a space and then a brief description or 
        % symbolic name of the group.
        % 
        % For describing the single properties of the group the rules from above
        % apply.
        % A blank line should be inserted between this description and the properties
        % of the group.
        
        % Multiline description of property 3.
        property3 = 3;
        property4 = 4; % Single line description of property 4
    end
    
    methods
        
        function [outputArg1,outputArg2] = example_function(inputArg1,inputArg2)  
            
        % function [outputArg1,outputArg2] = example_function(inputArg1,inputArg2)
        %
        % DESCRIPTION: Brief description
        % 
        %    Detailed description: Leading blank line and  intended by 4 spaces.
        %    Multiline description also intended by 4 spaces.
        %    
        %
        %
        % INPUT:
        %     inputArg1.....Description
        %     inputArg2.....If an input or output parameter description goes over
        %       multiple lines, the lines are intended by one tab (4 spaces)
        %
        % OUTPUT:
        %     outputArg1.....Description
        %     outputArg2.....Description
        %
        % REQUIRED FILES / DEPENDENCIES / CREATED_FILES / UPDATED_FILES:
        %     examples.mat.....Description of external file dependencies. Can go
        %         over multiple lines.
        %
        %
        % Author: <authors name>
        % e-mail: jhon.doe@blablabla.com
        %
        % Created on <date of creation>
        % Last modified <date of last modification> by <name and e-mail of editor>
        %
            outputArg1 = inputArg1;
            outputArg2 = inputArg2;
        end
    end
    
end

