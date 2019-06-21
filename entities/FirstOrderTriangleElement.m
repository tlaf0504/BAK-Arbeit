classdef FirstOrderTriangleElement
    %FIRSTORDERTRIANGLEELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        nodes_of_triangle_edges = [1,1,0;
                                   1,0,1;
                                   0,1,1];
    end
    
    methods
        function obj = FirstOrderTriangleElement()
            %FIRSTORDERTRIANGLEELEMENT Construct an instance of this class
            %   Detailed explanation goes here
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

