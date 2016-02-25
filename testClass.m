classdef testClass
    %testClass Model of a serial robot manipulator
    
    properties
        q
        T
    end
    
    properties (Dependent)
        
    end
    
    methods
        % Constructor 
        function obj = testClass(val)
            obj.q=sym('q',[1,val]);
            obj.T=cell(1,val);
            for i = 1 : val
                obj.T{i} = sym(strcat('T',int2str(i-1),'_',int2str(i)),[4,4]);
            end
        end
    end
    
end

