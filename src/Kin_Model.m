classdef Kin_Model
    %Kin_Model Model of a serial robot manipulator
    
    properties
        q
        DOF
        J
        inv_J
    end
    
    properties(Access = public)
        T_matrices
    end
    
    properties (Dependent)
        
    end
    
    methods
        % Constructor. Takes in Chain of Homogeneous Transformation
        % Matrices
        
        function value = T(obj,from,to)
            value = obj.T_matrices{from+1,to+1};
        end 
    end
    
    methods(Static)
        function obj = fromH_TransChain(varargin)
            obj = Kin_Model;
            obj = obj.init(size(DH_Matrix,1));
            for i = 1:obj.DOF
                obj.T_matrices{i,i+1}=vpa(varargin{i});
            end
            obj=obj.calculateFromChain();
        end
        
        function obj=fromDH(DH_Matrix)
            obj = Kin_Model;
            obj = obj.init(size(DH_Matrix,1));
            for i = 1:obj.DOF
                obj.T_matrices{i,i+1}=simplify(vpa(H_Trans.fromDH(DH_Matrix(i,:)).H));
            end
            obj=obj.calculateFromChain();
        end
        
    end
    
    methods(Access = public)
        function obj = init(obj,size)
            obj.DOF=size;
            obj.q=sym('q',[1,obj.DOF]);
            obj.T_matrices=cell(size+1,size+1);
            for i = 1:size+1
                obj.T_matrices{i,i} = eye(4);
            end
        end
        
        function obj = calculateFromChain(obj)
            for r = 1:size(obj.T_matrices,1)
               for c = 1:size(obj.T_matrices,2) 
                   if c>(r+1)
                       temp=eye(4);
                       for i=r:(c-1)
                           temp=temp*obj.T_matrices{i,i+1};
                       end
                       obj.T_matrices{r,c} = simplify(vpa(temp));
                   end
               end
            end
            
            for r = 1:size(obj.T_matrices,1)
               for c = 1:size(obj.T_matrices,2) 
                   if c<r
                       obj.T_matrices{r,c} = simplify(vpa(inv(obj.T_matrices{c,r})));
                   end
               end
            end
        end
        
        function obj = calculateFromAbsolute(obj)
            for i=2:size(obj.T_matrices,2)
                obj.T_matrices{i,1} = simplify(vpa(inv(obj.T_matrices{1,i})));
            end
            
            for r = 2:size(obj.T_matrices,1)
               for c = 2:size(obj.T_matrices,2)
                   if c~=r
                       obj.T_matrices{r,c} = simplify(vpa(obj.T_matrices{r,1}*obj.T_matrices{1,c}));
                   end
               end
            end
        end
        
        function obj = calculateJacobian(obj)
            obj.J=H_Trans(obj.T(0,obj.DOF)).getJacobian(obj.q);
        end
    end
    
end

