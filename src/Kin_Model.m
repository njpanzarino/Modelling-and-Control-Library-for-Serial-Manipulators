classdef Kin_Model
    %Kin_Model: Kinematic Model of a serial robot manipulator
    
    properties
    
    end
    
    properties(SetAccess = private)
        q
        n_frames
        J
        inv_J
        inv_J0
    end
    
    properties(Access = private)
        T_matrices
        lambda_val = 0.001  %small positive real value used to remove singularities in inverse Jacobian
        lambda_inv_J    %inverse Jacobian with symbolic lambda to quickly recalculate when lambda changes
        sym_lambda = sym('lambda');
    end
    
    properties (Dependent)
        lambda
    end
    
    methods

        function value = T(obj,from,to)
            value = obj.T_matrices{from+1,to+1};
        end 
        
        function value = forward_kin(obj,q)
            if nargin<2
                value = obj.T(0,obj.n_frames);
            else
                value = subs(obj.T(0,obj.n_frames),obj.q,q);
            end
        end
        
        function value = inv_kin(obj,H)
            error('function not implemented')
        end
        
        function obj = set.lambda(obj,value)
            obj.lambda_val = value;
            if ~isempty(obj.lambda_inv_J)
                obj.inv_J = simplify(vpa(subs(obj.lambda_inv_J,obj.sym_lambda,obj.lambda)));
            end
        end
        function value = get.lambda(obj)
            value = obj.lambda_val;
        end
    end
    
    methods(Static)
        function obj = fromH_TransChain(varargin)
            obj = Kin_Model;
            obj = obj.init(size(varargin{1},1));
            for i = 1:obj.n_Frames
                obj.T_matrices{i,i+1}=vpa(varargin{i});
            end
            obj=obj.calculateFromChain();
            obj=obj.calculateJacobian();
            obj=obj.calculatePesudoInvJacobian();
        end
        
        function obj=fromDH(DH_Matrix,jointVars)
            obj = Kin_Model;
            obj = obj.init(size(DH_Matrix,1));
            if nargin<2
                obj.q = sym('q',[obj.n_frames,1]);
            else
                obj.q = jointVars;
            end
            for i = 1:obj.n_frames
                obj.T_matrices{i,i+1}=simplify(vpa(H_Trans.fromDH(DH_Matrix(i,:)).H));
            end
            obj=obj.calculateFromChain();
            obj=obj.calculateJacobian();
            obj=obj.calculatePesudoInvJacobian();
        end
        
        function obj=fromPrompt()
            options.Interpreter = 'tex';
            options.Resize = 'on';
            dh_table = inputdlg('Enter the dh table for the robotic manipulator Use the following format: \theta, d, a, \alpha  \theta_1,d_1,a_1,\alpha_1; \theta_2,d_2,a_2,\alpha_2;                        ...','dh table',[3 50],{''},options);
            dh_data = cellstr(dh_table{1, 1});

            M = {};

            for i = 1:size(dh_data,1)
                M(i, 1:4) = strsplit(dh_data{i},',');
            end

            [h, w] = size(M);
            DH = sym('DH',[h w]);
            for r = 1:h
                for c = 1:w
                   DH(r,c) = sym(M{r,c});
                end
            end
            
            allvars = symvar(DH);
            
            arrayString=char(allvars);
            arrayString([1:8,end-1:end]) = [];
            cells=strsplit(arrayString(2:end-1),', ').';

            [S,~] = listdlg('ListString',cells,'PromptString','Select Joint Variables');

            if isempty(S)
                jointVars=allvars;
            else
                jointVars=allvars(S);
            end
            
            obj=Kin_Model.fromDH(DH,jointVars.');
        end
    end
    
    methods(Access = private)
        function obj = init(obj,frames)
            obj.n_frames=frames;
            obj.T_matrices=cell(frames+1,frames+1);
            for i = 1:frames+1
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
        
        function obj = calculateJacobian(obj)
            obj.J=H_Trans(obj.T(0,size(obj.q))).getJacobian(obj.q);
        end
        
        function obj = calculatePesudoInvJacobian(obj,lambda)
            if nargin<2
                lambda = obj.lambda;
            end
            
            obj.lambda_inv_J = simplify(vpa(obj.J.'/(obj.J*obj.J.' + (obj.sym_lambda^2).*eye(max(size(obj.J))))));
            obj.inv_J0 = simplify(vpa(subs(obj.lambda_inv_J,obj.sym_lambda,0)));
            obj.lambda = lambda;
        end
    end
    
end

