classdef Kin_Model
    %Kin_Model: Kinematic Model of a serial robot manipulator
    
    properties
		q_limit
    end
    
    properties(SetAccess = private)
        q
        
        n_frames
        
    end
    
    properties(Access = private)
        T_matrices
        
        sym_lambda = sym('lambda');

        f_kin_func
        J_func
        inv_JLambda_func
        inv_J0_func
    end
    
    properties (Dependent)
        
    end
    
    methods

        function value = T(obj,from,to)
            value = obj.T_matrices{from+1,to+1};
        end 
        
        function value = forward_kin(obj,q,useEuler)
            if nargin<2
                q=obj.q;
            end
            if nargin<3
                useEuler=false;
            end
            
            if useEuler
                T=obj.f_kin_func(q);
                value=[T.Trans;T.Euler];
            else
                value=obj.f_kin_func(q);
            end
        end
        
        function value = inverse_kin(obj,H,q0)
            if nargin<2
                H=sym('H',4);
            end
            if nargin<3
                q0=zeros(size(obj.q));
            end
            
            options = optimoptions(@fsolve,'Display','iter',...
                'Algorithm','trust-region-reflective',...
                'Jacobian','off');
            %baisc fsolve
            [value,~]=fsolve(@(q)(obj.forward_kin(q)-H),q0,options);
            
            %Inverse Jacobian Method (Hill Climb)
            
            %Cyclic Coordinate Descent
			%Simulated Annealing
			%Map areas separated by singularities
            
        end
        
        function value = J(obj,q)
            if nargin>1
                value=obj.J_func(q);
            else
                value=obj.J_func(obj.q);
            end
        end
        
        function value = inv_J(obj,q,lambda)
            switch nargin
                case 1
                    value=obj.inv_J0_func(obj.q);
                case 2
                    value=obj.inv_J0_func(q);
                case 3
                    if lamda<=0
                        value=obj.inv_J0_func(q);
                    else
                        value=obj.inv_JLambda_func(q,lambda);
                    end
            end
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
            jointVars=reshape(jointVars,numel(jointVars),1);
            if nargin<2
                obj.q = sym('q',[obj.n_frames,1]);
            else
                obj.q = jointVars;
            end
            assume(obj.q, 'real');
            obj.q_limit=cell([size(obj.q,1),2]);
            for i = 1:obj.n_frames
                obj.T_matrices{i,i+1}=H_Trans.fromDH(DH_Matrix(i,:));
            end
            obj=obj.calculateFromChain();
            obj=obj.calculateForwardKinematics();
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
    
%     methods(Static,Access=private)
%         function f = funcFromSym(expr,vars)
%             if size(setdiff(symvar(expr),vars))>0
%                 expr=simplify(vpa(expr));
%                 f = @(q)simplify(vpa(subs(expr,vars,q)));
%             else
%                 f = matlabFunction(expr,'Vars',{vars});
%             end
%         end
%     end
    
    methods(Access = private)
        function obj = init(obj,frames)
            obj.n_frames=frames;
            obj.T_matrices=cell(frames+1,frames+1);
            for i = 1:frames+1
                obj.T_matrices{i,i} = H_Trans();
            end
        end
        
        function obj = calculateFromChain(obj)
            for r = 1:size(obj.T_matrices,1)
               for c = 1:size(obj.T_matrices,2) 
                   if c>(r+1)
                       temp=eye(4);
                       for i=r:(c-1)
                           temp=temp*obj.T_matrices{i,i+1}.H;
                       end
                       obj.T_matrices{r,c} = H_Trans(vpa(temp));
                   end
               end
            end
            
            for r = 1:size(obj.T_matrices,1)
               for c = 1:size(obj.T_matrices,2) 
                   if c<r
                       obj.T_matrices{r,c} = obj.T_matrices{c,r}.inv;
                   end
               end
            end
        end
        
        function obj = calculateForwardKinematics(obj)
            T=obj.T(0,obj.n_frames);
            obj.f_kin_func = T.getFunction(obj.q);
        end
        
        function obj = calculateJacobian(obj)
            J=obj.T(0,obj.n_frames).getJacobian(obj.q);
            obj.J_func=H_Trans.createFunction(J,obj.q);
        end
        
        function obj = calculatePesudoInvJacobian(obj)
            J=obj.J;

            lambda_inv_J = J.'/(J*J.' + (obj.sym_lambda^2).*eye(size(J,1)));
            
            obj.inv_JLambda_func = H_Trans.createFunction(lambda_inv_J,{obj.q,sym('lambda')});
%             if size(setdiff(symvar(lambda_inv_J),[obj.q;obj.sym_lambda]))>0
%                 lambda_inv_J = simplify(vpa(lambda_inv_J));
%                 obj.inv_JLambda_func = @(q,lambda)simplify(vpa(subs(lambda_inv_J,[obj.q;obj.sym_lambda],[q;lambda])));
%             else
%                 obj.inv_JLambda_func = matlabFunction(lambda_inv_J,'Vars',{obj.q,obj.sym_lambda});
%             end
                
            inv_J0 = vpa(subs(lambda_inv_J,obj.sym_lambda,0));
            obj.inv_J0_func=H_Trans.createFunction(inv_J0,obj.q);
        end
        
%         function obj = train_anfis(obj,n)
%             disp('start training')
%             obj.anfis_net=cell(size(obj.q));
%             default_space=linspace(-pi,pi,n);
%             
%             input=cell(1,size(obj.q,1));
%             
%             for i=1:size(obj.q,1)
%                 if isempty(obj.q_limit{i})
%                     input{i}=default_space;
%                 else
%                     input{i}=linspace(obj.q_limit{i}(1),obj.q_limit{i}(2),n);
%                 end
%             end
%             G=cell(1,size(obj.q,1));
%             [G{:}]=ndgrid(input{:});
%             
%             all_H=obj.forward_kin(G.');
%             s=size(all_H);
%             r_sizes=ones(s(1)/4,1).*4;
%             c_sizes=ones(s(2)/4,1).*4;
%             all_H=mat2cell(all_H,r_sizes,c_sizes,ones(s(3),1));
% %             num=numel(G{1});
% %             for i=1:numel(G)
% %                 G{i}=reshape(G{i},[num,1]);
% %             end
% %             
% %             all_H=cell(numel(G{i}),1);
% %             for i=1:num
% %                 j=zeros(size(obj.q));
% %                 for n=1:numel(j)
% %                     j(n)=G{n}(i,1);
% %                 end
% %                 all_H{i}=obj.forward_kin(j);
% %             end
%         end
    end
    
end

