classdef Kin_Model < handle
    %Kin_Model: Kinematic Model of a serial robot manipulator
    
    properties
        q_limit
    end
    
    properties(SetAccess = private)
        n_frames
    end
    
    properties(Dependent, SetAccess=private)
        q    %joint variables
        d_q  %1st time derivative of joint variables
        dd_q %2nd time derivative of joint variables
    end
    
    properties(Access = private)
        val_q    %joint variables
        val_d_q  %1st time derivative of joint variables
        val_dd_q %2nd time derivative of joint variables
        
        T_matrices
        
        sym_lambda = sym('lambda','real');
        
        f_kin_func
        
        J_func
        inv_JLambda_func
        inv_J0_func
        d_J_func
        
        Ja_func
        inv_JaLambda_func
        inv_Ja0_func
        d_Ja_func
        
        trans_func %Function for determining points/transformations for each link when plotting
    end
    
    properties (Dependent)
        
    end
    
    methods
        
        function value = T(obj,from,to)
            value = obj.T_matrices{from+1,to+1};
            if isempty(value)
                obj=obj.calculateTransformation(from,to);
                value = obj.T_matrices{from+1,to+1};
            end
        end
        
        function value = forward_kin(obj,q,asWrench)
            if nargin<2
                q=obj.q;
            end
            if nargin<3
                asWrench=false;
            end
            
            if isempty(obj.f_kin_func)
                obj=obj.calculateForwardKinematics();
            end
            
            if asWrench
                T=H_Trans(obj.f_kin_func(q));
                value=T.Wrench;
            else
                value=obj.f_kin_func(q);
            end
        end
        function value = inverse_kin(obj,d_pose,q0)
            sz=size(d_pose);
            if isequal(sz,[4,4])
                asWrench=false;
            elseif isequal(sz,[6,1])
                asWrench=true;
            else
                error('POSE must be 4x4 homogeneous transformation matrix or 6x1 wrench');
            end
            
            if nargin<3
                q0=zeros(size(obj.q));
            end
            
            value=basicFSolve;
            
            %baisc fsolve
            function value = basicFSolve()
                options = optimoptions(@fsolve,'Display','none',...
                    'Algorithm','trust-region-reflective',...
                    'TolFun',1e-8,'Jacobian','off');
                [value,fval,flag]=fsolve(@(q)(obj.forward_kin(q,asWrench)-d_pose),q0,options);
                if flag>1
                    [value,fval,flag]=fsolve(@(q)(obj.forward_kin(q,asWrench)-d_pose),value,options);
                end
            end
            
            function value = wrenchFSolve()
                if ~asWrench
                    d_pose=H_Trans(d_pose).Wrench;
                    asWrench=true;
                end
                options = optimoptions(@fsolve,'Display','iter',...
                    'Algorithm','trust-region-reflective',...
                    'TolFun',1e-8,'Jacobian','off');
                [value,fval,flag]=fsolve(@(q)(obj.forward_kin(q,asWrench)-d_pose),q0,options);
                if flag>1
                    [value,fval,flag]=fsolve(@(q)(obj.forward_kin(q,asWrench)-d_pose),value,options);
                end
            end
            %Jacobian fsolve
%             W=H_Trans(H).Wrench;
%             options = optimoptions(@fsolve,'Display','iter',...
%                 'Algorithm','trust-region-reflective',...
%                 'Jacobian','on');
%             [value,~]=fsolve(@solve_fun_J,q0,options);
%             
%             function [F,J] = solve_fun_J(q)
%                 F=H_Trans(obj.forward_kin(q)).Wrench-W;
%                 J=obj.J(q);
%             end
            %Inverse Jacobian Method (Hill Climb)
            
            %Cyclic Coordinate Descent
            %Simulated Annealing
            %Map areas separated by singularities
            
        end
        
        function value = forward_vel(obj,q,d_q)
            value=obj.J(q)*d_q;
        end
        function value = inverse_vel(obj,q,wrench)
            value=obj.inv_J(q)*wrench;
        end
        
        function value = forward_acc(obj,q,d_q,dd_q)
            value=obj.J(q)*dd_q+obj.d_J(q,d_q)*d_q;
        end
        function value = inverse_acc(obj,q,d_q,wrench)
            value=obj.inv_J(q)*(wrench-obj.d_J(q,d_q)*d_q);
        end
        
        function value = J(obj,q)
            if isempty(obj.J_func)
                obj.calculateJacobian();
            end
            
            if nargin>1
                value=obj.J_func(q);
            else
                value=obj.J_func(obj.q);
            end
        end
        function value = inv_J(obj,q,lambda)
            if isempty(obj.inv_J0_func)
                obj.calculatePesudoInvJacobian();
            end
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
        function value = d_J(obj,q,d_q)
            if isempty(obj.d_J_func)
                obj.calculateJacobianDot();
            end
            
            if nargin>2
                value=obj.d_J_func(q,d_q);
            elseif nargin>1
                value=obj.d_J_func(q,obj.d_q);
            else
                value=obj.d_J_func(obj.q);
            end
        end
        
        function value = Ja(obj,q)
            if isempty(obj.Ja_func)
                obj.calculateAnalyticJacobian();
            end
            if nargin>1
                value=obj.Ja_func(q);
            else
                value=obj.Ja_func(obj.q);
            end
        end
        function value = inv_Ja(obj,q,lambda)
            if isempty(obj.inv_Ja0_func)
                obj.calculatePesudoInvAnalyticJacobian();
            end
            switch nargin
                case 1
                    value=obj.inv_Ja0_func(obj.q);
                case 2
                    value=obj.inv_Ja0_func(q);
                case 3
                    if lamda<=0
                        value=obj.inv_Ja0_func(q);
                    else
                        value=obj.inv_JaLambda_func(q,lambda);
                    end
            end
        end
        function value = d_Ja(obj,q,d_q)
            if isempty(obj.d_Ja_func)
                obj.calculateAnalyticJacobianDot();
            end
            
            if nargin>2
                value=obj.d_Ja_func(q,d_q);
            elseif nargin>1
                value=obj.d_Ja_func(q,obj.d_q);
            else
                value=obj.d_Ja_func(obj.q);
            end
        end
        
        function draw(obj,q,drawFrames,ax,plotArgs,viewArgs)
            
            if nargin<2 || isempty(q)
                q=zeros(size(obj.q));
            end
            
            if nargin<3 || isempty(drawFrames)
                drawFrames=0;
            end
            
            if nargin<4 || isempty(ax)
                ax=gca;
            end
            
            oldX=ax.XLim;
            oldY=ax.YLim;
            oldZ=ax.ZLim;
            
            if nargin<5 || isempty(plotArgs)
                plotArgs={};
            end
            
            if isempty(obj.trans_func)
                obj=obj.calculatePlot();
            end
            
            h=obj.trans_func(q);
            p=reshape(h(1:3,4,:),3,[]);
            scale=norm(p);
            
            plot3(ax,p(1,:),p(2,:),p(3,:),plotArgs{:});
            
            scale=scale/8;
            
            hchek = ishold;
            hold on
            
            H_Trans(h(:,:,1)).draw(scale,'0',ax);
            
            if drawFrames>0
                for i=1:obj.n_frames
                    if drawFrames>1
                        H_Trans(h(:,:,i+1)).draw(scale,num2str(i),ax);
                    else
                        H_Trans(h(:,:,i+1)).draw(scale,[],ax);
                    end
                end
            else
                h=plot3(ax,p(1,:),p(2,:),p(3,:),'.k',plotArgs{:});
                h.MarkerSize=h.MarkerSize+h.LineWidth*2;
            end
            
            if hchek == 0
                hold off
            end
            
            % Normalize Axes
            lX=ax.XLim;
            lY=ax.YLim;
            lZ=ax.ZLim;
            
            lX(1)=min([oldX(1),lX(1)]);
            lX(2)=max([oldX(2),lX(2)]);
            lY(1)=min([oldY(1),lY(1)]);
            lY(2)=max([oldY(2),lY(2)]);
            lZ(1)=min([oldZ(1),lZ(1)]);
            lZ(2)=max([oldZ(2),lZ(2)]);
            
            rX=lX(2)-lX(1);
            rY=lY(2)-lY(1);
            rZ=lZ(2)-lZ(1);
            maxRange=max([rX,rY,rZ]);
            
            %             ax.XLim=[ax.XLim(1)-(maxRange-rX)/2,ax.XLim(2)+(maxRange-rX)/2];
            %             ax.YLim=[ax.YLim(1)-(maxRange-rY)/2,ax.YLim(2)+(maxRange-rY)/2];
            %             ax.ZLim=[ax.ZLim(1)-(maxRange-rZ)/2,ax.ZLim(2)+(maxRange-rZ)/2];
            
            ax.XLim=[-(maxRange)/(1-lX(2)/lX(1)),(maxRange)/(1-lX(1)/lX(2))];
            ax.YLim=[-(maxRange)/(1-lY(2)/lY(1)),(maxRange)/(1-lY(1)/lY(2))];
            ax.ZLim=[-(maxRange)/(1-lZ(2)/lZ(1)),(maxRange)/(1-lZ(1)/lZ(2))];
            
            if nargin>5 && ~isempty(viewArgs)
                view(ax,viewArgs);
            end
        end
        
        function simulate(obj,q,delay,trace,drawArgs)
            if nargin<2 || isempty(q)
                q=zeros(size(obj.q));
            end
            
            if nargin<3 || isempty(delay)
                delay=0.001;
            end
            
            if nargin<4 || isempty(trace)
                trace=false;
            end
            
            if nargin<5 
                drawArgs={};
            end
            
            if nargin<5 || isempty(drawArgs)||isempty(drawArgs{2})
                ax=gca;
            end
            ax.NextPlot='replacechildren';
            
            n=size(q,2);
            for i=1:n
                obj.draw(q(:,i),drawArgs{:});
                pause(delay);
            end
        end
    end
    
    methods
        function set.q(obj,q)
            sz=size(q);
            
            obj.val_q=q;
            obj.val_d_q=sym('d_q',sz,'real');
            obj.val_dd_q=sym('dd_q',sz,'real');
            
            for i=1:sz(1)
                obj.val_d_q(i)=sym(strcat('d_',char(q(i))),'real');
                obj.val_dd_q(i)=sym(strcat('dd_',char(q(i))),'real');
            end
        end
        function value = get.q(obj)
            value=obj.val_q;
        end
        function value = get.d_q(obj)
            value=obj.val_d_q;
        end
        function value = get.dd_q(obj)
            value=obj.val_dd_q;
        end
    end
    
    methods(Static)
        function obj = fromH_TransChain(tForms,jointVars)
            obj = Kin_Model;
            obj = obj.init(numel(tForms));
            
            if nargin<2
                obj.q = sym('q',[obj.n_frames,1]);
            else
                jointVars=reshape(jointVars,numel(jointVars),1);
                obj.q = jointVars;
            end
            assume(obj.q, 'real');
            obj.q_limit=cell([size(obj.q,1),2]);
            for i = 1:obj.n_frames
                obj.T_matrices{i,i+1}=tForms{i};
            end

        end
        
        function obj=fromDH(DH_Matrix,jointVars)
            obj = Kin_Model;
            obj = obj.init(size(DH_Matrix,1));
            if nargin<2
                obj.q = sym('q',[obj.n_frames,1]);
            else
                jointVars=reshape(jointVars,numel(jointVars),1);
                obj.q = jointVars;
            end
            assume(obj.q, 'real');
            obj.q_limit=cell([size(obj.q,1),2]);
            for i = 1:obj.n_frames
                obj.T_matrices{i,i+1}=H_Trans.fromDH(DH_Matrix(i,:));
            end
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
            DH = sym('DH',[h w],'real');
            for r = 1:h
                for c = 1:w
                    DH(r,c) = sym(M{r,c},'real');
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
                obj.T_matrices{i,i} = H_Trans();
            end
        end
        
        function obj = calculateTransformation(obj,from,to)
            c=to+1;
            r=from+1;
            if c>(r+1)
                temp=eye(4);
                for i=r:(c-1)
                    temp=temp*obj.T(i-1,i).H;
                end
                obj.T_matrices{r,c} = H_Trans(vpa(temp));
            elseif c<r
                obj.T_matrices{r,c} = obj.T(c-1,r-1).inv;
            end
        end
        
        function obj = calculatePlot(obj)
            h=sym(zeros(4,4,obj.n_frames));
            frames=cell(obj.n_frames,1);
            
            for i=0:obj.n_frames
                frames{i+1}=obj.T(0,i);
                h(:,:,i+1)=frames{i+1}.H;
            end
            obj.trans_func=createFunction(h,obj.q);
        end
        
        function obj = calculateForwardKinematics(obj)
            T=obj.T(0,obj.n_frames);
            obj.f_kin_func = T.getFunction(obj.q);
        end
        
        function obj = calculateJacobian(obj)
            J=obj.T(0,obj.n_frames).getJacobian(obj.q);
            obj.J_func=createFunction(J,obj.q);
        end
        function obj = calculatePesudoInvJacobian(obj)
            J=obj.J;
            
            lambda_inv_J = J.'/(J*J.' + (obj.sym_lambda^2).*eye(size(J,1)));
            
            obj.inv_JLambda_func = createFunction(lambda_inv_J,{obj.q,sym('lambda')});
            
            inv_J0 = vpa(subs(lambda_inv_J,obj.sym_lambda,0));
            obj.inv_J0_func=createFunction(inv_J0,obj.q);
        end
        function obj = calculateJacobianDot(obj)
            d_J=obj.diffT(obj.J);
            obj.d_J_func=createFunction(d_J,{obj.q,obj.d_q});
        end
        
        function obj = calculateAnalyticJacobian(obj)
            Ja=obj.T(0,obj.n_frames).JgToJa(obj.J);
            obj.Ja_func=createFunction(Ja,obj.q);
        end
        function obj = calculatePesudoInvAnalyticJacobian(obj)
            Ja=obj.Ja;
            
            lambda_inv_Ja = Ja.'/(Ja*Ja.' + (obj.sym_lambda^2).*eye(size(Ja,1)));
            
            obj.inv_JaLambda_func = createFunction(lambda_inv_Ja,{obj.q,sym('lambda')});
            
            inv_Ja0 = vpa(subs(lambda_inv_Ja,obj.sym_lambda,0));
            obj.inv_Ja0_func=createFunction(inv_Ja0,obj.q);
        end
        function obj = calculateAnalyticJacobianDot(obj)
            d_Ja=obj.diffT(obj.Ja);
            obj.d_Ja_func=createFunction(d_Ja,{obj.q,obj.d_q});
        end
        
        function equ=diffT(obj,equ)
            sz=size(obj.q);
            
            t=sym('t','real');
            
            q_t=sym('q_t',sz,'real');
            d_q_t=sym('d_q_t',sz,'real');
            dd_q_t=sym('dd_q_t',sz,'real');
            
            for i=1:sz
                q_t(i) = sym(strcat(char(obj.q(i)),'(',char(t),')'),'real');
                d_q_t(i) = diff(q_t(i),t);
                dd_q_t(i) = diff(q_t(i),t,2);
            end
            
            equ=subs(equ,obj.q,q_t);
            equ=subs(equ,obj.d_q,d_q_t);
            equ=subs(equ,obj.dd_q,dd_q_t);
            
            equ=diff(equ,t);
            
            equ=subs(equ,dd_q_t,obj.dd_q);
            equ=subs(equ,d_q_t,obj.d_q);
            equ=subs(equ,q_t,obj.q);
        end
    end
    
end