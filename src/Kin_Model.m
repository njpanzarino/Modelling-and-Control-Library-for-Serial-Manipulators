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
        
        trans_func %Function for determining points/transformations for each link when plotting
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
            
            %baisc fsolve
            options = optimoptions(@fsolve,'Display','iter',...
                'Algorithm','trust-region-reflective',...
                'Jacobian','off');
            [value,~]=fsolve(@(q)(obj.forward_kin(q)-H),q0,options);
            
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
        
        function simulate(obj,q,drawFrames,ax,plotArgs,viewArgs)
            if nargin<2 || isempty(q)
                q=zeros(size(obj.q));
            end
            
            if nargin<3 || isempty(drawFrames)
                drawFrames=0;
            end
            
            if nargin<4 || isempty(ax)
                ax=gca;
            end
            ax.NextPlot='replacechildren';
            
            if nargin<5 || isempty(plotArgs)
                plotArgs={};
            end
            
            if nargin<6 || isempty(viewArgs)
                viewArgs=[];
            end
            
            if drawFrames>0
                p=0.01;
            else
                p=0.0001;
            end
            
            n=size(q,2);
            for i=1:n
                obj.draw(q(:,i),drawFrames,ax,plotArgs,viewArgs);
                pause(p);
            end
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
            
            obj=obj.calculateFromChain();
            obj=obj.calculatePlot();
            obj=obj.calculateForwardKinematics();
            obj=obj.calculateJacobian();
            obj=obj.calculatePesudoInvJacobian();
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
            obj=obj.calculateFromChain();
            obj=obj.calculatePlot();
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
        
    end
    
end