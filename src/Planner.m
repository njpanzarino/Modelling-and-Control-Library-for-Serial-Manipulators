classdef Planner
% Planner is used to create paths/trajectories
%   Each static function in the class should have the following format:
%       function d_func = XXXXXX(args)
%   
%   d_func should be an anonymous function of the form:
%       function desired = d_func(t)
%
%   desired should be of the form [q,d_q,dd_q]

    methods(Static)
        function d_func = fromSym(sym_q,sym_t)
            if nargin<2
                vars=symvar(sym(sym_q));
                if numel(vars)>0
                    sym_t=vars(1);
                else
                    sym_t=sym('t','real');
                end
            end
            
            desired_t(:,1)=sym_q;
            desired_t(:,2)=diff(sym_q,sym_t);
            desired_t(:,3)=diff(desired_t(:,2),sym_t);

            d_func=matlabFunction(desired_t,'Vars',sym_t);
        end
        
        function d_func = fromList(desired_list,t0,step)
            n=numel(desired_list);
            
            d_func=@fromList;
            function desired = fromList(t)
                index=(t-t0)/step;
                i1=floor(index);
                i2=ceil(index);
                if i1==i2
                    desired=desired_list(index);
                elseif i2>n
                    desired=desired_list(n);
                else
                    desired=Planner.interpolate(t,i1*step,desired_list(i1),i2*step,desired_list(i2));
                end
            end
        end
        
        function d_func = spline(desired_list,t_range)
            n=numel(desired_list);
            sz=size(desired_list{1});
            funcs=cell(sz(1),1);
%             t_range=t0:step:n*step;
            
            for q_i=1:sz(1)
                q_d=zeros(n,1);
                
                for i=1:n
                    q_d(i)=desired_list{i}(q_i,1);
                end
                
                if size(desired_list{1},2)>1
                    vi=desired_list{1}(q_i,2);
                else
                    vi=0;
                end
                
                if size(desired_list{n},2)>1
                    vf=desired_list{n}(q_i,2);
                else
                    vf=0;
                end
                
                pp=spline(t_range,[vi;q_d;vf]);
                d_pp=deriv_pp(pp);
                dd_pp=deriv_pp(d_pp);
                funcs{q_i}=@(t)[ppval(pp,t),ppval(d_pp,t),ppval(dd_pp,t)];
            end
            
            d_func = Planner.join(funcs);
            
            function d_pp=deriv_pp(pp)
                [breaks,coefs,l,k,d] = unmkpp(pp);
                % make the polynomial that describes the derivative
                d_pp = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
            end
        end
        
        function [d_func, sym_a] = trajectory(final, initial, tf, t0)
            if nargin<4 || isempty(t0)
                t0=0;
            end
            
            sym_t=sym('t','real');

            [q0,d_q0]=Planner.interpretInput(initial);
            [qf,d_qf]=Planner.interpretInput(final);

            sym_q=sym(zeros(size(q0)));
            sym_a=sym(zeros(size(q0,1),4));

            for i=0:size(sym_q,1);
                [sym_q(i),sym_a(i,:)]=trajectory_single(q0(i),d_q0(i),qf(i),d_qf(i),t0,tf,sym_t);
            end
            
            d_func = matlabFunction(sym_q,'Vars',t);
            
            function [sym_q,sym_a] = trajectory_single(q0,v0,qf,vf,t0,tf,t)
                b=[q0;v0;qf;vf];

                M=[1,t0, t0^2,t0^3;
                    0,1,2*t0,3*t0^2;
                    1,tf,tf^2,tf^3;
                    0,1,2*tf,3*tf^2];

                a=M\b;

                sym_a=a.';
                sym_q=a(1)+a(2)*t+a(3)*t^2+a(4)*t^3;
            end
        end
        
        function [d_func,tf,tm,tu] = trapezoid(initial,final,v_max,a_max)
            sz=size(initial,1);
            
            if numel(v_max)==1
                v_max=v_max.*ones(sz,1);
            end
            if numel(a_max)==1
                a_max=a_max.*ones(sz,1);
            end
            
            t=sym('t','real');
            
            if sz==1
                [d_func,tf,tm,tu] = trapezoid_single(final,initial,v_max,a_max);
            else
                tf=zeros(size(sz,1));
                tm=zeros(size(sz,1));
                tu=zeros(size(sz,1));
                funcs=cell(sz,1);
                for i=1:sz
                    [funcs{i},tf(i),tm(i),tu(i)] = trapezoid_single(final(i,:),initial(i,:),v_max(i),a_max(i));
                end
                d_func=Planner.join(funcs);
            end
            
            function [desired,tf,tm,tu] = trapezoid_single(final,initial,v_max,a_max)
                v_max=abs(v_max);
                a_max=abs(a_max);
            
                [qi,vi]=Planner.interpretInput(initial);
                [qf,vf]=Planner.interpretInput(final);
                
                if abs(vi)>v_max
                    error('vi > v_max');
                end
                if abs(vf)>v_max
                    error('vf > v_max');
                end

                dir=sign(qf-qi);
                vel=dir*v_max;
                acc=dir*a_max;
                
                tc = -(qi - qf + (vel^2 - vf^2)/(2*acc) + (vel^2 - vi^2)/(2*acc))/vel;
            
                if tc<0
                    vel=dir*abs((2^(1/2)*(vf^2 + vi^2 + 2*acc*qf - 2*acc*qi)^(1/2))/2);
                    tc=0;
                end

                tu=(vel-vi)/acc;
                td=(vf-vel)/-acc;
                tm=tu+tc;
                tf=tu+tc+td;

                q1=(vel^2-vi^2)/(2*acc);
                q2=q1+vel*tc;

                sym_u=qi+vi*t+(1/2)*acc*t^2;
                sym_c=qi+q1+vel*(t-tu);
                sym_d=qi+q2+vel*(t-tm)-(1/2)*acc*(t-tm)^2;
                sym_b=qf+vf*(t-tf);
                
%                 sympref('HeavisideAtOrigin',1);
%                 heaviside(sym(0))
                
                f_u=matlabFunction([sym_u,diff(sym_u,t),diff(sym_u,t,2)],'Vars',t);
                f_c=matlabFunction([sym_c,diff(sym_c,t),diff(sym_c,t,2)],'Vars',t);
                f_d=matlabFunction([sym_d,diff(sym_d,t),diff(sym_d,t,2)],'Vars',t);
                f_b=matlabFunction([sym_b,diff(sym_b,t),diff(sym_b,t,2)],'Vars',t);
                
                desired = @trapezoid_single;
                
                function desired = trapezoid_single(t)
                    if t<tu
                        desired=f_u(t);
                    elseif t<(tm)
                        desired=f_c(t);
                    elseif t<(tf)
                        desired=f_d(t);
                    else
                        desired=f_b(t);
                    end
                end
            end
        end
        
        function d_func = toJointSpace(d_func,equ_q,vars,t_range)
            equ_q=equ_q(:,1);
            vars=reshape(vars,[],1);
            
            map_pos=mapPos;
            map_vel=mapVel;
            
%             t_range=t0:t_step:tf;
            n=numel(t_range);
            desired=cell(n,1);
            for i=1:n
                workspace_d=d_func(t_range(i));
                desired{i}=map_pos(workspace_d(:,1));
                if i==1||i==n
                    desired{i}=[desired{i},map_vel(workspace_d(:,1:2))];
                end
            end
            
            d_func=Planner.spline(desired,t_range);
            
            function func = mapPos
                sym_q=equ_q;
                
                func=matlabFunction(sym_q,'Vars',{vars});
            end
            function func = mapVel
                sym_q=equ_q;
                [sym_d_q,d_vars,~]=Planner.diffT(sym_q,vars);
%                 [sym_dd_q,~,~]=Planner.diffT(sym_d_q,vars);
                
                func=matlabFunction(sym_d_q,'Vars',{[vars,d_vars]});
            end
        end
        
        function d_func = join(funcs)
            n=numel(funcs);
%             funcs=varargin{:};
            
            desired=zeros(0,3);
            for t=1:n
                desired=[desired;funcs{t}(0)];
            end
            sz_d=size(desired);
            
            d_func=@join;
            function desired = join(t)
                desired=zeros(sz_d);
                for i=1:n
                    g=funcs{i}(t);
                    for r=1:size(g,1)
                        desired(i+r-1,:)=g(r,:);
                    end
                end
            end
        end
    end
    
    methods(Static, Access=private)
        function [q,d_q,dd_q] = interpretInput(input)
            sz=size(input);
            
            q=input(:,1);
            
            if sz(2)>1
                d_q=input(:,2);
            else
                d_q=zeros(size(q));
            end
            
            if sz(2)>2
                dd_q=input(:,3);
            else
                dd_q=zeros(size(q));
            end
        end
        
        function q = interpolate(t,t1,q1,t2,q2)
            q=(q2-q1).*((t-t1)/(t2-t1))+q1;
        end
        
        function [equ,d_q,dd_q]=diffT(equ,vars,t)
            sz=size(vars);
            if nargin<3||isempty(t)
                t=sym('t','real');
            end
            
            d_q=sym('d_q',sz);
            dd_q=sym('dd_q',sz);
            
            q_t=sym('q_t',sz);
            d_q_t=sym('d_q_t',sz);
            dd_q_t=sym('dd_q_t',sz);
            for i=1:sz
                d_q(i)=sym(strcat('d_',char(vars(i))),'real');
                dd_q(i)=sym(strcat('dd_',char(vars(i))),'real');
                
                q_t(i) = sym(strcat(char(vars(i)),'(',char(t),')'),'real');
                d_q_t(i) = diff(q_t(i),t);
                dd_q_t(i) = diff(q_t(i),t,2);
            end
            
            equ=subs(equ,vars,q_t);
            equ=subs(equ,d_q,d_q_t);
            equ=subs(equ,dd_q,dd_q_t);
            
            equ=diff(equ,t);
            
            equ=subs(equ,dd_q_t,dd_q);
            equ=subs(equ,d_q_t,d_q);
            equ=subs(equ,q_t,vars);
        end
        
    end
end

