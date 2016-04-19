%% Workspace Planning
clear all; clc;

t=sym('t','real');
t_range=0:0.02:10;

q=sym('q',[6,1],'real');

disp('Calculating Kinematics');
kin=Kin_Model.fromDH([q(1),1,0,pi/2;
                      q(2),-0.1,0,0;
                      0,0,1,0;
                      q(3),0.1,0,0;
                      0,0,1,pi/2;
                      q(4),0.1,0,-pi/2;
                      q(5),0.1,0,-pi/2;
                      q(6),0.1,0,pi/2],q);


r=1;
p0=[1;1;1];
z_t=r-r/t_range(numel(t_range))*t+p0(3);
r_t=r-(z_t-p0(3))^2;
p_d=[p0(1)+r_t*cos(10*t);p0(2)+r_t*sin(10*t);z_t];
x_d=[NaN,NaN,NaN,p_d(1);
     NaN,NaN,NaN,p_d(2);
     NaN,NaN,NaN,p_d(3);
     NaN,NaN,NaN,NaN];

 pts=subs(p_d,t_range);
plot3(pts(1,:),pts(2,:),pts(3,:));

plan=Planner.fromSym(H_Trans(x_d).H);
plan=Planner.toJointSpace_func(plan,t_range,[],@kin.inverse_kin);

D=evalf(plan,t_range.');
kin.simulate(D(:,:,1).',0.01,true,{0,[],{'linewidth',2},[1,1,1]});