close all; clear all; clc;

%% Kinematics

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

desired=[0;1;1;0;0;0]
test=kin.forward_kin(desired)
result=kin.inverse_kin(test)

kin.forward_kin(result)

kin.draw(desired);
figure
kin.draw(result);

%% Dynamics

model=Dyn_Model(kin);

% model=model.add(H_Trans([-0.5,0,0]),5,10,1);
% model=model.add(H_Trans([-0.5,0,0]),5,10,2);
% model=model.add(H_Trans([-0.5,0,0]),5,10,4);
% model=model.add(H_Trans,1,diag([1,1,1]),8);

% model.b=0.1.*ones(size(model.b));

model.g_dir=[0;0;-1];

disp('Calculating Dynamics');
model = model.calculateDynamics();

%% Simulation
t=model.t;
q_d = [t/6;0;0;t;2*t;3*t];

x0=[ 0      ,0;
     0      ,0;
     0      ,0;
     0      ,0;
     0      ,0;
     0      ,0];

Kp=[10;10;10;1000;1000;1000];
Kv=[10;10;10;100;100;100];
Ki=[1;1;1;1;1;1];

plan=Planner.fromSym(q_d);

control=Controller.ComputedTorque(@model.inverse_dyn,Kp,Kv);
% control=Controller.PID(Kp,Ki,Kv);
% noise=@(a,t)0.02*randn(size(a));
% noise=[];

options = odeset('RelTol',1e-2,'AbsTol',1e-2.*ones(numel(model.q)*2,1));

figure
D=evalf(plan,(0:.02:10).');
kin.simulate(D(:,:,1).');

% model.simulate({plan,control,noise},0:.1:10,x0,options,{0.01,false,{0,[],{'linewidth',2},[1,1,1]}});

