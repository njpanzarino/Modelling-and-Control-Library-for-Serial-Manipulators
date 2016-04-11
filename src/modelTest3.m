close all; clear all; clc;

%% Modelling

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
% 

result=kin.inverse_kin(test)

kin.forward_kin(result)

kin.draw(desired);
figure
kin.draw(result);

model=Dyn_Model(kin);

% model=model.add(H_Trans([-0.5,0,0]),5,10,1);
% model=model.add(H_Trans([-0.5,0,0]),5,10,2);
% model=model.add(H_Trans([-0.5,0,0]),5,10,3);

model.g_dir=[0;0;-1];

disp('Calculating Dynamics');
model = model.calculateDynamics();

%% Simulation
t=model.t;
q_d = [0;0;0;0;0;0];

x0=[ 0      ,0;
     0      ,0;
     0      ,0;
     0      ,0;
     0      ,0;
     0      ,0];

Kp=[10;10;10;10;10;10];
Kv=[10;10;10;10;10;10];

plan=Planner.fromSym(q_d);

control=Controller.ComputedTorque(@model.inverse_dyn,Kp,Kv);
noise=@(a,t)0.02*randn(size(a));

options = odeset('RelTol',1e-4,'AbsTol',1e-4.*ones(numel(model.q)*2,1));

model.simulate({plan,control,noise},0:.02:10,x0,options,{0,[],{'linewidth',2},[1,1,1]});

