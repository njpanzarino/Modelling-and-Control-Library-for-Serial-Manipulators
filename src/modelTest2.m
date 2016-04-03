close all; clear all; clc;

%% Modelling

q=sym('q',[3,1],'real');

disp('Calculating Kinematics');
kin=Kin_Model.fromDH([q(1),1,.1,-pi/2;
                      q(2),0,1,0;
                      q(3),0,1,0],q);

model=Dyn_Model(kin);

model=model.add(H_Trans([-0.5,0,0]),5,10,1);
model=model.add(H_Trans([-0.5,0,0]),5,10,2);
model=model.add(H_Trans([-0.5,0,0]),5,10,3);

model.g_dir=[0;0;-1];

disp('Calculating Dynamics');
model = model.calculateDynamics();

%% Simulation
t=model.t;
q_d = [t/8;0.2;sin(2*t)];

x0=[ 0      ,0;
    -0.5    ,0.1;
     0.2    ,0.1];

Kp=[10;10;10];
Kv=[10;10;10];

plan1=Planner.fromSym(q_d(1));
plan2=Planner.trapezoid(x0(2,:),[q_d(2),0],1,1);
plan3=Planner.fromSym(q_d(3));
plan=Planner.join({plan1,plan2,plan3});

control=Controller.ComputedTorque(model,Kp,Kv);
noise=@(a,t)0.02*randn(size(a));

options = odeset('RelTol',1e-4,'AbsTol',1e-4.*ones(numel(model.q)*2,1));

model.simulate({plan,control,noise},0:.02:10,x0,options,{0,[],{'linewidth',2},[1,1,1]});

