close all; clear all; clc;

%% Modelling

%model=Kin_Model.fromPrompt()

q=sym('q',[2,1],'real');

kin=Kin_Model.fromDH([q(1),0,1,0;
                      q(2),0,1,0],q);

model=Dyn_Model(kin);

model=model.add(H_Trans([-0.5,0,0]),5,10,1);
model=model.add(H_Trans([-0.5,0,0]),5,10,2);

model.g_dir=[0;-1;0];

model = model.calculateDynamics();

%% Simulation

t=model.t;

q_d = [0.2;sin(2*t)];
Kp=[10;10];
Kv=[10;10];

plan=Planner.fromSym(q_d);
control=Controller.ComputedTorque(model,Kp,Kv);
noise=@(a,t)0.02*randn(size(a));

x0=[[-0.5;0.2],[0.1;0.1]];
options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4]);

model.simulate({plan,control,noise},0:.02:10,x0,options,{0,[],{'linewidth',2},[0,0,1]});

