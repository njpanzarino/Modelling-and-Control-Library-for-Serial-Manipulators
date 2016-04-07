close all; clear all; clc;

%% Modelling

q=sym('q',[2,1],'real');

disp('Calculating Kinematics');
kin=Kin_Model.fromDH([q(1),0,1/2,0;
                      q(2),0,1,0],q);

model=Dyn_Model(kin);

model=model.add(H_Trans(),8,0,1);
model=model.add(H_Trans(),8,0,2);

model.g_dir=[0;-1;0];

disp('Calculating Dynamics');
model = model.calculateDynamics();

%% Simulation
t=model.t;

x0=[ pi/4    ,0;
     pi/2    ,0];

Kp=[10;100];
Kv=[10;100];

plan=Planner.fromSym(x0(:,1));

control=Controller.ComputedTorque(@model.inverse_dyn,Kp,Kv);
control=@(desired,actual)[0;1].*control(desired,actual);

% noise=[];
noise=@(a,t)0.02*randn(size(a));

options = odeset('RelTol',1e-4,'AbsTol',1e-4.*ones(numel(model.q)*2,1));

model.simulate({plan,control,noise},0:.02:10,x0,options,{0,[],{'linewidth',2},[0,0,1]});

