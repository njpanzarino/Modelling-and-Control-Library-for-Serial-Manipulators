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
Kp=10*ones(size(q_d));
Kv=10*ones(size(q_d));

c=Controller.ComputedTorque(model,Kp,Kv);
noise=@(a,t)0.02*randn(size(a));

x0=[[-0.5;0.2],[0.1;0.1]];
options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4]);

model.simulate({q_d,c,noise},0:.02:10,x0,options,{0,[],{'linewidth',2},[0,0,1]});
% 
% [ode,getTorques]=model.ODE(q_d,c);
% 
% disp('Beginning Simulation');
% 
% [T,Y]=ode45(ode,[0 10],x0);
% [torques,times]=getTorques();
% 
% disp('Simulation Complete');
% 
% a=[Y(:,1),Y(:,2)].';
% kin.simulate(a,0,[],{'linewidth',2});
% 
% plot(T, Y(:,1),'r-');
% hold on
% plot(T, q_d(1)*ones(size(T,1),1),'b-');
% figure('Name','Theta_2 under Computed Torque Control');
% plot(T, Y(:,2),'r--');
% hold on
% plot(T, sin(2*T),'b-');
% 
% figure('Name','Input_Computed Torque Control');
% plot(times, torques(1:size(times,1),1),'-' );
% hold on
% plot(times, torques(1:size(times,1),2),'r--');
% 
