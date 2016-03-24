close all; clear all; clc;

%model=Kin_Model.fromPrompt()

q=sym('q',[2,1]);
d=sym('d',[3,1]);
kin=Kin_Model.fromDH([%q(1),3,0,pi/2;
                      q(1),0,1,0;
                      q(2),0,1,0],q);

% kin.q

% H1= kin.forward_kin([.1;.1;.1]);
% H2= kin.forward_kin([pi/3;pi/4;-pi/3]);

model=Dyn_Model(kin);

% model=model.add(H_Trans([0,0,-1.5]),3,0,1);
model=model.add(H_Trans([-0.5,0,0]),5,10,1);
model=model.add(H_Trans([-0.5,0,0]),5,10,2);

model.g_dir=[0;-1;0];

model = model.calculateDynamics();

syms t;

d_q = [0.2;sin(2*t)];
Kp=10*ones(size(d_q));
Kv=10*ones(size(d_q));

x0=[[-0.5;0.2],[0.1;0.1]];
options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4]);

[ode,getTorques]=ModelODE(model,t,d_q,'ComputedTorque',{model,Kp,Kv});

disp('Beginning Simulation');

[T,Y]=ode45(ode,[0 10],x0);
[torques,times]=getTorques();

disp('Simulation Complete');

plot(T, Y(:,1),'r-');
hold on
plot(T, d_q(1)*ones(size(T,1),1),'b-');
figure('Name','Theta_2 under Computed Torque Control');
plot(T, Y(:,2),'r--');
hold on
plot(T, sin(2*T),'b-');

figure('Name','Input_Computed Torque Control');
plot(times, torques(1:size(times,1),1),'-' );
hold on
plot(times, torques(1:size(times,1),2),'r--');

