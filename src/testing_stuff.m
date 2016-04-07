clear all; close all; clc;

syms t1 t2 t3
assume(t1,'real');
assume(t2,'real');
assume(t3,'real');
H=H_Trans();
H.Euler=[t1;t2;t3]

H.Euler
simplify(H.Euler)
H.Euler=H.Euler;
H.Euler=H.Euler;
simplify(H.Euler)

q=sym('q',[3,1],'real');
disp('Calculating Kinematics');
kin=Kin_Model.fromDH([q(1),1,.1,-pi/2;
                      q(2),0,1,0;
                      q(3),0,1,0],q);
                  
% test1=kin.forward_kin([0;0;0]);
% test2=kin.forward_kin([pi/4;pi/3;3*pi/4]);
% test3=kin.forward_kin([1;1;1]);
% 
% kin.inverse_kin(test1)
% kin.inverse_kin(test2)
% kin.inverse_kin(test3)

H=H_Trans(kin.forward_kin);

Jg=H.getJacobian(kin.q);
Ja=H.getAnalyticJacobian(kin.q);



