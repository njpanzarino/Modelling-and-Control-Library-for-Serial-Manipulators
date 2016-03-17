close all; clear all; clc;

%model=Kin_Model.fromPrompt()

q=sym('potato',[3,1]);
d=sym('d',[3,1]);
kin=Kin_Model.fromDH([q(1),d(1),0,pi/2;
                        q(2),0,d(2),0;
                        q(3),0,d(3),0],q);

kin.J

kin.q

model=Dyn_Model(kin)