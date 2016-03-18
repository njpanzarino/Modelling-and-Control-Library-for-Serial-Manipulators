close all; clear all; clc;

%model=Kin_Model.fromPrompt()

q=sym('potato',[3,1]);
d=sym('d',[3,1]);
kin=Kin_Model.fromDH([q(1),3,0,pi/2;
                        q(2),0,2,0;
                        q(3),0,1,0],q);

kin.q

H1= kin.forward_kin([.1;.1;.1])
H2= kin.forward_kin([pi/3;pi/4;-pi/3])

model=Dyn_Model(kin);

kin.inverse_kin(H2)
