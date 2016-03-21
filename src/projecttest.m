close all; clear all; clc;

%model=Kin_Model.fromPrompt()

q=sym('potato',[3,1]);
d=sym('d',[3,1]);
kin=Kin_Model.fromDH([q(1),3,0,pi/2;
                        q(2),0,2,0;
                        q(3),0,1,0],q);

kin.q

H1= kin.forward_kin([.1;.1;.1]);
H2= kin.forward_kin([pi/3;pi/4;-pi/3]);

model=Dyn_Model(kin);

model=model.add(H_Trans([2,3,4]),1,.5,1);
model=model.add(H_Trans([3,4,5]),2,.4,2);
model=model.add(H_Trans([4,5,6]),3,.3,3);

testEqu=sin(model.q(1)+model.q(2))+model.q(3);

model.diff(testEqu);

model.calculateDynamics();
