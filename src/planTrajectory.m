function [sym_q,sym_t,sym_a] = planTrajectory( desired, actual, tf, t0 )

if nargin<4 || isempty(t0)
    t0=0;
end

sym_t=sym('t','real');

q0=actual(:,1);
d_q0=actual(:,2);

qf=desired(:,1);
d_qf=desired(:,2);


sym_q=sym(zeros(size(q0)));
sym_a=sym(zeros(size(q0,1),4));

for i=0:size(sym_q,1);
    [sym_q(i),sym_a(i,:)]=plan(q0(i),d_q0(i),qf(i),d_qf(i),t0,tf,sym_t);
end

end

function [sym_q,sym_a] = plan(q0,v0,qf,vf,t0,tf,t)

    b=[q0;v0;qf;vf];

    M=[1,t0, t0^2,t0^3;
        0,1,2*t0,3*t0^2;
        1,tf,tf^2,tf^3;
        0,1,2*tf,3*tf^2];

    a=M\b;
    
    sym_a=a.';
    sym_q=a(1)+a(2)*t+a(3)*t^2+a(4)*t^3;
end