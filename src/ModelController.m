function [ control_output ] = ModelController(desired, actual, type, args)
%ModelController allows access to a variety of different control schemes which
%utilize the Dyn_Model object for calculations
%   desired: desired joint variables. 3xn matrix in the form [q,d_q,dd_q]
%       where n is the number of joint variables
%   actual: same format as desired, but with current joint values (as
%       measured)
%   type: string representation of desired control method
%   varargin: the remaining inputs will be passed to the function of the
%       chosen control scheme
    switch(type)
        case 'ComputedTorque'
            control_output=ComputedTorque(desired,actual,args{:});
        case 'PD'
            control_output=PD(desired,actual,args{:});
        case 'P'
            control_output=P(desired,actual,args{:});
        otherwise
            control_output=P(desired,actual,args{:});
    end
end

function tau=ComputedTorque(desired,actual,model,Kp,Kv)
    q=actual(:,1);
    d_q=actual(:,2);

    q_d=desired(:,1);
    d_q_d=desired(:,2);
    dd_q_d=desired(:,3);

    e1=q-q_d;
    e2=d_q-d_q_d;
    
    tau=model.inverse_dyn(q,d_q,dd_q_d-Kp.*e1-Kv.*e2);
end

function tau=PD(desired,actual,Kp,Kd)
    q=actual(:,1);
    d_q=actual(:,2);

    q_d=desired(:,1);
    d_q_d=desired(:,2);

    e1=q-q_d;
    e2=d_q-d_q_d;
    
    tau=Kp.*e1+Kd.*e2;
end

function tau=P(desired,actual,Kp)
    q=actual(:,1);

    q_d=desired(:,1);

    e1=q-q_d;
    
    tau=Kp.*e1;
end