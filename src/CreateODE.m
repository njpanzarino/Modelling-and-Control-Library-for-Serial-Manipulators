function [ ode_func ] = CreateODE(desired_func,control_func,response_func,noise_func)
%CreateODE Returns derivative of state x at time t, given functions for
%desired state, controller output, and expected system response
%   desired_func: desired = desired_func(t)
%   control_func: tau = control_func(desired,actual)
%   response_func: diff(actual,t) = response_func(actual,tau)
%   noise_func: noise = noise_func(actual,tau)  **(Optional)**
%       noise = actual response - expected response
%       noise_func is an optional parameter which may be added to model
%           uncertainty in the system. Similar functionality can be
%           achieved by building noise into the response_fun instead of
%           providing a noise_func
%
%   actual = [q,d_q]
%   desired = [q_d,d_q_d,dd_q_d]
%   tau = column vector of joint torques/forces. same size as q
    
    if nargin>3 && ~isempty(noise_func)
        response=@(actual,tau)(response_func(actual,tau)+noise_func(actual,tau));
    else
        response=response_func;
    end
    
    ode_func = @ode;
    
    function [d_x,tau] = ode(t,x)
        sz=size(x);
        desired=desired_func(t);
        actual=reshape(x,[],2);
        
        tau=control_func(desired, actual);
        
        d_x=response(actual,tau);
        
        d_x=reshape(d_x,sz);
    end

end

