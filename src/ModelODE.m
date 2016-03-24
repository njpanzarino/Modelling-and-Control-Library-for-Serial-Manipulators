function [ ode_func, tau_func ] = ModelODE(model,sym_t,sym_x, type, args)
%ModelODE Summary of this function goes here
%   Detailed explanation goes here

    desired_t(:,1)=sym_x;
    desired_t(:,2)=diff(sym_x,sym_t);
    desired_t(:,3)=diff(desired_t(:,2),sym_t);
    
    func_x=matlabFunction(desired_t);
    
    save_tau=[];
    save_time=[];
    
    ode_func = @constructFunc;
    tau_func = @getTau;
    
    function d_x = constructFunc(t,x)
        sz=size(x);
        desired=func_x(t);
%         actual=x;
        actual=reshape(x,[sz(1)/2,2]);
        tau=ModelController(desired, actual, type, args);
        
        save_tau=[save_tau;tau.'];
        save_time=[save_time;t];
        
        d_x(:,1)=actual(:,2);
        d_x(:,2)=model.forward_dyn(actual(:,1),actual(:,2),tau);
        
        d_x=reshape(d_x,sz);
    end

    function [torques,times] = getTau()
        torques=save_tau;
        times=save_time;
    end
end

