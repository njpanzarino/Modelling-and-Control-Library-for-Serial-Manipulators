classdef Dyn_Model
    %Dyn_Model Representation of a Dynamical Model of a Serial Robotic
    %Manipulator
    %   Detailed explanation goes here
    
    properties(SetAccess = private)
        kin %Kinematic Model of a serial Robot Manipulator
        
        I %Inertias relative to base frame
        
        m_val %mass values
        m_loc %mass locations relative to base frame
        
        k %stiffness for each joint 
        b %damping for each joint
        
    end
    
    properties(Dependent, SetAccess=private)
        q    %joint variables
        d_q  %1st time derivative of joint variables
        dd_q %2nd time derivative of joint variables
    end
    
    properties(Access=private)
        val_q    %joint variables
        val_d_q  %1st time derivative of joint variables
        val_dd_q %2nd time derivative of joint variables
    end
    
    methods
        function obj = Dyn_Model(model)
            if nargin>0
                obj=Dyn_Model.fromKin_Model(model);
            end
        end
        
        function obj = set.q(obj,q)
            n=size(q);
            
            obj.val_q=q;
            obj.val_d_q=sym('d_q',n);
            obj.val_dd_q=sym('dd_q',n);
            
            for i=1:n
                obj.val_d_q(i)=sym(strcat('d_',char(q(i))));
                obj.val_dd_q(i)=sym(strcat('dd_',char(q(i))));
            end
        end
        function value = get.q(obj)
            value=obj.val_q;
        end
        function value = get.d_q(obj)
            value=obj.val_d_q;
        end
        function value = get.dd_q(obj)
            value=obj.val_dd_q;
        end

    end
    
    methods(Static)
        function obj = fromKin_Model(model)
            obj = Dyn_Model;
			obj.kin=model;
            obj.q=obj.kin.q;
        end
    end
end

