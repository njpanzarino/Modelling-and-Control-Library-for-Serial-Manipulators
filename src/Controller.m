classdef Controller
% Controller allows access to various control methods
%   Each static function in the class should have the following format:
%       function c_func = XXXXXX(args)
%   
%   c_func should be an anonymous function of the form:
%       function tau = c_func(desired,actual,...)

    methods(Static)
        function tau = ComputedTorque(inv_dyn_func,Kp,Kv)

            tau = @ComputedTorque;

            function tau = ComputedTorque(desired,actual)
                [q,d_q]=Controller.interpretInput(actual);
                [q_d,d_q_d,dd_q_d]=Controller.interpretInput(desired);

                ep=q_d-q;
                ev=d_q_d-d_q;

                tau=inv_dyn_func(q,d_q,dd_q_d+Kp.*ep+Kv.*ev);
            end
        end

        function tau = PID(Kp,Ki,Kd)

            sum=zeros(size(Kp));
            tau = @PID;

            function tau = PID(desired,actual)
                [q,d_q]=Controller.interpretInput(actual);
                [q_d,d_q_d]=Controller.interpretInput(desired);;

                ep=q_d-q;
                ev=d_q_d-d_q;
                
                sum=sum+ep;

                tau=Kp.*ep+Kd.*ev+Ki.*sum;
            end
        end
        
        function tau = PI(Kp,Ki)

            sum=zeros(size(Kp));
            tau = @PI;

            function tau = PI(desired,actual)
                [q]=Controller.interpretInput(actual);
                [q_d]=Controller.interpretInput(desired);

                ep=q_d-q;
                
                sum=sum+ep;

                tau=Kp.*ep+Ki.*sum;
            end
        end
        
        function tau = PD(Kp,Kd)

            tau = @PD;

            function tau = PD(desired,actual)
                [q,d_q]=Controller.interpretInput(actual);
                [q_d,d_q_d]=Controller.interpretInput(desired);

                ep=q_d-q;
                ev=d_q_d-d_q;

                tau=Kp.*ep+Kd.*ev;
            end
        end

        function tau = P(Kp)

            tau = @P;

            function tau = P(desired,actual)
                [q]=Controller.interpretInput(actual);
                [q_d]=Controller.interpretInput(desired);

                ep=q_d-q;

                tau=Kp.*ep;
            end
        end
        
        function tau = I(Ki)

            sum=zeros(size(Kp));
            tau = @I;

            function tau = I(desired,actual)
                [q]=Controller.interpretInput(actual);
                [q_d]=Controller.interpretInput(desired);

                ep=q_d-q;
                sum=sum+ep;

                tau=Ki.*sum;
            end
        end
    end
    
    methods(Static, Access=private)
        function [q,d_q,dd_q] = interpretInput(input)
            sz=size(input);
            
            q=input(:,1);
            
            if sz(2)>1
                d_q=input(:,2);
            else
                d_q=zeros(size(q));
            end
            
            if sz(2)>2
                dd_q=input(:,3);
            else
                dd_q=zeros(size(q));
            end
        end
    end
end

