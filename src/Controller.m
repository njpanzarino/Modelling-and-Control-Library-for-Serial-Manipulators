classdef Controller
% Controller allows access to various control methods
%   Each static function in the class should have the following format:
%       function c_func = XXXXXX(args)
%   
%   c_func should be an anonymous function of the form:
%       function tau = c_func(desired,actual,...)

    methods(Static)
        function tau = ComputedTorque(model,Kp,Kv)

            if nargin<2 || isempty(Kp)
                Kp=ones(size(model.q));
            end

            if nargin<2 || isempty(Kv)
                Kv=ones(size(model.q));
            end

            tau = @ComputedTorque;

            function tau = ComputedTorque(desired,actual)
                q=actual(:,1);
                d_q=actual(:,2);

                q_d=desired(:,1);
                d_q_d=desired(:,2);
                dd_q_d=desired(:,3);

                e1=q-q_d;
                e2=d_q-d_q_d;

                tau=model.inverse_dyn(q,d_q,dd_q_d-Kp.*e1-Kv.*e2);
            end
        end

        function tau = PID(Kp,Ki,Kd)

            sum=zeros(size(Kp));
            tau = @PID;

            function tau = PID(desired,actual)
                q=actual(:,1);
                d_q=actual(:,2);

                q_d=desired(:,1);
                d_q_d=desired(:,2);

                e1=q-q_d;
                e2=d_q-d_q_d;
                
                sum=sum+e1;

                tau=Kp.*e1+Kd.*e2+Ki.*sum;
            end
        end
        
        function tau = PI(Kp,Ki)

            sum=zeros(size(Kp));
            tau = @PI;

            function tau = PI(desired,actual)
                q=actual(:,1);

                q_d=desired(:,1);

                e1=q-q_d;
                
                sum=sum+e1;

                tau=Kp.*e1+Ki.*sum;
            end
        end
        
        function tau = PD(Kp,Kd)

            tau = @PD;

            function tau = PD(desired,actual)
                q=actual(:,1);
                d_q=actual(:,2);

                q_d=desired(:,1);
                d_q_d=desired(:,2);

                e1=q-q_d;
                e2=d_q-d_q_d;

                tau=Kp.*e1+Kd.*e2;
            end
        end

        function tau = P(Kp)

            tau = @P;

            function tau = P(desired,actual)
                q=actual(:,1);
                q_d=desired(:,1);

                e1=q-q_d;

                tau=Kp.*e1;
            end
        end
    end
end

