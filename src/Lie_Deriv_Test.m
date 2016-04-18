function output_args = Lie_Deriv_Test( )

    a=sym('a');
    x=sym('x',[2,1],'real');

    function d = desired(~,actual)
        d=actual;
    end

    function func = control(k)
        func=@control;
        T1=x(1);
        T2=L(f(x),T1,x);

        v=-[T1,T2]*k;
        sym_tau = 1/(L(g,T2,x))*(v-L(f(x),T2,x));
        function tau = control(~,actual)
            tau=[0;subs(sym_tau,x,actual.')];
        end
    end

    function r = response(actual,tau)
        r=f(actual)+g.*tau;
    end

    function lie_deriv = L(f,h,x)        
        lie_deriv = jacobian(h,x)*f;
    end
    
    a=3;
    k=[5;5];
    x0=[.1 .4];
    
    f = @(x)[a*sin(x(2)); -(x(1))^2];
    g = [0; 1];
    
    ode=CreateODE(@desired,control(k),@response);
    
    ode45(ode,[0 5],x0)
end

