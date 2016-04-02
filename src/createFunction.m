function func = createFunction(expr,input)
    function vars = extractVars(in)
        vars=[];
        in=reshape(in,numel(in),1);
        if iscell(in)
            for i=1:numel(in)
                in{i}=reshape(in{i},numel(in{i}),1);
                vars=union(vars,in{i});
            end
        else
            vars=reshape(in,numel(in),1);
        end
    end
    vars=extractVars(input);

    if size(setdiff(symvar(expr),symvar(vars)))>0
        expr=vpa(expr);
        func = @(varargin)vpa(subs(expr,vars,extractVars(varargin)));
    else
        if iscell(input)
            func = matlabFunction(expr,'Vars',input);
        else
            func = matlabFunction(expr,'Vars',{input});
        end
    end
end

