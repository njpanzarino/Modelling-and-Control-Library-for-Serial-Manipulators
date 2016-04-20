function func = createFunction(expr,input)
%CREATEFUNCTION creates an anonymous function from the given symbolic expression.
%   expr: a symbolic expression to create a function for
%   input: the variables contained in the expression, in the format desired
%       for the anonymous function. If there are variables in the expression
%       that are not in the input, the subs function will be used (which is
%       slower than matlabFunction)
    function vars = extractVars(in)
        vars=[];
        in=reshape(in,numel(in),1);
        if iscell(in)
            for i=1:numel(in)
                in{i}=reshape(in{i},numel(in{i}),1);
                vars=union(in{i},vars);
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

