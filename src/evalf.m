function vals = evalf( func, inputs )
%EVALF Summary of this function goes here

    sz_in=size(inputs);

    val1=func(inputs(1,:));
    n=numel(val1);
    sz_out=size(val1);

    vals=zeros([sz_in(1),sz_out]);

    vals(1,:)=reshape(val1,[1,n]);
    for ti=2:sz_in(1)
        current=func(inputs(ti,:));
        vals(ti,:)=reshape(current,[1,n]);
    end

end

