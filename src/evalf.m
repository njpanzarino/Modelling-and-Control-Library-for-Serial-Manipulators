function vals = evalf( func, varargin )
%EVALF Summary of this function goes here
    inputs=varargin;
    
    n_in=nargin-1;
    
    sz_in=cell(n_in,1);
    for i=1:n_in
        sz_in{i}=size(inputs{i});
    end

    in=getInputs(1);
    val1=func(in{:});
    n=numel(val1);
    sz_out=size(val1);

    vals=zeros([sz_in{1}(1),sz_out]);

    vals(1,:)=reshape(val1,[1,n]);
    for ti=2:sz_in{1}(1)
        in=getInputs(ti);
        current=func(in{:});
        vals(ti,:)=reshape(current,[1,n]);
    end

    function in=getInputs(index)
        in=cell(n_in,1);
        for i2=1:n_in
            in{i2}=zeros(sz_in{i2}(2:end));
            in{i2}(:)=inputs{i2}(index,:);
        end
    end
end

