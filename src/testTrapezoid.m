function testTrapezoid(initial,final,v_max,a_max)
    
    if nargin < 4 || isempty(a_max)
        a_max=1;
    end
    if nargin < 3 || isempty(v_max)
        v_max=2;
    end
    if nargin < 2 || isempty(final)
        final=[5,0;2,0];
    end
    if nargin < 1 || isempty(initial)
        initial=[0,0;-1,0];
    end

    [plan,tf]=Planner.trapezoid(initial,final,v_max,a_max);

    range=0:0.01:tf;

    pos=zeros(2,numel(range));
    vel=zeros(2,numel(range));
    acc=zeros(2,numel(range));

    for i=1:numel(range);
        d=plan(range(i));
        pos(:,i)=d(:,1);
        vel(:,i)=d(:,2);
        acc(:,i)=d(:,3);
    end

    figure
    subplot(1,3,1);
    plot(range,pos);
    title('position');
    subplot(1,3,2);
    plot(range,vel);
    title('velocity');
    subplot(1,3,3);
    plot(range,acc);
    title('acceleration');
end


