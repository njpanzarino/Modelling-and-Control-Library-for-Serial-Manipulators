function [equ,q,d_q,dd_q]=diffT(equ,q,d_q)
    sz=size(q);
    
    t=sym('t',real);

    if nargin<3 || isempty(d_q)
        d_q=sym('d_q',sz);
        for i=1:sz
            d_q(i)=sym(strcat('d_',char(q(i))),'real');
        end
    end
    dd_q=sym('dd_q',sz);
    
    q_t=sym('q_t',sz);
    d_q_t=sym('d_q_t',sz);
    dd_q_t=sym('dd_q_t',sz);

    for i=1:sz
        dd_q(i)=sym(strcat('dd_',char(q(i))),'real');
        
        q_t(i) = sym(strcat(char(q(i)),'(',char(t),')'),'real');
        d_q_t(i) = diff(q_t(i),t);
        dd_q_t(i) = diff(q_t(i),t,2);
    end

    equ=subs(equ,q,q_t);
    equ=subs(equ,d_q,d_q_t);
    equ=subs(equ,dd_q,dd_q_t);

    equ=diff(equ,obj.t);

    equ=subs(equ,dd_q_t,dd_q);
    equ=subs(equ,d_q_t,d_q);
    equ=subs(equ,q_t,q);
end
