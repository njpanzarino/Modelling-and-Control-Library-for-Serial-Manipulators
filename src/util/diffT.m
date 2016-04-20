function [equ,q,d_q,dd_q]=diffT(equ,q,d_q)
%DIFFT Used to differentiate a symbolic expression with respect to time, 
% even if the variables are not expressed with respect to time.
%   equ: the expression to differentiate
%   q: vector of variables which vary over time 
%   d_q: (OPTIONAL) Time derivative of variables in q. specify only if
%       known or performing the 2nd derivative.
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
