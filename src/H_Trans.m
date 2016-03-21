classdef H_Trans
    %H_Trans Provides helpful methods for creating and decomposing
    %Homogeneous Transformation matrices
    %   Detailed explanation goes here
    
    properties
        H = sym(eye(4));
    end
    
    properties (Dependent)
        Trans
        Rot
        Euler   %ZYX Euler Angles
    end
    
    methods
        function obj = H_Trans(varargin)
			if nargin>0
				obj.H=H_Trans.single(varargin{1});
				if nargin>1
					for i=2:length(varargin)
						obj.H=obj.H*H_Trans.single(varargin{i});
					end
				end
			end
        end
        
        function obj = set.Rot(obj,value)
            obj.H(1:3,1:3) = value(1:3,1:3);
        end
        function value = get.Rot(obj)
            value = obj.H(1:3,1:3);
        end
        
        function obj = set.Trans(obj,value)
            obj.H(1:3,4) = value;
        end
        function value = get.Trans(obj)
            value = obj.H(1:3,4);
        end
		
		function obj = set.Euler(obj,value)
            T=(H_Trans.rotZ(value(3))*H_Trans.rotY(value(2))*H_Trans.rotX(value(1)));
            obj.Rot=T.Rot;
        end
        function value = get.Euler(obj)
            R=obj.Rot;
            value=[atan2(R(3,2),R(3,3));
                   atan2(-R(1,3),sqrt(R(3,2)^2+R(3,3)^2));
                   atan2(R(2,1),R(1,1))];
        end
        
        function value = getRotVel(obj,var)
           w=simplify(diff(obj.Rot,var)*obj.Rot.');
           value=[w(3,2);w(1,3);w(2,1)];
        end
		
        function value = getJacobian(obj,q)
            value = sym(zeros(6,size(q,1)));
            for i=1:size(q,1)
                value(1:3,i) = simplify(diff(obj.Trans,q(i)));
                value(4:6,i) = simplify(obj.getRotVel(q(i)));
            end 
        end
        
        function func = getFunction(obj,input,expr)
            if nargin<3
                expr=obj.H;
            end
            
            func=H_Trans.createFunction(expr,input);
        end
        
        function obj = inv(obj)
			R=obj.Rot.';
			obj.H = [ R  -R*obj.Trans
					 0 0 0       1     ];
		end
    end
    
    methods
        function  c = mtimes(a,b)
            c=H_Trans(mtimes(a.H,b.H));
        end
        
        function  c = mrdivide(a,b)
            a=a.inv();
            c=H_Trans(mtimes(b.H,a.H));
        end
        
        function  c = mldivide(a,b)
            a=a.inv();
            c=H_Trans(mtimes(a.H,b.H));
        end
    end
    
    methods(Static)
        function obj = fromDH(DH)
            obj = H_Trans();
            for i = 1:size(DH,1)
                obj.H=obj.H*H_Trans.fromDH_single(DH(i,1),DH(i,2),DH(i,3),DH(i,4));
            end
            obj.H=simplify(vpa(obj.H));
        end
        
        function obj = rotX(theta)
            obj = H_Trans([1,0,0,0;
                    0,cos(theta),-sin(theta),0;
                    0,sin(theta),cos(theta),0;
                    0,0,0,1]);
        end
        function obj = rotY(theta)
            obj = H_Trans([cos(theta),0,sin(theta),0;
                    0,1,0,0;
                    -sin(theta),0,cos(theta),0;
                    0,0,0,1]);
        end
        function obj = rotZ(theta)
            obj = H_Trans([cos(theta),-sin(theta),0,0;
                    sin(theta),cos(theta),0,0;
                    0,0,1,0;
                    0,0,0,1]);
        end
        
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
                func = @(q)simplify(vpa(subs(expr,vars,extractVars(q))));
            else
                if iscell(input)
                    func = matlabFunction(expr,'Vars',input);
                else
                    func = matlabFunction(expr,'Vars',{input});
                end
            end
        end
    end
    
    methods (Static,Access=private)
        
        function [output_args] = fromDH_single( theta, d, a, alpha )
        theta=sym(theta);d=sym(d);a=sym(a);alpha=sym(alpha);
        output_args=[cos(theta),-sin(theta)*cos(alpha),sin(theta)*sin(alpha),a*cos(theta);...
            sin(theta),cos(theta)*cos(alpha),-cos(theta)*sin(alpha),a*sin(theta);...
            0,sin(alpha),cos(alpha),d;...
            0,0,0,1];
        end
        
        function M = single( input_args )
        M=sym(eye(4));
            if isequal(size(input_args),[4,4]),
                M = input_args;
			elseif isequal(size(input_args),[3,3]),
                M(1:3,1:3) = input_args(1:3,1:3);
			elseif isequal(size(input_args),[1,3]),
                M(1:3,4) = input_args(1,:);
            elseif isequal(size(input_args),[3,1]),
                M(1:3,4) = input_args(:,1).';
            end
        end
    end
end

