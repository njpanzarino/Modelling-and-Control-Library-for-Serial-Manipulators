function rot = rotFuncs
    rot.rotx=@rotx;
    rot.roty=@roty;
    rot.rotz=@rotz;
end

function [ output_args ] = rotx( input_args )
%returns 3x3 rotation matrix by angle about x axis
output_args=[1,0,0;
    0,cos(input_args),-sin(input_args);
    0,sin(input_args),cos(input_args)];

end

function [ output_args ] = roty( input_args )
%returns 3x3 rotation matrix by angle about y axis
output_args=[cos(input_args),0,sin(input_args);
    0,1,0;
    -sin(input_args),0,cos(input_args)];

end

function [ output_args ] = rotz( input_args )
%returns 3x3 rotation matrix by angle about z axis
output_args=[cos(input_args),-sin(input_args),0;
    sin(input_args),cos(input_args),0;
    0,0,1];
end