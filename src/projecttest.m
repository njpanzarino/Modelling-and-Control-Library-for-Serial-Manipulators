n = inputdlg('Enter the DOF (degrees of freedom) of the robotic manipulator you wish to model. [Up to 6 DOF]','Number of DOFs')
options.Interpreter = 'tex'
options.Resize = 'on'
dh_table = inputdlg('Enter the dh table for the robotic manipulator Use the following format: \theta, d, a, \alpha  \theta_1,d_1,a_1,\alpha_1; \theta_2,d_2,a_2,\alpha_2;                        ...','dh table',[str2num(n{1}) 50],{''},options)
dh_data = cellstr(dh_table{1, 1})
n = str2num(cell2mat(n))
M = {}
i = 1
for i = 1:n
    M(i, 1:4) = strsplit(dh_data{i},',')
end
j = 1
k = 1
[h, w] = size(M)
DH = sym('DH',[n 4])
for j = 1:h
    for k = 1:w
       DH(j,k) = sym(M{j,k})
    end
end

Kin_Model.DH