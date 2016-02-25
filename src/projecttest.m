n = inputdlg('Enter the DOF (degrees of freedom) of the robotic manipulator you wish to model. [Up to 6 DOF]','Number of DOFs')
options.Interpreter = 'tex'
options.Resize = 'on'
dh_table = inputdlg('Enter the dh table for the robotic manipulator Use the following format: \theta, d, a, \alpha  \theta_1,d_1,a_1,\alpha_1; \theta_2,d_2,a_2,\alpha_2;                        ...','dh table',[str2num(n{1}) 50],{''},options)
dh_data = str2num(dh_table{:})