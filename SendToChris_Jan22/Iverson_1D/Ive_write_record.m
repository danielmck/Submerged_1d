% Writes to csv file to store the parameters used in the simulation. They
% are:
% data_file - File Name
% sim_type - Type of simulation, the types are:
% dilatancy - "dil"
% diffusion only - "pdriv"
% constant pressure profile - "pcon"
% constant velocity profile - "ucon"
% N - number of points in z
% h - Dimensional flow height
% d - Dimensional particle diameter
% reg_param - mu regularisation parameter
% rho_r - Density ratio
% alpha - Nondimensional compressubility
% phi_c - Critical volume fraction
% theta - Slope angle
% eta_f - Nondimensional fluid viscosity

% NOTE: I have no idea why I did this but the first two values (h and d)
% are dimensional while the rest of the parameters aren't. This should
% probably be changed at some point to get rid of h and replace d with the
% dimensionless version.

function Ive_write_record(data_file,sim_type,N,h,d,reg_param,rho_r,alpha,phi_c,theta,eta_f)
    T = readtable('Results/result_record.csv');
    in_table = strcmp(T.Name, data_file);
    row_data = {data_file,sim_type,N,h,d,reg_param,rho_r,alpha,phi_c,theta,eta_f};
    if (sum(in_table) == 0)
        T = [T;row_data];
    else
        row_num = find(in_table,1);
        T = [T(1:row_num-1,:);row_data;T(row_num+1:end,:)];
    end
    writetable(T,'Results/result_record.csv');
end