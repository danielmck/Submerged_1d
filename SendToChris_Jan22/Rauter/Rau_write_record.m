% Writes to csv file to store the parameters used in the simulation. They
% are:
% data_file - File Name
% N - number of points in z
% h - Dimensional flow height
% d - Dimensional particle diameter
% reg_param - mu regularisation parameter
% rho_r - Density ratio
% phi_c - Critical volume fraction
% theta - Slope angle
% eta_f - Nondimensional fluid viscosity
% a - Nondimensional compressubility
% phi_rcp - The random close packing fraction
% phi_rlp - The random loose packing fraction
% S0 - The shear limit

% NOTE: I have no idea why I did this but the first two values (h and d)
% are dimensional while the rest of the parameters aren't. This should
% probably be changed at some point to get rid of h and replace d with the
% dimensionless version.
function Rau_write_record(data_file,N,h,d,reg_param,rho_r,phi_c,theta,eta_f,a,phi_rcp,phi_rlp,S0)
    T = readtable('Results/result_record.csv');
    in_table = strcmp(T.Name, data_file);
    row_data = {data_file,N,h,d,reg_param,rho_r,phi_c,theta,eta_f,a,phi_rcp,phi_rlp,S0};
    if (sum(in_table) == 0)
        T = [T;row_data];
    else
        row_num = find(in_table,1);
        T = [T(1:row_num-1,:);row_data;T(row_num+1:end,:)];
    end
    writetable(T,'Results/result_record.csv');
end