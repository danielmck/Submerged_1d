function write_record(data_file,sim_type,N,h,d,reg_param,rho_r,alpha,phi_c,theta,eta_f,t_step)
    T = readtable('Results/result_record.csv');
    in_table = strcmp(T.Name, data_file);
    row_data = {data_file,sim_type,N,h,d,reg_param,rho_r,alpha,phi_c,theta,eta_f,t_step};
    if (sum(in_table) == 0)
        T = [T;row_data];
    else
        row_num = find(in_table,1);
        T(row_num,:) = row_data;
    end
    writetable(T,'Results/result_record.csv');
end