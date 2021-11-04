function EOS_write_record(data_file,N,h,d,reg_param,rho_r,phi_c,theta,eta_f,a,phi_rcp,phi_rlp,t_step,creep_type,S0)
    T = readtable('EqnOfState_Results/result_record.csv');
    in_table = strcmp(T.Name, data_file);
    row_data = {data_file,N,h,d,reg_param,rho_r,phi_c,theta,eta_f,a,phi_rcp,phi_rlp,t_step,creep_type,S0};
    if (sum(in_table) == 0)
        T = [T;row_data];
    else
        row_num = find(in_table,1);
        T = [T(1:row_num-1,:);row_data;T(row_num+1:end,:)];
    end
    writetable(T,'EqnOfState_Results/result_record.csv');
end