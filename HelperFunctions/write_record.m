function write_record(db_file,data_file,params, text_params)
    T = readtable(db_file);
    in_table = strcmp(T.Name, data_file);
    if exist('text_params','var')
        row_data = array2table([data_file,text_params,params]);
    else
        row_data = array2table([data_file,params]);
    end
    row_data.Properties.VariableNames = T.Properties.VariableNames;
    if (sum(in_table) == 0)
        T = [T;row_data];
    else
        row_num = find(in_table,1);
        T = [T(1:row_num-1,:);row_data;T(row_num+1:end,:)];
    end
    writetable(T,db_file);
end