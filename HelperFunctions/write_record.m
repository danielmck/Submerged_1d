function write_record(db_file,data_file,params)
    T = readtable(db_file);
    in_table = strcmp(T.Name, data_file);
    row_data = [{data_file},params];
    if (sum(in_table) == 0)
        T = [T;row_data];
    else
        row_num = find(in_table,1);
        T = [T(1:row_num-1,:);row_data;T(row_num+1:end,:)];
    end
    writetable(T,db_file);
end