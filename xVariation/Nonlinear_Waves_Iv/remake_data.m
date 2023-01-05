function remake_data
    record = readtable('Results/wave_record.csv');
    in_table = strcmp(record.wave_type, "full");
    for i=19:size(in_table)
        if in_table(i) == 1
            fname = record.Name(i);
            if ~strcmp(fname,"master_wave_full.txt")
                fname
                theta = record.theta(i); 
                lambda = record.lambda(i);
                Fr = record.Fr(i);
                nu = record.nu(i);
                tau0 = record.tau0(i);
                d = record.d(i);
                alpha = record.alpha(i);
                params = [Fr,theta,lambda,nu,alpha,d,tau0];
                [xi_final,y_final] = bvp_full_from_master(params);
                out_vec = vertcat(xi_final,y_final);
                save("Results/"+fname,"out_vec","-ascii")
                write_record("Results/wave_record.csv",fname,{"full","water",Fr,theta,lambda,nu,alpha,d,tau0})
            end
        end
    end
end