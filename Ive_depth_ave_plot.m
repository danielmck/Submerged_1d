function Ive_depth_ave_plot
    n_times = 5001;
    h_vals = zeros(n_times,71);
    u_vals = zeros(n_times,71);
    pe_vals = zeros(n_times,71);
    phi_vals = zeros(n_times,71);
    for l=1:70
        cd DA_Results
        q=fix(l/10);
        r=mod(l,10);
        fname = "Ive_da_"+num2str(q);
        if (r==0)
            fname = append(fname,"_deep.txt");
        else
            fname = append(fname,"_"+num2str(r)+"_deep.txt");
        end
        data_file = load(fname);
        theta = 0.1*l;
        cd ..
        h_vals(:,l+1) = data_file(:,1);
        phi_vals(:,l+1) = data_file(:,2);
        u_vals(:,l+1) = data_file(:,3);
        pe_vals(:,l+1) = data_file(:,4)-cosd(theta)*h_vals(:,l+1);
%         depth_ave_vals(:,l+1) = depth_average(u_p',N,n_times);
    end
    save('DA_flow_height_angles.txt', 'h_vals','-ascii')
    save('DA_phi_ave_angles.txt', 'phi_vals','-ascii')
    save('DA_u_ave_angles.txt', 'u_vals','-ascii')
    save('DA_pe_base_angles.txt', 'pe_vals','-ascii')
end