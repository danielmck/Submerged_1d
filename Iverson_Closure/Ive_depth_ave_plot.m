function Ive_depth_ave_plot
    n_times = 5001;
    N = 200;
    h_vals = zeros(n_times,70);
    depth_ave_vals = zeros(n_times,70);
    base_shear_vals = zeros(n_times,70);
    cd Results
    for l=1:70        
        q=fix(l/10);
        r=mod(l,10);
        fname = "Ive_comp_"+num2str(q);
        if (r==0)
            fname = append(fname,"_deep.txt");
        else
            fname = append(fname,"_"+num2str(r)+"_deep.txt");
        end
        data_file = load(fname);
%         u_p = data_file(:,3*N+1:end);
%         dz = 1/(N-0.5);
%         dupdz = diff(u_p,1,2)/dz;
%         dupdz = horzcat(dupdz,zeros(n_times,1));
%         base_shear_vals(:,l+1) = dupdz(:,1);
%         is_flow = ((u_p > 5e-4) | (dupdz > 2e-3));
%         h_vals(:,l) = sum(is_flow,2)/N;
        depth_ave_vals(:,l) = data_file(:,1);
    end
%     save('Ive_flow_height_angles.txt', 'h_vals','-ascii')
%     save('Ive_shear_base_angles.txt', 'base_shear_vals','-ascii')
    save('Ive_basal_pe_angles.txt', 'depth_ave_vals','-ascii')
    cd ..
end