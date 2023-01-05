function plot_full_waveforms
% Plots waveforms for the variation of different parameters. Outputs pdfs of 
% waves in h, u, phi, pb and pe. Parameters that can be varied are Fr,
% theta, lambda, alpha and d.
    n=100;
    rho_f = 1000;
    rho_p = 2500;
    phi_c = 0.585;
    
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    chi = (rho_f+3*rho)/(4*rho);
    
    % Default parameter values
    n_var = 6;
    Fr_eq = 0.8;
    theta = 12;
    lambda = 12;
    nu = 1.13e-4;
    alpha = 1e-5;
    d = 1e-4;
    params = [Fr_eq,theta,lambda,nu,alpha,d,0];
    var_list = [[0.8,1,1.2];[10,12,14];[8,12,16];[5e-5,1e-4,5e-4];[1e-6,1e-5,1e-4];[1e-5,1e-4,1e-3]];

    % Loads the master wave as the starter
%     master_file = load("Results/master_wave_full.txt");
%     master_xi = master_file(1,:);
%     master_y = master_file(2:end,:);
%     master_cond = readtable("Results/wave_record.csv");
%     master_param = [master_cond.Fr_eq,master_cond.theta,master_cond.lambda,master_cond.nu,master_cond.alpha,master_cond.d];
    
    var_names = ["Fr","\theta","\lambda","\nu","\alpha","d"];
    var_file = ["Fr","theta","lambda","nu","alpha","d"];
    names = ["h","u","phi","pb","pe"];
    ax_title = ["h","\bar{u}","\bar{\phi}","p_b","p_e"];
    num_type = ["%.2g","%.2g","%.2g","%.2e","%.0e","%.0e"];
    
    % Loops over parameter values
    for k = 4:4
        param_in = params;
        n_vals = size(var_list(k,:)~=0,2);
        xi_out = cell(n_vals,1);
        h_out = cell(n_vals,1);
        u_out = cell(n_vals,1);
        phi_out = cell(n_vals,1);
        pb_out = cell(n_vals,1);
        pe_out = cell(n_vals,1);
        C = viridis(n_vals+1);
        for i = 1:n_vals
            param_in(k) = var_list(k,i); % Set the parameter you want to vary here
            % Solves for the waveform
            [xi_temp,y_out] = bvp_full_from_master(param_in);
            xi_out{i,1} = xi_temp;
            %         [xi_out(i,:),y_out] = bvp_non_pe_to_full(true,Fr_eq,nu,theta,lambda,alpha,d);,master_y,master_xi,master_param
            h_out{i,1} = y_out(3,:);
            u_out{i,1} = y_out(1,1) - y_out(2,:)./h_out{i,:};
            phi_out{i,1} = y_out(6,:)./y_out(2,:);
            pb_out{i,1} = y_out(7,:) + rho/rho_f*chi.*h_out{i,:};
            pe_out{i,1} = pb_out{i,:} - h_out{i,:};
        end
        % Makes all of the plots
        waveforms = {h_out,u_out,phi_out,pb_out,pe_out};
        title_str = "$";
        for m = 1:n_var
            if m~=k
                val_str = num2str(params(m),num_type(m));
                title_str=title_str+var_names(m)+"="+val_str;
                title_str=title_str+"$, $";
            end
            if m==n_var
                title_str=extractBefore(title_str,strlength(title_str)-2);
            end
        end
        clf;
        for j = 1:5
            SetPaperSize(8,8)
            plot_vec = waveforms{j};
            hold on
            for l=1:n_vals
                plot(xi_out{l,1},plot_vec{l,:}, "DisplayName","$"+var_names(k)+"="+num2str(var_list(k,l),num_type(k))+"$","color",C(l,:))
            end
            xlabel("$\xi$")
            ylabel("$"+ax_title(j)+"$")
            % $\lambda = "+num2str(lambda)
            [t,s] = title(title_str);
            t.FontSize=t.FontSize-2;
            legend("Location","best")
            figname = "full_wave_"+names(j)+"_comp_"+var_file(k)+".pdf";
            exp_graph(gcf,figname)
            movefile(figname,"../Figures/FullWaveforms")
            clf;
        end
    end
    close all
end