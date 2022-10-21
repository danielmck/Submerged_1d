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
    Fr_eq = 0.8;
    theta = 12;
    lambda = 13;
    nu = 1.13e-3;
    alpha = 1e-5;
    d = 1e-4; 
    var_list = [8,13,18];
    n_vals = size(var_list,2);
    xi_out = zeros(n_vals,n);
    h_out = zeros(n_vals,n);
    u_out = zeros(n_vals,n);
    xi_out = zeros(n_vals,n);
    phi_out = zeros(n_vals,n);
    pb_out = zeros(n_vals,n);
    pe_out = zeros(n_vals,n);
    C = viridis(n_vals+1);
    
    % Loads the master wave as the starter
    master_file = load("master_wave_full.txt");
    master_xi = master_file(1,:);
    master_y = master_file(2:end,:);
    master_cond = readtable("master_wave_full_cond.csv");
    master_param = [master_cond.Fr_eq,master_cond.theta,master_cond.lambda,master_cond.nu,master_cond.alpha,master_cond.d];
    
    % Loops over parameter values
    for i = 1:n_vals
        lambda = var_list(i); % Set the parameter you want to vary here
        % Solves for the waveform
        [xi_temp,y_out] = bvp_full_from_master([Fr_eq,theta,lambda,nu,alpha,d],master_y,master_xi,master_param);
        xi_out(i,:) = xi_temp;
        %         [xi_out(i,:),y_out] = bvp_non_pe_to_full(true,Fr_eq,nu,theta,lambda,alpha,d);
        master_param(1) = Fr_eq;
        h_out(i,:) = y_out(3,:);
        u_out(i,:) = y_out(1,1) - y_out(2,:)./h_out(i,:);
        phi_out(i,:) = y_out(6,:)./y_out(2,:);
        pb_out(i,:) = y_out(7,:) + rho/rho_f*chi.*h_out(i,:);
        pe_out(i,:) = pb_out(i,:) - h_out(i,:);
    end
    % Makes all of the plots
    waveforms = {h_out,u_out,phi_out,pb_out,pe_out};
    names = ["h","u","phi","pb","pe"];
    ax_title = ["h","\bar{u}","\bar{\phi}","p_b","p_e"];
    for j = 1:5
        SetPaperSize(8,8)
        plot_vec = waveforms{j};
        hold on
        for k=1:n_vals
            plot(xi_out(k,:),plot_vec(k,:), "DisplayName","$\lambda="+num2str(var_list(k))+"$","color",C(k,:))
        end
        xlabel("$\xi$")
        ylabel("$"+ax_title(j)+"$")
        % $\lambda = "+num2str(lambda)
        title("$Fr = "+num2str(Fr_eq)+"$, $\theta = "+num2str(theta)+"$, $\alpha ="+num2str(lambda)+"$, $d = "+num2str(d)+"$")
        legend("Location","best")
        figname = "full_wave_"+names(j)+"_comp_lambda2.pdf";
        exp_graph(gcf,figname)
        movefile(figname,"../Figures/FullWaveforms")
        clf;
    end
    close all
end