function [xi_final,y_final] = no_pe_full_depo_wave(custom_init,reverse,params,provide_init,master_xi,master_y)
% Does not work, an attempt to convert the no pe wave with no viscosity
% into the full wave
% Likely doesn't work as there is no way to find the flow height gradient
% at the critical point.
    [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
    P = (rho-rho_f)/rho;
    if ~exist("custom_init","var")
        custom_init = false;
    end
    if ~exist("reverse","var")
        reverse = false;
    end
    if ~exist("provide_init","var")
        provide_init = false;
    end
    % If there are custom initial parameters, the function may have to call
    % viscous_Iv_bvp_from_master to make the input
    if custom_init
        param_cell = num2cell(params);
        [Fr,theta,h_alt,alpha,d,tau0] = param_cell{:};
        if ~provide_init
            if reverse
                [master_xi,master_y] = bvp_full_from_master(true,params(1:4));
            else
                [master_xi,master_y] = viscous_Iv_bvp_from_master(true,params(1:4));
            end
        end
    else
        master_name = "no_pe_no_vis_convert.txt";
        master_file = load("Results/"+master_name);
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        record = readtable('Results/wave_record.csv');

        in_table = strcmp(record.Name, master_name);
        wave_type = record.wave_type(in_table);
        theta = record.theta(in_table);
        h_alt = record.h_alt(in_table);
        Fr = record.Fr(in_table);
        tau0 = record.tau0(in_table);
        stat_len = 0;
        if reverse
            alpha = record.alpha(in_table);
            d = record.d(in_table);
        else
            alpha = 1e-5;
            d = 1e-4;
        end
    end
    
    
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
    
    [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0);
    u_eq = Fr*sqrt(g*cosd(theta)*h0);
    phi_eq = phi_c/(1+sqrt(crit_Iv));
    p_tot = rho*g*cosd(theta);
    crit_pb = rho_f*g*cosd(theta)*h0;
    
    z_scale = h0;
    v_scale = u_eq;
    p_scale = crit_pb;
    t_scale = z_scale/v_scale;

    tau0_dl = tau0/p_scale;
    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    g_dl = g*t_scale/v_scale; 
    
    u_eq_dl = u_eq/v_scale;
    p_tot_grad_dl = p_tot/p_scale*z_scale;
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale; 
    rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
    d_dl = d/z_scale;
    
    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
    chi = (rho_f+3*rho)/(4*rho);
    
    delta_vals = logspace(-6,0,50);
    delta_vals = horzcat(0,delta_vals);
    if reverse
        delta_vals = wrev(delta_vals(2:end));
    end
%     stable = viscous_stability(theta,Fr,nu,lambda);

    if ~reverse
        y6_in = phi_eq*master_y(2,:);
        y7_in = master_y(4,:) - rho_dl*g_dl*cosd(theta)*chi.*master_y(4,:);
        y_in = vertcat(master_y,y6_in,y7_in);
%         y_in=master_y;
    else
        y_in = master_y;
    end
    [xi_out,y_out] = run_bvp_iter(delta_vals, master_xi,y_in);
    if reverse
        no_pe_param = [Fr,theta,h_alt,tau0];
        [xi_out,y_out] = viscous_Iv_bvp_from_master(true,no_pe_param,true,xi_out,y_out(1:5,:),no_pe_param);
    end
    if ~custom_init && ~reverse && ~provide_init
        out_vec = vertcat(xi_out,y_out);
        filename = "no_vis_full_master.txt";
        save("Results/"+filename,"out_vec","-ascii")
        write_record("Results/wave_record.csv",filename,{"full_no_vis","water",Fr,theta,h_alt,alpha,d,tau0})
    end
    
    function [xi_in,y_in] = run_bvp_iter(delta_vals, xi_in, y_in, tol,counter)
        if ~exist('counter','var')
            counter = 1;
        end
        if ~exist('tol','var')
            tol = 1e-4;
        end
        if counter > 10
            error("Max iteration depth reached, non convergence")
        end
        n_step = size(delta_vals,2);
        for i = 2:n_step
            delta_in = delta_vals(i);
            solInit1=bvpinit(xi_in,@bvp_guess);
            opts = bvpset('RelTol',tol,'NMax',200*counter);
            try
                solN1 = bvp4c(@viscous_syst,@bc_vals,solInit1,opts);
                resid = solN1.stats.maxres;
            catch ME
                switch ME.identifier
                    case 'MATLAB:bvp4c:SingJac'
                        warning('Singular Jacobian encountered, reducing step size');
                        resid = tol+1;
                    otherwise
                        rethrow(ME)
                end    
            end
            if resid < tol
                h_wave = solN1.y(4,:);
                h_diff = max(h_wave)-min(h_wave);
                if (h_diff>1e-4)
                    y_in = solN1.y;
                    xi_in = solN1.x;
                else
                    [xi_in,y_in] = run_bvp_iter(linspace(delta_vals(i-1),delta_in,3)...
                    ,xi_in, y_in, tol,counter+1);
                end
            else
                [xi_in,y_in] = run_bvp_iter(linspace(delta_vals(i-1),delta_in,3)...
                    ,xi_in, y_in, tol,counter+1);
            end
        end
        
        function guess = bvp_guess(xi)
            % Initial guess function from the ode soln
            guess_ind = sum(xi_in<xi)+1;
            if guess_ind == max(size(xi_in))
                guess = y_in(:,end);
            else
                gap = xi_in(guess_ind+1)-xi_in(guess_ind);
                guess = y_in(:,guess_ind)*(xi-xi_in(guess_ind))/gap + y_in(:,guess_ind+1)*(xi_in(guess_ind+1)-xi)/gap;
            end  
        end
        
        function dydxi = viscous_syst(xi,y)
            % The system from the viroulet paper with a different friction
            % law. Need to solve for u_w and Q1 too.
            u_w = y(1);
            Q1 = y(2);
            lambda = y(3);
            h = y(4);
            u = (-Q1 + h.*u_w)./h;
            m = y(5);
            phi = y(6)./Q1;
            pb = y(7) + rho_dl*g_dl*cosd(theta)*chi.*h;
            
            zeta = 3/(2*alpha_dl*h) + P/4;
            p_p = p_tot_grad_dl*h-pb;
            Iv = 3*eta_f_dl*abs(u)/h/p_p;
            D = -2/beta_dl/h*(pb-h);
            
            denom = (h.^3/Fr^2-Q1.^2);
            fric_coeff = delta_in*p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv)+(1-delta_in)*P*mu_Iv_fn(crit_Iv*abs(u)/h^2);
            fb_val = tand(theta)-fric_coeff-tau0_dl*rho_f/rho/h;
            dhdxi = 1/Fr^2.*h.^3.*fb_val./denom;
            
            R_w3 = -phi*rho_f_dl/rho_dl*D;
            R_w4 = (-P*chi+zeta)*D - 2*3/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));
            
            dQdxi = -delta_in*P*D;
            dmdxi = h/(lambda+stat_len)*u;
            dy6dxi = -delta_in*R_w3;
            dy7dxi = delta_in*R_w4/(u-u_w);
            dydxi = [0,dQdxi,0,dhdxi,dmdxi,dy6dxi,dy7dxi]'*lambda;
        end
        
        function resid = bc_vals(ya,yb)
            h_crit = (ya(2)*Fr)^(2/3);
            h_start = ya(2)/ya(1)+h_alt;
            h_stop = (-h_start + sqrt(h_start.^2+8*ya(2)^2./(h_start/Fr^2)))/2;
            resid = [ya(2)-yb(2), ya(4)-h_start, yb(4)-h_stop, ya(5), yb(5)-1, ya(6)-yb(6), ya(7)-yb(7)]; 
        end
    end
end