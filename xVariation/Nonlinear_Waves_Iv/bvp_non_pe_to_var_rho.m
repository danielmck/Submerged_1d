function [xi_out,y_out] = bvp_non_pe_to_var_rho(custom_init,reverse,params,provide_init,master_xi,master_y)
% Converts the wave with no excess pressure into a waveform that satisfies
% the whole system. Achieves this by varying the paramter delta from 0 to 1
% to change between systems.

% Can be run with no inputs and just with the parameter values set in the
% code

% Can be run with custom_init set to true and with params as a 6 entry
% list [Fr,theta,lambda,nu,alpha,d]. If it is run with provide_init set to
% false then viscous_Iv_bvp_from_master is called to make input.

% When run with provide_init, the input xi and y must be provided.

% If reverse is set to true, a solution to the full system is converted to
% a solution to the non-pe case by running in reverse. In this case the
% initial y and xi must be provided.

    phi_c=0.585; % Volume fraction

    g=9.81; % m/s^2
    
%     eta_f = 1.18e-5;
%     rho_f = 1;
    
    rho_f = 1000;
    eta_f = 0.0010016; % Pa s
    rho_p = 2500;
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
        if (size(param_cell,2) == 7)
            [Fr,theta,lambda,nu,alpha,d,tau0] = param_cell{:};
            rel_flux = 1;
            pres_h = 0;
        elseif (size(param_cell,2) == 9)
            [Fr,theta,lambda,nu,alpha,d,tau0,rel_flux,pres_h] = param_cell{:};
        end
        if ~provide_init
            if reverse
                [master_xi,master_y] = bvp_full_from_master(true,params(1:4));
            else
                [master_xi,master_y] = viscous_Iv_bvp_from_master(true,params(1:4));
            end
        end
    else
        
        if ~provide_init
            if reverse
                master_name="master_wave_var_rho.txt";
            else
%                 master_name = "master_wave_no_pe.txt";
                master_name = "lambda_10_no_pe.txt";
            end
            master_file = load(strcat("Results/",master_name));
            master_xi = master_file(1,:);
            master_y = master_file(2:end,:);
        end
        record = readtable('Results/wave_record.csv');
        in_table = strcmp(record.Name, master_name);
        wave_type = record.wave_type(in_table);
        theta = record.theta(in_table); 
        lambda = record.lambda(in_table);
        Fr = record.Fr(in_table);
        nu_ratio = record.nu(in_table);
        tau0 = record.tau0(in_table);
        rel_flux = master_y(5,end);
        if reverse
            alpha = record.alpha(in_table);
            d = record.d(in_table);
        else
            alpha = 1e-5;
            d = 1e-4;
        end
        pres_h=strcmp(wave_type{1}(max(end-6,1):end),"_pres_h");
    end
    
    rho = rho_p*phi_c+rho_f*(1-phi_c);
%     P = (rho-rho_f)/rho;
    
    delta_vals = logspace(-6,0,100);
    delta_vals = horzcat(0,delta_vals);
    if reverse
        delta_vals = wrev(delta_vals(2:end));
    end
    % Checks that the full case is unstable at the requested parameter
    % values
    if ~reverse
        [h0_init, crit_Iv_init] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0);
        phi_init = phi_c/(1+sqrt(crit_Iv_init));
        u_init = Fr*sqrt(g*cosd(theta)*h0_init);

        crit_pb_init = rho_f*g*cosd(theta)*h0_init;
        nu_init = 3/2*mu_Iv_fn(crit_Iv_init)*eta_f/crit_Iv_init/rho*nu_ratio;
        nu_dl = nu_init/(u_init*h0_init);

        z_scale_init = h0_init;
        v_scale_init = u_init;
        p_scale_init = crit_pb_init;
        t_scale_init = z_scale_init/v_scale_init;

        g_dl_init = g*t_scale_init/v_scale_init; 
        rho_f_dl_init = rho_f*v_scale_init^2/p_scale_init;
        rho_p_dl_init = rho_p*v_scale_init^2/p_scale_init;
        rho_init = rho_p_dl_init*phi_init + rho_f_dl_init*(1-phi_init);
        chi_init = (rho_f_dl_init+3*rho_init)/(4*rho_init);
            
        k = 2*pi/lambda;

        A_mat = make_A_mat(k,rho_p,rho_f,theta,eta_f,d,alpha,Fr,crit_Iv_init,nu_init);
        A_eig = eigs(A_mat);
        [~, idx] = sort(imag(A_eig),'descend');
        A_eig = A_eig(idx);

        if (imag(A_eig(1))<0)
            error("Wave is not unstable, try a different value")
        end
        y6_in = phi_init*master_y(3,:);
        y7_in = master_y(3,:) - rho_init*g_dl_init*cosd(theta)*chi_init.*master_y(3,:);
        y_in = vertcat(master_y,y6_in,y7_in);
%         y_in=master_y;
    else
        y_in = master_y;
    end
    [xi_out,y_out] = run_bvp_iter(delta_vals, master_xi,y_in);
    if reverse
        no_pe_param = [Fr,theta,lambda,nu,tau0];
        [xi_out,y_out] = viscous_Iv_bvp_from_master(true,no_pe_param,true,xi_out,y_out(1:5,:),no_pe_param);
    end
    if ~custom_init && ~reverse && ~provide_init
        out_vec = vertcat(xi_out,y_out);
        filename = "lambda_20_var_rho.txt";
%         filename = "master_wave_var_rho.txt";
        save("Results/"+filename,"out_vec","-ascii")
        if pres_h
            type = "var_rho_pres_h";
        else
            type = "var_rho";
        end
        write_record("Results/wave_record.csv",filename,{type,"water",Fr,theta,lambda,nu_ratio,alpha,d,tau0})
    end

    u_w_f = y_out(1,1);
    Q1_f = y_out(2,:);
    h_f = y_out(3,:);
    u_f = u_w_f - Q1_f./h_f;
    n_f = y_out(4,:);
    m_f = y_out(5,:);
    if ~reverse
        phi_f = y_out(6,:)./Q1_f;
%         rho_fin = rho_p_dl*phi_f + rho_f_dl*(1-phi_f);
%         chi_f = (rho_f+3*rho_fin)/(4*rho_fin);
%         pb_f = y_out(7,:) + rho_fin*g_dl*cosd(theta)*chi_f.*h_f;
    end
    if ~custom_init
%         plot(xi_out,u)
        hold on
        plot(xi_out,h_f)
        plot(master_xi,master_y(3,:))
    end
    
    % Recursively runs if the solution does not converge or if it converges
    % to the unpeturbed solution
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
            [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0,false,true,delta_in);
    
            u_eq = Fr*sqrt(g*cosd(theta)*h0);
            phi_eq = phi_c/(1+sqrt(crit_Iv));
            rho_eq = phi_eq*rho_p+(1-phi_eq)*rho_f;
            
            p_tot = rho*g*cosd(theta);
            crit_pb = rho_f*g*cosd(theta)*h0;

            nu = 3/4*mu_Iv_fn(crit_Iv)*eta_f/crit_Iv/rho_eq*nu_ratio;
            nu_dl = nu/(u_eq*h0);

            z_scale = h0;
            v_scale = u_eq;
            p_scale = crit_pb;
            t_scale = z_scale/v_scale;

            tau0_dl = tau0/p_scale;
            eta_f_dl = eta_f/(p_scale*t_scale);
            alpha_dl = alpha*p_scale;
            g_dl = g*t_scale/v_scale; 

            u_eq_dl = u_eq/v_scale;
        %     p_tot_grad_dl = p_tot/p_scale*z_scale;
            rho_f_dl = rho_f*v_scale^2/p_scale;
            rho_p_dl = rho_p*v_scale^2/p_scale;
            rho_eq = rho_p_dl*phi_c + rho_f_dl*(1-phi_c);
            d_dl = d/z_scale;

            beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
            chi = (rho_f+3*rho)/(4*rho);
            
            solInit1=bvpinit(xi_in,@bvp_guess);
            opts = bvpset('RelTol',tol,'NMax',800*counter);
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
                h_wave = solN1.y(3,:);
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
        
        function dydxi = viscous_syst(xi,y)
            u_w = y(1);
            Q1 = y(2);
            h = y(3);
            u = u_w - Q1/h;
            n = y(4);
            m = y(5);
            phi = y(6)./Q1;
            rho_dl = delta_in*(rho_p_dl*phi+rho_f_dl*(1-phi))+(1-delta_in)*rho_eq;
            p_tot_grad_dl = rho_dl*g_dl*cosd(theta);
            chi = (rho_f_dl+3*rho_dl)/(4*rho_dl);
            pb = y(7) + rho_dl*g_dl*cosd(theta)*chi.*h;
            P = (rho_dl-rho_f_dl)/rho_dl;

            zeta = 3/(2*alpha_dl*h) + P/4;
            p_p = p_tot_grad_dl*h-pb;
            D = -2/beta_dl/h*(pb-h);
            Iv = 3*eta_f_dl*abs(u)/h/p_p;
            R_w3 = -phi*rho_f_dl/rho_dl*D;
            R_w4 = (-P/4+zeta)*D - 9/2/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));

            dhdxi = n;
            n_coeff = h/Fr^2-Q1^2/h^2;
            
            mu_val = delta_in*p_p/rho_dl*mu_Iv_fn(Iv)+(1-delta_in)*P/Fr^2*mu_Iv_fn(crit_Iv*abs(u)/h^2)*h+tau0_dl*rho_f/rho/h; 
            force_bal = h*tand(theta)/Fr^2-sign(u).*mu_val-delta_in*(u_w-u)*P*D;
            n_eq = (force_bal)./n_coeff;
            dQdxi = -delta_in*P*D;
            
%             dndxi = 1/(2*h)*n^2 + h^(3/2)/Fr^2/nu_dl/Q1*n_coeff*(n-n_eq);
            dmdxi = rho_dl/rho_eq*h/lambda*u^(1-pres_h);

            dy6dxi = -delta_in*R_w3;
            dy7dxi = delta_in*R_w4/(u-u_w);
            dphidxi = (dy6dxi-phi*dQdxi)/Q1;
            dPdxi = rho_f_dl*(rho_p_dl-rho_f_dl)./rho_dl^2*dphidxi;
            dpbdxi = dy7dxi + rho_dl*g_dl*cosd(theta)*chi*dhdxi + 3/4*g_dl*cosd(theta)*(rho_p_dl-rho_f_dl)*dphidxi;
            dDdxi = -2/beta_dl*(dpbdxi/h-pb/h^2*dhdxi);
%             dndxi = 1/Q1/nu_dl*(n_coeff*n-force_bal+delta_in*u*dQdxi)-delta_in*h/Q1*(P*dDdxi+D*dPdxi);
%             dndxi_old = 1/Q1/nu_dl*(n_coeff*n - force_bal-nu_dl*(h*D*dPdxi-P*2/beta_dl*dpbdxi));
            dndxi = 1/h*n^2-delta_in*h*D/Q1*dPdxi-delta_in*h*P/Q1*dDdxi+delta_in*P*D/Q1*n+h/Q1/nu_dl*(n_coeff*n - force_bal); %
            dydxi = [0,dQdxi,dhdxi,dndxi,dmdxi,dy6dxi,dy7dxi]';  
        end
        
        function guess = bvp_guess(xi)
            guess_ind = sum(xi_in<xi)+1;
            if guess_ind == max(size(xi_in))
                guess = y_in(:,end);
            else
                gap = xi_in(guess_ind+1)-xi_in(guess_ind);
                guess = y_in(:,guess_ind)*(xi-xi_in(guess_ind))/gap + y_in(:,guess_ind+1)*(xi_in(guess_ind+1)-xi)/gap;
            end
        end
    end
    
    function resid = bc_vals(ya,yb)
       resid = [ya(3)-yb(3), ya(4), yb(4), ya(5), yb(5)-rel_flux, ya(6)-yb(6), ya(7)-yb(7)]; 
    end
end