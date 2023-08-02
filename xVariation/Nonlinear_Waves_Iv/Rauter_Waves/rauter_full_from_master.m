function [xi_final,y_final] = rauter_full_from_master(params,master_y,master_xi,master_params)
    % Uses the master solution to paramfrom_masterfrom_mastereter step to the desired values
    % specified in params. If the params vector is not set then the default
    % values specified are used.
    
    % The parameters that can be set are Fr, theta, lambda, alpha and d.
    
    % If desired, a different master wave can be specified in order to
    % speed up computation. The master y, xi and parameter values must be
    % specified.
    
    Res_dir = "~/Documents/MATLAB/1D_System/xVariation/Nonlinear_Waves_Iv/Rauter_Waves/Results/";

    phi_c=0.585; % Volume fraction

    g=9.81; % m/s^2
    
%     eta_f = 1.18e-5;
%     rho_f = 1;
    
    rho_f = 1000;
    eta_f = 0.0010016; % Pa s
    rho_p = 2500;
    phi_rcp = 0.63;
    phi_rlp = 0.53;
    
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    chi = (rho_f+3*rho)/(4*rho);
    P = (rho-rho_f)/rho;
    
    % Defines parameters if not specified
    if ~exist('params','var')
        custom_init = false;
        wave_type = "full_pres_h";
%         Fr_eq = 0.8;
        
        theta = 12;
        lambda = 12;
        nu = 2e-4;
        a = 130;
        d = 1e-4;
        tau0 = 0; % Pa
        rel_flux = 1;
        pres_h = (wave_type == "full_pres_h");
        
        h0 = 0.0075;
        [Fr_eq, ~] = crit_Iv_rauter_h(theta, rho_p, rho_f, eta_f, h0, a, phi_rlp, phi_rcp, tau0);
        params = [Fr_eq,theta,lambda,nu,a,d,tau0,rel_flux,pres_h];
        filename = "Comp_to_Ive.txt";
    else
        custom_init = true;
        param_cell = num2cell(params);
        if (size(param_cell,2) == 7)
            [Fr_eq,theta,lambda,nu,a,d,tau0] = param_cell{:};
            rel_flux = 1;
            pres_h = 0;
        elseif (size(param_cell,2) == 9)
            [Fr_eq,theta,lambda,nu,a,d,tau0,rel_flux,pres_h] = param_cell{:};
        end
    end
    
    if ~exist('master_y','var')
        master_name = "master_pres_h.txt";
        master_file = load(Res_dir+master_name);
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        record = readtable(Res_dir+'full_record.csv');

        in_table = strcmp(record.Name, master_name);
        master_wave_type = record.wave_type(in_table);
        master_pres_h = (master_wave_type == "full_pres_h");
        master_theta = record.theta(in_table); 
        master_lambda = record.lambda(in_table);
        master_Fr = record.Fr(in_table);
        master_nu = record.nu(in_table);
        master_tau0 = record.tau0(in_table);
%         if strcmp(wave_type,'full')
        master_d = record.d(in_table);
        master_a = record.a(in_table);
        master_u_w = record.u_w(in_table);
        master_rf = 1;
%         end
        master_params = [master_Fr,master_theta,master_lambda,master_nu,master_a,master_d,master_tau0];
    else
        
        master_temp = num2cell(master_params);
        if (size(param_cell,2) == 7)
            [master_Fr,master_theta,master_lambda,master_nu,master_a,master_d,master_tau0] = master_temp{:};
            rel_flux = 1;
            pres_h = 0;
        elseif (size(param_cell,2) == 9)
            [master_Fr,master_theta,master_lambda,master_nu,master_a,master_d,master_tau0,master_rf,master_pres_h] = master_temp{:};
        end
    end
    
    
    
    if (pres_h ~= master_pres_h)
        m_init = zeros(size(master_xi));
        m_val = 0;
        u_flux = (master_y(1,1) - master_y(2,:)./master_y(3,:)).^(1-pres_h);
        for j = 1:size(master_xi,2)
            m_init(j) = m_val;
            if j < size(master_xi,2)
                m_val = m_val + 1/master_lambda* (master_y(3,j)*u_flux(j)+master_y(3,j+1)*u_flux(j+1))/2*(master_xi(j+1)-master_xi(j));
            end
        end
        master_y = vertcat(master_u_w*ones(size(master_xi)),master_y(1:3,:),m_init,master_y(5,:));
        master_rf = m_val;
    end

    Fr_ratio = max(Fr_eq/master_Fr,master_Fr/Fr_eq);
    nu_ratio = max(nu/master_nu,master_nu/nu);
    theta_ratio = max(theta/master_theta,master_theta/theta);
    lambda_ratio = max(lambda/master_lambda,master_lambda/lambda);
    a_ratio = 1+max(log(a/master_a),log(master_a/a));
    d_ratio = 1+max(log(d/master_d),log(master_d/d));
    tau0_ratio = abs(tau0-master_tau0)/(max(master_tau0,tau0)+1e-6);
    rf_ratio = abs(rel_flux - master_rf);
    
    ratio_sum = Fr_ratio-1+nu_ratio-1+50*(theta_ratio-1)+d_ratio-1+a_ratio-1+tau0_ratio+35*(lambda_ratio-1)+rf_ratio;
    if ratio_sum > 1e-6 %
        n_steps = max(min(ceil(30*(ratio_sum-1)),400),90);
    else
        n_steps=2;
    end
    
%     if lambda_ratio>1
%         [xi_nope,y_nope] = bvp_non_pe_to_full(true,true,master_params,true,master_xi,master_y);
%         y_nope = y_nope(1:5,:);
%         new_params = [master_params(1:2),lambda,master_params(4),master_params(7)];
%         [xi_nope,y_nope] = viscous_Iv_bvp_from_master(true,new_params,true,xi_nope,y_nope,master_params);
%         [xi_lambda,y_lambda] = bvp_non_pe_to_full(true,false,[new_params(1:4),master_params(5:7)],true,xi_nope,y_nope);
%     else
        xi_lambda = master_xi;
        y_lambda = master_y;
%     end
    % Sets lists for parameter stepping between master and designated
    % values
    Fr_list = linspace(master_Fr,Fr_eq,n_steps);
    nu_list = linspace(master_nu,nu,n_steps);
    theta_list = linspace(master_theta,theta,n_steps);
    lambda_list = linspace(master_lambda,lambda,n_steps);
%     lambda_list = lambda*ones(1,n_steps);
    a_list = logspace(log10(master_a),log10(a),n_steps);
    d_list = logspace(log10(master_d),log10(d),n_steps);
    tau0_list = linspace(master_tau0,tau0,n_steps);
    rf_list = linspace(master_rf, rel_flux, n_steps);

    [xi_final,y_final] = run_bvp_step(Fr_list, theta_list, lambda_list, nu_list, a_list, d_list, tau0_list, rf_list, xi_lambda, y_lambda);
%     [xi_nope,y_nope] = viscous_Iv_bvp_from_master(true,[Fr_eq,theta,lambda,nu]);
    
    if ~custom_init
        u_w = y_final(1,1);
        Q1 = y_final(2,:);
        h = y_final(3,:);
        u = u_w - Q1./h;
        n = y_final(4,:);
        m = y_final(5,:);
        phi = y_final(6,:)./Q1;
        out_vec = vertcat(xi_final,y_final(2:6,:));
        
        write_record("Results/full_record.csv",filename,{"full","water",Fr_eq,theta,nu,a,phi_rlp,phi_rcp,d,tau0,y_final(1,1),lambda,0})
        save(Res_dir+filename,"out_vec","-ascii")
%         plot(xi_out,u)

        hold on
        plot(xi_final,u)
%         plot(master_xi,master_y(3,:))
    end
    
    hold on
%     plot(xi_nope,y_nope(3,:))
    
    function [xi_out, y_out] = run_bvp_step(Fr_vals, theta_vals, lambda_vals, nu_vals, a_vals, d_vals, tau0_vals, rf_vals, xi_in, y_in, tol,counter)
        % If not specified, the depth counter is set to 1.
        if ~exist('counter','var')
            counter = 1;
        end
        % Tolerance default
        if ~exist('tol','var')
            tol = 1e-3;
        end
        % If depth reaches 5, attempts to convert to non pe case and step
        % that way.
        if counter > 5
            out_vec = vertcat(xi_in,y_in(2:6,:));
            fail_name = filename;
            write_record("Results/full_record.csv",filename,{"full","water",Fr_eq,theta,nu,a,phi_rlp,phi_rcp,d,tau0,y_in(1,1),lambda,0})
            save(Res_dir+fail_name,"out_vec","-ascii")
            error("Max iteration depth reached, non convergence. Trying no excess pressure method")
            y_in = -1;
            xi_in = -1;
        end             
        n_step = size(Fr_vals,2);
        % Iterates over the number of steps
        for i = 2:n_step
            % Only carries out if not at max depth
            if ~isequal(y_in,-1)
                theta_in = theta_vals(i);
                lambda_in = lambda_vals(i);
                Fr_in = Fr_vals(i);
                tau0_in = tau0_vals(i);
                rf_in = rf_vals(i);
                a_in = a_vals(i);
                xi_run = xi_in/lambda_vals(i-1)*lambda_vals(i);
%                 xi_run = horzcat(linspace(0,3*lambda_vals(i)/4,70),linspace(3*lambda_vals(i)/4, lambda_vals(i),70));
                [h0, crit_Iv] = crit_Iv_rauter(theta_in, rho_p, rho_f, eta_f, Fr_in, a_in, phi_rlp, phi_rcp, tau0_in);
                u_eq = Fr_in*sqrt(g*cosd(theta_in)*h0);
                pp0 = 3*eta_f*u_eq/h0/crit_Iv; 
                rho_eq = pp0/h0/g/cosd(theta_in)+rho_f;
                
                phi_eq = (rho_eq-rho_f)/(rho_p-rho_f);
                crit_pb = rho_f*g*cosd(theta_in)*h0;
                nu_dl = nu_vals(i)/(u_eq*h0);
                
%                 h_stop = tau_0/rho*g*cosd(theta)/(rho/(rho-rho_f)-mu1_Iv);
                
                z_scale = h0;
                v_scale = u_eq;
                p_scale = crit_pb;
                t_scale = z_scale/v_scale;

                eta_f_dl = eta_f/(p_scale*t_scale);
                a_dl = a_in/p_scale;
                g_dl = g*t_scale/v_scale; 

                tau0_dl = tau0_in/p_scale;
                rho_f_dl = rho_f*v_scale^2/p_scale;
                rho_p_dl = rho_p*v_scale^2/p_scale;
                rho_eq_dl = rho_eq*v_scale^2/p_scale; 
                d_dl = d_vals(i)/z_scale;

                beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);

                % Solves the stepped system
                solInit1=bvpinit(xi_run,@bvp_guess);
                opts = bvpset('RelTol',tol,'NMax',500*counter);
                try
                    solN1 = bvp4c(@rauter_syst,@bc_vals,solInit1,opts);
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
                    % Ensures the residual is below the tolerance
                if resid < tol
                    h_wave = solN1.y(3,:);
                    h_diff = max(h_wave)-min(h_wave);
                    % Ensures it is not the base state solution
                    if (h_diff>1e-4)
                        xi_in = horzcat(linspace(0,3*lambda_in/4,70),linspace(3*lambda_in/4+0.001, lambda_in,70));
                        y_in = interp1(solN1.x,solN1.y',xi_in)';
                    else
                        [xi_in,y_in,non_pe] = recurse_manage();
                        if non_pe
                            break
                        end
                    end
                   
                else
                    [xi_in,y_in,non_pe] = recurse_manage();
                    % (i,Fr_vals, nu_vals, theta_vals, lambda_vals, alpha_vals, d_vals, xi_in, y_in, tol,counter+1)
                    
                    if non_pe
                        break
                    end
                end
            end
        end
        xi_out = xi_in;
        y_out = y_in;
        
        
        function [xi_return,y_return,non_pe] = recurse_manage()
            % Runs the recursive calling of the stepping function
            [xi_return,y_return] = run_bvp_step(linspace(Fr_vals(i-1),Fr_in,3)...
                         ,linspace(theta_vals(i-1),theta_vals(i),3),linspace(lambda_vals(i-1),lambda_in,3),...
                         linspace(nu_vals(i-1),nu_vals(i),3),linspace(a_vals(i-1)...
                        ,a_vals(i),3),linspace(d_vals(i-1),d_vals(i),3),linspace(tau0_vals(i-1),tau0_in,3),...
                        linspace(rf_vals(i-1),rf_in,3),xi_in, y_in, tol,counter+1);
            non_pe = false;
            if counter == 1
                if isequal(y_return,-1)                
                    curr_params = [Fr_vals(i-1),theta_vals(i),lambda_vals(i-1),nu_vals(i-1),alpha_vals(i-1),d_vals(i-1),tau0_vals(i-1),rf_vals(i-1),pres_h];
                    no_pe_curr_params = [curr_params(1:4),curr_params(7:end)];
                    no_pe_params = [params(1:4),params(7:end)];
                    [xi_no_pe,y_no_pe] = bvp_non_pe_to_full(true,true,curr_params,true,xi_in,y_in);
                    y_no_pe = y_no_pe(1:5,:);
                    [xi_no_pe,y_no_pe] = viscous_Iv_bvp_from_master(true,no_pe_params,true,xi_no_pe,y_no_pe,no_pe_curr_params);
                    [xi_return,y_return] = bvp_non_pe_to_full(true,false,params,true,xi_no_pe,y_no_pe);
                    non_pe = true;
                end
            end
                    
        end
        
        function guess = bvp_guess(xi)
            guess_ind = sum(xi_run<xi)+1;
            if guess_ind == max(size(xi_run))
                guess = y_in(:,end);
            else
                gap = xi_run(guess_ind+1)-xi_run(guess_ind);
                guess = y_in(:,guess_ind)*(xi-xi_run(guess_ind))/gap + y_in(:,guess_ind+1)*(xi_in(guess_ind+1)-xi)/gap;
            end
        end
            
        function dydxi = rauter_syst(xi,y)
            u_w = y(1);
            Q1 = y(2);
            h = y(3);
            u = u_w - Q1/h;
            n = y(4);
            m = y(5);
            phi = y(6)./Q1;

            rho_dl = phi*rho_p_dl+(1-phi)*rho_f_dl;
            p_tot_grad_dl = rho_dl*g_dl*cosd(theta_in);
            p_p_contact = max(a_dl*(phi-phi_rlp)/(phi_rcp-phi),0);
            Iv_phi = (phi_c/phi-1)^2;
            p_p_shear = 3*u*eta_f_dl/Iv_phi/h;
            p_p = (p_p_contact+p_p_shear);
            pb = p_tot_grad_dl*h-p_p;
            D = -2/beta_dl/h*(pb-h);
            R_w3 = -phi*rho_f_dl/rho_dl*D;

            dhdxi = n;
            n_coeff = 1-Q1.^2.*Fr_eq^2/h^3;
            Iv = 3*eta_f_dl*abs(u)/h/p_p;
            mu_val = p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv)+tau0_dl*rho_f_dl/rho_dl/h;
            force_bal = h*(tand(theta_in)-sign(u).*mu_val)/Fr_eq^2+u*P*D;
            n_eq = (force_bal)./n_coeff;
            dQdxi = -P*D;
            
            dmdxi = h*rho_dl/rho_eq_dl/lambda*u^(1-pres_h);
            dy5dxi =-R_w3;
            dphidxi = (dy5dxi-phi*dQdxi)/Q1;
            dPdxi = (rho_p_dl-rho_f_dl)*rho_f_dl/rho_dl^2*dphidxi;
            dppdh = 3*eta_f_dl/Iv_phi/h*(Q1/h^2)-3*u*eta_f_dl/Iv_phi/h^2;
            dppdQ = 3*eta_f_dl/Iv_phi/h*(-1/h);
            dppdphi = a_dl*(phi_rcp-phi_rlp)/(phi_rcp-phi)^2+6*u*eta_f_dl/h/(phi_c/phi-1)^3*phi_c/phi^2;
            dppdxi = dppdh.*dhdxi+dppdQ.*dQdxi+dppdphi.*dphidxi;
            dDdxi = -2/beta_dl*((rho_p_dl-rho_f_dl)-dppdxi/h+p_p/h^2*dhdxi);
            dndxi = 1/Q1/nu_dl*(n_coeff*n-force_bal+u*dQdxi)-h/Q1*(P*dDdxi+D*dPdxi);
            dydxi = [0,dQdxi,dhdxi,dndxi,dmdxi,dy5dxi];
        end
        
        function resid = bc_vals(ya,yb)
           resid = [ya(3)-yb(3), ya(4), yb(4), ya(5), yb(5)-rel_flux, ya(6)-yb(6)]; 
        end
    end      

end