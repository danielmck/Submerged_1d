function [xi_final,y_final] = bvp_full_from_master_s(params,master_y,master_xi,master_params)
    % Uses the master solution to parameter step to the desired values
    % specified in params. If the params vector is not set then the default
    % values specified are used.
    
    % The parameters that can be set are Fr, theta, lambda, alpha and d.
    
    % If desired, a different master wave can be specified in order to
    % speed up computation. The master y, xi and parameter values must be
    % specified.
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction

    g=9.81; % m/s^2
    
%     eta_f = 1.18e-5;
%     rho_f = 1;
    
    rho_f = 1000;
    eta_f = 0.0010016; % Pa s
    rho_p = 2500;
    
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    chi = (rho_f+3*rho)/(4*rho);
    P = (rho-rho_f)/rho;
    
    % Defines parameters if not specified
    if ~exist('params','var')
        custom_init = false;
        Fr_eq = 0.8;
        theta = 8.35;
        lambda = 50;
        nu = 2e-4;
        alpha = 1e-5;
        d = 1e-4;
        tau0 = 0; % Pa
        s = 0.1;
        params = [Fr_eq,theta,lambda,nu,alpha,d,tau0];
        filename = "theta_8_35_s_point_1_lambda_50.txt";
    else
        custom_init = true;
        param_cell = num2cell(params);
        [Fr_eq,theta,lambda,nu,alpha,d] = param_cell{:};  
    end
    
    if ~exist('master_y','var')
        master_name = "theta_8_35_s_point_1.txt";
        master_file = load("Results/"+master_name);
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        record = readtable('Results/wave_record.csv');

        in_table = strcmp(record.Name, master_name);
        wave_type = record.wave_type(in_table);
        master_theta = record.theta(in_table); 
        master_lambda = record.lambda(in_table);
        master_Fr = record.Fr(in_table);
        master_nu = record.nu(in_table);
        master_tau0 = record.tau0(in_table);
        master_s=0.1;
        if strcmp(wave_type,'full')
            master_d = record.d(in_table);
            master_alpha = record.alpha(in_table);
        end
    else
        master_temp = num2cell(master_params);
        [master_Fr,master_theta,master_lambda,master_nu,master_alpha,master_d] = master_temp{:};
    end
    
%     
%     max_eig = zeros(1,100);
%     for lin_lambda = 1:100
%         lin_k = 2*pi/(0.1*lin_lambda);
%         A_mat = make_A_mat(lin_k,rho_p,rho_f,theta,eta_f,d,alpha,Fr_eq,newt_solve_crit_Iv(theta,rho_p,rho_f),nu);
%         A_eig = eigs(A_mat);
%         [~, idx] = sort(imag(A_eig),'descend');
%         A_eig = A_eig(idx);
%         max_eig(lin_lambda)=imag(A_eig(1));
%     end
    k = 2*pi/lambda;
    A_mat = make_A_mat(k,rho_p,rho_f,theta,eta_f,d,alpha,Fr_eq,newt_solve_crit_Iv(theta,rho_p,rho_f),nu);
    A_eig = eigs(A_mat);
    [~, idx] = sort(imag(A_eig),'descend');
    A_eig = A_eig(idx);
    % Checks that the parameters lead to linearly unstable wave perts
    if (imag(A_eig(1))<0)
        error("Wave is not unstable, try a different value")
    end
    Fr_ratio = max(Fr_eq/master_Fr,master_Fr/Fr_eq);
    nu_ratio = max(nu/master_nu,master_nu/nu);
    theta_ratio = max(theta/master_theta,master_theta/theta);
    lambda_ratio = max(lambda/master_lambda,master_lambda/lambda);
    alpha_ratio = 1+max(log(alpha/master_alpha),log(master_alpha/alpha));
    d_ratio = 1+max(log(d/master_d),log(master_d/d));
    tau0_ratio = abs(tau0-master_tau0)/(max(master_tau0,tau0)+1e-6);
    s_diff = abs(master_s-s);
    ratio_sum = Fr_ratio-1+nu_ratio-1+50*(theta_ratio-1)+35*(lambda_ratio-1)+d_ratio-1+alpha_ratio-1+tau0_ratio+s_diff;
    if ratio_sum > 1e-6
        n_steps = max(min(ceil(30*(ratio_sum-1)),400),90);
    else
        n_steps=2;
    end

    % Sets lists for parameter stepping between master and designated
    % values
    Fr_list = linspace(master_Fr,Fr_eq,n_steps);
    nu_list = linspace(master_nu,nu,n_steps);
    theta_list = linspace(master_theta,theta,n_steps);
    lambda_list = linspace(master_lambda,lambda,n_steps);
    alpha_list = logspace(log10(master_alpha),log10(alpha),n_steps);
    d_list = logspace(log10(master_d),log10(d),n_steps);
    tau0_list = linspace(master_tau0,tau0,n_steps);
    s_list = linspace(master_s,s,n_steps);

    [xi_final,y_final] = run_bvp_step(Fr_list, theta_list, lambda_list, nu_list, alpha_list, d_list, tau0_list, s_list, master_xi, master_y);
%     [xi_nope,y_nope] = viscous_Iv_bvp_from_master(true,[Fr_eq,theta,lambda,nu]);
    
    if ~custom_init
        u_w = y_final(1,1);
        Q1 = y_final(2,:);
        h = y_final(3,:);
        u = u_w - Q1./h;
        n = y_final(4,:);
        m = y_final(5,:);
        phi = y_final(6,:)./Q1;
        pb = y_final(7,:) + rho/rho_f*chi.*h;
        out_vec = vertcat(xi_final,y_final);
        
        write_record("Results/wave_record.csv",filename,{"full","water",Fr_eq,theta,lambda,nu,alpha,d,tau0})
        save("Results/"+filename,"out_vec","-ascii")
%         plot(xi_out,u)

        hold on
        plot(xi_final,u)
%         plot(master_xi,master_y(3,:))
    end
    
    hold on
%     plot(xi_nope,y_nope(3,:))
    
    function [xi_out, y_out] = run_bvp_step(Fr_vals, theta_vals, lambda_vals, nu_vals, alpha_vals, d_vals, tau0_vals, s_vals, xi_in, y_in, tol,counter)
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
        if counter > 15
            out_vec = vertcat(xi_in,y_in);
            fail_name = "theta_8_35_lambda_long.txt";
            s_vals(1)
            write_record("Results/wave_record.csv",fail_name,{"full","water",Fr_vals(1),theta_vals(1),lambda_vals(1),nu_vals(1),alpha_vals(1),d_vals(1),tau0_vals(1)})
            save("Results/"+fail_name,"out_vec","-ascii")
            error("Max iteration depth reached, non convergence. Outputting final state to "+fail_name)
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
                s_in = s_vals(i);
                xi_run = xi_in/lambda_vals(i-1)*lambda_vals(i);
%                 xi_run = horzcat(linspace(0,3*lambda_vals(i)/4,70),linspace(3*lambda_vals(i)/4, lambda_vals(i),70));
                if tau0_in == 0
                    crit_Iv = newt_solve_crit_Iv(theta_in, rho_p, rho_f,s_in);
                    pp_eq_grad = (rho_p-rho_f)*g*phi_c*cosd(theta_in);
                    u_const = crit_Iv/eta_f/2*pp_eq_grad*(1-s_in);
                    h0 = ((Fr_in*sqrt(g*cosd(theta_in)))./u_const)^(2/3);  
                else
                    [h0, crit_Iv] = crit_Iv_tau0(theta_in, rho_p, rho_f, eta_f, Fr_in, tau0_in);
                end
                u_eq = Fr_in*sqrt(g*cosd(theta_in)*h0);
                phi_eq = phi_c/(1+sqrt(crit_Iv));
                p_tot = rho*g*cosd(theta_in);
                crit_pb = rho_f*g*cosd(theta_in)*h0;
                nu_dl = nu_vals(i)/(u_eq*h0);
                
%                 h_stop = tau_0/rho*g*cosd(theta)/(rho/(rho-rho_f)-mu1_Iv);
                
                z_scale = h0;
                v_scale = u_eq;
                p_scale = crit_pb;
                t_scale = z_scale/v_scale;

                eta_f_dl = eta_f/(p_scale*t_scale);
                alpha_dl = alpha_vals(i)*p_scale;
                g_dl = g*t_scale/v_scale; 

                tau0_dl = tau0_in/p_scale;
                
                p_tot_grad_dl = p_tot/p_scale*z_scale;
                rho_f_dl = rho_f*v_scale^2/p_scale;
                rho_p_dl = rho_p*v_scale^2/p_scale; 
                rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
                d_dl = d_vals(i)/z_scale;
                pb_eq_grad = 1+s_in*(rho_dl-rho_f_dl)*g_dl*cosd(theta);

                beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);

                % Solves the stepped system
                solInit1=bvpinit(xi_run,@bvp_guess);
                opts = bvpset('RelTol',tol,'NMax',200*counter);
                try
                    solN1 = bvp4c(@full_syst,@bc_vals,solInit1,opts);
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
                         linspace(nu_vals(i-1),nu_vals(i),3),linspace(alpha_vals(i-1)...
                        ,alpha_vals(i),3),linspace(d_vals(i-1),d_vals(i),3),linspace(tau0_vals(i-1),tau0_in,3)...
                        ,linspace(s_vals(i-1),s_in,3),xi_in, y_in, tol,counter+1);
            non_pe = false;
            if counter == 1
                if isequal(y_return,-1)                
                    curr_params = [Fr_vals(i-1),theta_vals(i),lambda_vals(i-1),nu_vals(i-1),alpha_vals(i-1),d_vals(i-1),tau0_vals(i-1)];
                    [xi_no_pe,y_no_pe] = bvp_non_pe_to_full(true,true,curr_params,true,xi_in,y_in);
                    y_no_pe = y_no_pe(1:5,:);
                    [xi_no_pe,y_no_pe] = viscous_Iv_bvp_from_master(true,params,true,xi_no_pe,y_no_pe,[curr_params(1:4),curr_params(7)]);
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
            
        function dydxi = full_syst(xi,y)
            u_w = y(1);
            Q1 = y(2);
            h = y(3);
            u = u_w - Q1/h;
            n = y(4);
            m = y(5);
            phi = y(6)./Q1;
            pb = y(7) + rho_dl*g_dl*cosd(theta_in)*chi.*h;

            zeta = 3/(2*alpha_dl*h) + P/4;
            p_p = p_tot_grad_dl*h-pb;
            p_eq = h*pb_eq_grad;
            D = -2/beta_dl/h*(pb-p_eq);
            Iv = abs(2*eta_f_dl*u/h/p_p);
            R_w3 = -phi*rho_f_dl/rho_dl*D;
            R_w4 = (-P*chi+zeta)*D - 2*3/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));

            dhdxi = n;
            n_coeff = 1-Q1.^2.*Fr_in^2/h^3;
            mu_val = p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv)+tau0_dl*rho_f/rho/h;
            n_eq = (tand(theta_in)-sign(u).*mu_val)./n_coeff;
            dQdxi = -P*D;
            dndxi = 1/(2*h)*n^2 + h^(3/2)/Fr_in^2/nu_dl/Q1*n_coeff*(n-n_eq);
            dmdxi = h/lambda_in*u;

            dy6dxi = -R_w3;
            dy7dxi = R_w4/(u-u_w);
            dydxi = [0,dQdxi,dhdxi,dndxi,dmdxi,dy6dxi,dy7dxi]';
        end
        
        function resid = bc_vals(ya,yb)
            resid = [ya(3)-yb(3), ya(4), yb(4), ya(5), yb(5)-1, ya(6)-yb(6), ya(7)-yb(7)]; 
        end
    end

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
end