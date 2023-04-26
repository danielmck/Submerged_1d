function [xi_out,y_out] = bvp_full_to_var_rho(custom_init,reverse,params,provide_init,master_xi,master_y)
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
            master_name="master_wave_full.txt";
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
        nu = record.nu(in_table);
        tau0 = record.tau0(in_table);
        rel_flux = master_y(5,end);
        if reverse
            alpha = record.alpha(in_table);
            d = record.d(in_table);
        else
            alpha = 1e-5;
            d = 1e-4;
        end
        pres_h=strcmp(wave_type(max(end-7,1):end),"_pres_h");
    end
    
    
    
    delta_vals = linspace(0,1,50);
    delta_vals = horzcat(0,delta_vals);
    if reverse
        delta_vals = wrev(delta_vals(2:end));
    end
    % Checks that the full case is unstable at the requested parameter
    % values

    [xi_out,y_out] = run_bvp_iter(delta_vals, master_xi,master_y);
    if ~custom_init
        out_vec = vertcat(xi_out,y_out);
        filename = "master_wave_full_var_rho.txt";
        save("Results/"+filename,"out_vec","-ascii")
        write_record("Results/wave_record.csv",filename,{"full","water",Fr,theta,lambda,nu,alpha,d,tau0})
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
            phi_param = delta_vals(i);

            if tau0 == 0
                crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f,true,phi_param);
                u_const = crit_Iv/eta_f/3*(rho_p-rho_f)*g*phi_c*cosd(theta);
                h0 = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3);  
            else
                [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0,false,true,phi_param);
            end
            
            u_eq = Fr*sqrt(g*cosd(theta)*h0);
            phi_eq = phi_c/(1+sqrt(crit_Iv));
            crit_pb = rho_f*g*cosd(theta)*h0;
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
            rho_f_dl = rho_f*v_scale^2/p_scale;
            rho_p_dl = rho_p*v_scale^2/p_scale; 
            d_dl = d/z_scale;

            beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
            
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
            rho_dl = phi_param*(rho_p_dl*phi+rho_f_dl*(1-phi))+(1-phi_param)*(rho_p_dl*phi_c+rho_f_dl*(1-phi_c));
            P = (rho_dl-rho_f_dl)/rho_dl;
            chi = (rho_f_dl+3*rho_dl)/(4*rho_dl);
            p_tot_grad_dl = rho_dl*g_dl*cosd(theta);
            pb = y(7) + rho_dl*g_dl*cosd(theta)*chi.*h;
            
            
            zeta = 3/(2*alpha_dl*h) + P/4;
            p_p = p_tot_grad_dl*h-pb;
            D = -2/beta_dl/h*(pb-h);
            Iv = abs(2*eta_f_dl*u/h/p_p);
            R_w3 = -phi*rho_f_dl/rho_dl*D;
            R_w4 = (rho_f_dl*g_dl*cosd(theta)*(-P/4)*(phi_param)-rho_dl*g_dl*cosd(theta)*P*chi*(1-phi_param)+zeta)*D - 2*3/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));

            dhdxi = n;
            n_coeff = 1-Q1.^2.*Fr^2/h^3;
            Iv = 3*eta_f_dl*abs(u)/h/p_p;
            mu_val = p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv)+tau0_dl*rho_f_dl/rho_dl/h; 
            n_eq = (tand(theta)-sign(u).*mu_val+P*D)./n_coeff;
            dQdxi = -P*D;
            dndxi = 1/(2*h)*n^2 + h^(3/2)/Fr^2/nu_dl/Q1*n_coeff*(n-n_eq);
            dmdxi = h/lambda*u^(1-pres_h);

            dy6dxi = R_w3;
            dy7dxi = R_w4/(u-u_w);
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