function [xi_out,y_out] = bvp_non_pe_to_rauter(custom_init,reverse,params,provide_init,master_xi,master_y)
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
            [Fr,theta,lambda,nu,a,d,tau0] = param_cell{:};
            rel_flux = 1;
            pres_h = 0;
        elseif (size(param_cell,2) == 9)
            [Fr,theta,lambda,nu,a,phi_rlp,phi_rcp,d,tau0,rel_flux,pres_h] = param_cell{:};
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
                master_name="master_wave_full.txt";
            else
                master_name = "rauter_convert.txt";
%                 master_name = "time_d_comp_vvshort.txt";
            end
            master_file = load(strcat("../Results/",master_name));
            master_xi = master_file(1,:);
            master_y = master_file(2:end,:);
        end
        record = readtable('../Results/wave_record.csv');
        in_table = strcmp(record.Name, master_name);
        wave_type = record.wave_type(in_table);
        theta = record.theta(in_table); 
        lambda = record.lambda(in_table);
        Fr = record.Fr(in_table);
        nu = record.nu(in_table);
        tau0 = record.tau0(in_table);
        rel_flux = master_y(5,end);
        if reverse
            a = record.alpha(in_table);
            d = record.d(in_table);
            phi_rlp = record.phi_rlp(in_table);
            phi_rcp = record.phi_rcp(in_table);
        else
            a = 130;
            d = 1e-4;
            phi_rcp = 0.63;
            phi_rlp = 0.53;
        end
        pres_h=0;
    end
    
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
    
    [~, init_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0);
    phi_eq = phi_c/(1+sqrt(init_Iv));
    
    delta_vals = logspace(-6,0,50);
    delta_vals = horzcat(0,delta_vals);
    if reverse
        delta_vals = wrev(delta_vals(2:end));
    end
    % Checks that the full case is unstable at the requested parameter
    % values
    if ~reverse
        y6_in = phi_eq*master_y(2,:);
        y_in = vertcat(master_y,y6_in);
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
        out_vec = vertcat(xi_out,y_out(2:end,:));
        filename = "rauter_wave_full.txt";
        save("Results/"+filename,"out_vec","-ascii")
        write_record("Results/full_record.csv",filename,{"full","water",Fr,theta,nu,a,phi_rlp,phi_rcp,d,tau0,y_out(1,1),lambda,0})
    end

    u_w = y_out(1,1);
    Q1 = y_out(2,:);
    h = y_out(3,:);
    u = u_w - Q1./h;
    n = y_out(4,:);
    m = y_out(5,:);
    if ~reverse
        phi = y_out(6,:)./Q1;
    end
    if ~custom_init
%         plot(xi_out,u)
        hold on
        plot(xi_out,h)
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
            [h0, crit_Iv] = crit_Iv_rauter(theta, rho_p, rho_f, eta_f, Fr, a, phi_rlp, phi_rcp, tau0, false, delta_in);
            u_eq = Fr*sqrt(g*cosd(theta)*h0);
            pp0 = 3*eta_f*u_eq/h0/crit_Iv; 
            rho_eq = pp0/h0/g/cosd(theta)+rho_f;
            phi_eq = (rho_eq-rho_f)/(rho_p-rho_f);
            
            p_tot = rho_eq*g*cosd(theta);
            crit_pb = rho_f*g*cosd(theta)*h0;
            nu_dl = nu/(u_eq*h0);

            z_scale = h0;
            v_scale = u_eq;
            p_scale = crit_pb;
            t_scale = z_scale/v_scale;

            tau0_dl = tau0/p_scale;
            eta_f_dl = eta_f/(p_scale*t_scale);
            a_dl = a/p_scale;
            g_dl = g*t_scale/v_scale; 

            u_eq_dl = u_eq/v_scale;
            p_tot_grad_dl = p_tot/p_scale*z_scale;
            rho_f_dl = rho_f*v_scale^2/p_scale;
            rho_p_dl = rho_p*v_scale^2/p_scale; 
            rho_eq_dl = rho_eq*v_scale^2/p_scale;
            d_dl = d/z_scale;

            beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
            solInit1=bvpinit(xi_in,@bvp_guess);
            opts = bvpset('RelTol',tol,'NMax',400*counter);
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
        
        function dydxi = rauter_syst(xi,y)
            u_w = y(1);
            Q1 = y(2);
            h = y(3);
            u = u_w - Q1/h;
            n = y(4);
            m = y(5);
            phi = y(6)./Q1;

            rho_dl = delta_in*(phi*rho_p_dl+(1-phi)*rho_f_dl)+(1-delta_in)*(phi_c*rho_p_dl+(1-phi_c)*rho_f_dl);
            rho_eq_merge = delta_in*rho_eq_dl+(1-delta_in)*(phi_c*rho_p_dl+(1-phi_c)*rho_f_dl);
            p_p_contact = max(a_dl*(phi-phi_rlp)/(phi_rcp-phi),0);
            Iv_phi = (phi_c/phi-1)^2;
            p_p_shear = 3*u*eta_f_dl/Iv_phi/h;
            p_p = delta_in*(p_p_contact+p_p_shear)+(1-delta_in)*(rho_dl-rho_f_dl)*g_dl*cosd(theta)*h;
            p_tot_grad_dl = rho_dl*g_dl*cosd(theta);
            pb = p_tot_grad_dl*h-p_p;
            D = -2/beta_dl/h*(pb-h);
            Iv = abs(3*eta_f_dl*u/h/p_p);
            R_w3 = -phi*rho_f_dl/rho_dl*D;

            dhdxi = n;
            n_coeff = 1-Q1.^2.*Fr^2/h^3;
            Iv = 3*eta_f_dl*abs(u)/h/p_p;
            mu_val = delta_in*p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv)+(1-delta_in)*P*mu_Iv_fn(crit_Iv*abs(u)/h^2)+tau0_dl*rho_f/rho/h; 
            force_bal = h*(tand(theta)-sign(u).*mu_val)/Fr^2+delta_in*u*P*D;
            n_eq = (force_bal)./n_coeff;
            dQdxi = -delta_in*P*D;
            
            dmdxi = h*rho_dl/rho_eq_merge/lambda*u^(1-pres_h);
            dy5dxi = delta_in*phi*rho_f_dl/rho_dl*D;
            dphidxi = (dy5dxi-phi*dQdxi)/Q1;
            dPdxi = delta_in*(rho_p_dl-rho_f_dl)*rho_f_dl/rho_dl^2*dphidxi;
            dppdh = 3*eta_f_dl/Iv_phi/h*(Q1/h^2)-3*u*eta_f_dl/Iv_phi/h^2;
            dppdQ = 3*eta_f_dl/Iv_phi/h*(-1/h);
            dppdphi = a_dl*(phi_rcp-phi_rlp)/(phi_rcp-phi)^2+6*u*eta_f_dl/h/(phi_c/phi-1)^3*phi_c/phi^2;
            dppdxi = dppdh.*dhdxi+dppdQ.*dQdxi+dppdphi.*dphidxi;
            dDdxi = -2/beta_dl*((rho_p_dl-rho_f_dl)*g_dl*cosd(theta)*dphidxi-dppdxi/h+p_p/h^2*dhdxi);
            dndxi = 1/Q1/nu_dl*(n_coeff*n-force_bal+delta_in*u*dQdxi)-delta_in*h/Q1*(P*dDdxi+D*dPdxi);
            dydxi = [0,dQdxi,dhdxi,dndxi,dmdxi,dy5dxi];
        end
        
        function dydxi = no_pe_syst(xi,y)
            % The system from the viroulet paper with a different friction
            % law. Need to solve for u_w and Q1 too.
            u_w = y(1);
            Q1 = y(2);
            h = y(3);
            u = u_w - Q1/h;
            n = y(4);
            m = y(5);

            dhdxi = n;
            n_coeff = h/Fr^2-Q1^2/h^2;
%             Fr = Fr*abs(u)/h;
            Iv = crit_Iv*abs(u)/h^2;
            force_bal = h*(tand(theta)-sign(u).*P*mu_Iv_fn(Iv))/Fr^2;
            n_eq = (force_bal)./n_coeff;
            dndxi = 1/Q1/nu_dl*(n_coeff*n-force_bal);
            dmdxi = h*u^(1-pres_h)/lambda;
            dydxi = [0,0,dhdxi,dndxi,dmdxi]';
        end        
        
        function resid = bc_vals(ya,yb)
           resid = [ya(3)-yb(3), ya(4), yb(4), ya(5), yb(5)-rel_flux,ya(6)-yb(6)]; % 
        end
    end   
end