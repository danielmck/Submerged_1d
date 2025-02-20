function viscous_Iv_bvp_from_ode_yield
% Converts ode solution from viscous_wave_replica_Iv into a solution to the
% viscous two equation bvp.
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction

    g=9.81; % m/s^2
    theta = 12;

%     eta_f = 1.18e-5;
%     rho_f = 1;
    
    rho_f = 1000;
    eta_f = 0.0010016; % Pa s
    
    rho_p = 2500;
    
    Fr_eq = 5; 
    lambda = 10;
    tau0=0.1;
    
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
    
%     crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);
    
    [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr_eq, tau0);
    u_const = crit_Iv/eta_f/3*(rho_p-rho_f)*g*phi_c*cosd(theta);
    u_eq = u_const.*h0^2;
    
    nu = 1.13e-4;
    
    z_scale = h0;
    v_scale = u_eq;
    t_scale = z_scale/v_scale;

    nu_dl = nu/z_scale/v_scale;
    tau0_dl = tau0/(rho_f*g*cosd(theta)*h0);
    R = u_eq*h0/nu;
    
    % Gets the waveform from the ode solver
    [xi_wave, y_wave, uw_wave] = viscous_wave_replica_yield(theta, rho_f, rho_p, eta_f, nu, Fr_eq, lambda,tau0);
    Q1_wave = uw_wave - 1;
    lambda_init = xi_wave(end);
    lambda_ratio = lambda/lambda_init;
    n_step = min(max(ceil(30*(lambda_ratio-1)),10),300);
    lambda_med = linspace(lambda_init,lambda,n_step);
    
    uw_init_vec = uw_wave*ones(size(xi_wave));
    Q1_init_vec = Q1_wave*ones(size(xi_wave));
    u_init = uw_wave - Q1_wave./y_wave(:,2);
    m_init = zeros(size(xi_wave));
    m_val = 0;
    % Defines the average flux m that is needed to solve the bvp
    for j = 1:size(xi_wave,1)
        m_init(j) = m_val;
        if j < size(xi_wave,1)
            m_val = m_val + 1/lambda_init* (y_wave(j,2)*u_init(j)+y_wave(j+1,2)*u_init(j))/2*(xi_wave(j+1)-xi_wave(j));
        end
    end
    % Creates the initial vector including wave speed and the mass balance
    % constant Q1
    y_init = vertcat(uw_init_vec',Q1_init_vec',y_wave(:,2:3)',m_init');
    
    [xi_final,y_final] = run_bvp_step(lambda_init, lambda, n_step, xi_wave, y_init,1e-3);
    h_final = y_final(3,:);
    plot(xi_final,h_final)
    [~,mindex] = min(h_final);
%     y_final = horzcat(y_final(:,mindex:end),y_final(:,1:mindex-1));
%     xi_final = mod(horzcat(xi_final(mindex:end),xi_final(1:mindex-1))-xi_final(mindex),lambda);
    out_final = vertcat(xi_final,y_final);
    filename = "master_wave_yield.txt";
    save(filename,"out_final","-ascii")
    write_record("Results/wave_record.csv",filename,{"no_pe","water",Fr_eq,theta,lambda,nu,0,0,tau0})
    
    function [xi_out, y_out] = run_bvp_step(lambda_init, lambda_fin, nstep, xi_in, y_in, tol)
        % Can account for a change in wavelength but should really use
        % viscous_Iv_bvp_from_master for that.
        lambda_in = lambda_init;
        for i = 1:nstep
            lambda_old = lambda_in;
            lambda_in = lambda_init + (lambda_fin-lambda_init)/n_step*i;
            xi_in = xi_in/lambda_old*lambda_in;
            solInit1=bvpinit(xi_in,@bvp_guess);
            try
                solN1 = bvp4c(@viscous_syst,@bc_vals,solInit1);
            catch ME
                switch ME.identifier
                    case 'MATLAB:UndefinedFunction'
                        warning('Function is undefined.  Assigning a value of NaN.');
                    otherwise
                        rethrow(ME)
                end
                
            end
            % Solves the 5 value bvp for u_w, Q1, h, n and m.
            resid = solN1.stats.maxres;
            h_wave = solN1.y(3,:);
            h_diff = max(h_wave)-min(h_wave);
            if h_diff>1e-4
                if resid < tol
                    y_in = solN1.y;
                    xi_in = solN1.x;
                else
                    [xi_in,y_in] = run_bvp_step(lambda_old, lambda_in, 2, xi_in*lambda_old/lambda_in, y_in, tol);
                end
            else
                [xi_in,y_in] = run_bvp_step(lambda_old, lambda_in, 2, xi_in*lambda_old/lambda_in, y_in, tol);
            end
        end
        y_out = solN1.y;
        xi_out = solN1.x;
        
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
            h = y(3);
            u = u_w - Q1/h;
            n = y(4);
            m = y(5);

            dhdxi = n;
            n_coeff = 1-Q1^2.*Fr_eq^2/h^3;
            Fr = Fr_eq*abs(u)/h;
            Iv = crit_Iv*abs(u)/h^2;
            n_eq = (tand(theta)-sign(u).*P*mu_Iv_fn(Iv)-tau0_dl*rho_f/rho/h)./n_coeff;
            dndxi = 1/(2*h)*n^2 + h^(3/2)/Fr_eq^2*R/Q1*n_coeff*(n-n_eq);
            dmdxi = h/lambda_in*u;
            dydxi = [0,0,dhdxi,dndxi,dmdxi]';
        end
    end
    
    function resid = bc_vals(ya,yb)
       resid = [ya(3)-yb(3), ya(4), yb(4), ya(5), yb(5)-1]; 
    end

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
end