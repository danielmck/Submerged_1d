function viscous_Iv_bvp_from_ode
% Converts ode solution from viscous_wave_verify into a solution to the
% viscous two equation bvp. Verifies with the work from Gray and Edwards.
% However, u_w and Q1 do not match up exactly as they are not linked by
% Q1=u_w-1 but are independent.
% Works as of 20/10/22. Don't change!

    mu1 = tand(20.9);
    mu2 = tand(32.76);
    beta = 0.136;
    L = 8.25e-4;
    

    g=9.81; % m/s^2
    theta = 29;
    gamma = (mu2-tand(theta))/(tand(theta)-mu1);
    Fr_eq = 1.02; 
    nu = 1.13e-3;
    
    lambda = 32;
    
    h0 = Fr_eq*L*gamma/beta;
    u_eq = Fr_eq*sqrt(g*h0*cosd(theta));
    
    z_scale = h0;
    v_scale = u_eq;
    t_scale = z_scale/v_scale;
    

    g_dl = g*t_scale/v_scale; 
    nu_dl = nu/z_scale/v_scale;
    L_dl = L/z_scale;
    R = u_eq*sqrt(h0)/nu;
    
    % Gets the waveform from the ode solver
    [xi_wave, y_wave, uw_wave] = viscous_wave_verify(theta, Fr_eq, nu);
    Q1_wave = uw_wave - 1;
    lambda_init = xi_wave(end);
    lambda_ratio = lambda/lambda_init;
    n_step = min(max(ceil(30*(lambda_ratio-1)),10),300);
    
    uw_init_vec = uw_wave*ones(size(xi_wave));
    Q1_init_vec = Q1_wave*ones(size(xi_wave));
    u_init = uw_wave - Q1_wave./y_wave(:,1);
    m_init = zeros(size(xi_wave));
    m_val = 0;
    % Defines the average flux m that is needed to solve the bvp
    for j = 1:size(xi_wave,1)
        m_init(j) = m_val;
        if j < size(xi_wave,1)
            m_val = m_val + 1/lambda_init* (y_wave(j,1)*u_init(j)+y_wave(j+1,1)*u_init(j))/2*(xi_wave(j+1)-xi_wave(j));
        end
    end
    % Creates the initial vector including wave speed and the mass balance
    % constant Q1 
    y_init = vertcat(uw_init_vec',Q1_init_vec',y_wave',m_init');
    
    [xi_final,y_final] = run_bvp_step(lambda_init, lambda, n_step, xi_wave, y_init,1e-3);
    h_final = y_final(3,:);
    plot(xi_final,h_final)
    [~,mindex] = min(h_final);
    out_final = vertcat(xi_final,y_final);
    
    function [xi_out, y_out] = run_bvp_step(lambda_init, lambda_fin, nstep, xi_in, y_in, tol)
        % Can account for a change in wavelength but should really use
        % viscous_Iv_bvp_from_master for that.
        lambda_in = lambda_init;
        for i = 1:nstep
            lambda_old = lambda_in;
            lambda_in = lambda_init + (lambda_fin-lambda_init)/n_step*i;
            xi_in = xi_in/lambda_old*lambda_in;
            solInit1=bvpinit(xi_in,@bvp_guess);
            solN1 = bvp4c(@viscous_syst,@bc_vals,solInit1);
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
            
            dy2dx = n;
            n_coeff = 1-Q1^2.*Fr_eq^2/h^3;

            mu_b = mu1 + (mu2-mu1)*(1-u_w+u_w*h)/(1-u_w+u_w*h+h^(5/2)*gamma);
            n_eq = (tand(theta)-mu_b)./n_coeff;

            dy3dx = 1/(2*h)*n^2 + h^(3/2)/Fr_eq^2*R/Q1*n_coeff*(n-n_eq);
            dmdxi = h/lambda_in*u;
            dydxi = [0,0,dy2dx,dy3dx,dmdxi]';
        end
    end
    
    function resid = bc_vals(ya,yb)
       resid = [ya(3)-yb(3), ya(4), yb(4), ya(5), yb(5)-1];
    end

    
end