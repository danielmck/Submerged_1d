function bvp_non_pe_to_full
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
    
    alpha = 1e-5;
    rho_p = 2500;
    d = 1e-4;
    
    master_file = load("master_wave_no_pe.txt");
    master_xi = master_file(1,:);
    master_y = master_file(2:end,:);
    master_cond = readtable("master_wave_no_pe_cond.csv");
    master_Fr = master_cond.Fr_eq;
    master_nu = master_cond.nu;
    master_theta = master_cond.theta;
    master_lambda = master_cond.lambda;
    
    Fr = master_Fr; 
    lambda = master_lambda;
    theta = master_theta;
    nu = master_nu;
    
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
    
    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);
    u_const = crit_Iv/eta_f/2*(rho_p-rho_f)*g*phi_c*cosd(theta);
    h0 = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3);  
    u_eq = u_const.*h0^2;
    phi_eq = phi_c/(1+sqrt(crit_Iv));
    p_tot = rho*g*cosd(theta);
    crit_pb = rho_f*g*cosd(theta)*h0;
    R = u_eq*sqrt(h0)/nu;
    
    z_scale = h0;
    v_scale = u_eq;
    p_scale = crit_pb;
    t_scale = z_scale/v_scale;

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
    delta_vals = logspace(-6,0,100);
    delta_vals = horzcat(0,delta_vals);
    
    y6_in = phi_eq*master_y(2,:);
    y7_in = master_y(3,:) - rho_dl*g_dl*cosd(theta)*chi.*master_y(3,:);
    
    k = 2*pi/lambda;
    
    A_mat = make_A_mat(k,rho_p,rho_f,theta,eta_f,d,alpha,Fr,crit_Iv);
    A_eig = eigs(A_mat);
    [~, idx] = sort(imag(A_eig),'descend');
    A_eig = A_eig(idx);
    
    if (imag(A_eig(1))<0)
        error("Wave is not unstable, try a different value")
    end
    y_in = vertcat(master_y,y6_in,y7_in);
    [xi_out,y_out] = run_bvp_iter(delta_vals, master_xi, y_in);
    out_vec = vertcat(xi_out,y_out);
    save("master_wave_full.txt","out_vec","-ascii")
    u_w = y_out(1,1);
    Q1 = y_out(2,:);
    h = y_out(3,:);
    u = u_w - Q1./h;
    n = y_out(4,:);
    m = y_out(5,:);
    phi = y_out(6,:)./Q1;
    pb = y_out(7,:) + rho_dl*g_dl*cosd(theta)*chi.*h;
    plot(xi_out,u)
    hold on
    plot(xi_out,h)
    
    function [xi_in,y_in] = run_bvp_iter(delta_vals, xi_in, y_in, tol,counter)
        if ~exist('counter','var')
            counter = 1;
        end
        if ~exist('tol','var')
            tol = 1e-3;
        end
        if counter > 10
            error("Max iteration depth reached, non convergence")
        end
        n_step = size(delta_vals,2);
        for i = 2:n_step
            delta_in = delta_vals(i);
            solInit1=bvpinit(xi_in,@bvp_guess);
            solN1 = bvp4c(@viscous_syst,@bc_vals,solInit1);
            resid = solN1.stats.maxres;
            if resid < tol
                h_wave = solN1.y(3,:);
                h_diff = max(h_wave)-min(h_wave);
                if (h_diff>1e-4)
                    y_in = interp1(solN1.x,solN1.y',master_xi)';
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
            pb = y(7) + rho_dl*g_dl*cosd(theta)*chi.*h;

            zeta = 3/(2*alpha_dl*h) + P/4;
            p_p = p_tot_grad_dl*h-pb;
            D = -2/beta_dl/h*(pb-rho_f_dl*g_dl*cosd(theta)*h);
            Iv = abs(2*eta_f_dl*u/h/p_p);
            R_w3 = -phi*rho_f_dl/rho_dl*D;
            R_w4 = (-P*rho_dl*g_dl*cosd(theta)*chi+zeta)*D - 2*3/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));

            dhdxi = n;
            n_coeff = 1-Q1.^2.*Fr^2/h^3;
            Iv = 2*eta_f_dl*abs(u)/h/p_p;
            mu_val = delta_in*p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv)+(1-delta_in)*P*mu_Iv_fn(crit_Iv*abs(u)/h^2);
            n_eq = (tand(theta)-sign(u).*mu_val)./n_coeff;
            dQdxi = delta_in*P*D;
            dndxi = 1/(2*h)*n^2 + h^(3/2)/Fr^2*R/Q1*n_coeff*(n-n_eq);
            dmdxi = h/lambda*u;

            dy6dxi = delta_in*R_w3;
            dy7dxi = delta_in*R_w4/(u-u_w);
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
       resid = [ya(3)-yb(3), ya(4), yb(4), ya(5), yb(5)-1, ya(6)-yb(6), ya(7)-yb(7)]; 
    end

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
end