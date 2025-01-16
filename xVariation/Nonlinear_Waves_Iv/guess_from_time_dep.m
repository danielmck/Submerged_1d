function guess_from_time_dep
    guess_file = "bvp_guess.txt";
    td_guess=load("Results/"+guess_file);
    h_deriv = td_guess(4,:);
    npts=size(h_deriv,2);
    
    Fr = 3.0;
    theta = 12;
    tau0 = 0;
    lambda=300;
    rf_val = 1;
    alpha=1e-5;
    d=1e-4;
    pres_h=true;
    xi_guess = linspace(0,lambda,npts);
    
    [phi_c,rho_f,rho_p,~,eta_f,g] = get_params_water();
    [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0,false,true);
    u_eq = Fr*sqrt(g*cosd(theta)*h0);
    phi_eq = phi_c/(1+sqrt(crit_Iv));
    rho_eq = phi_eq*rho_p+(1-phi_eq)*rho_f;
    crit_pb = rho_f*g*cosd(theta)*h0;
    nu = 3/4*mu_Iv_fn(crit_Iv)*eta_f/crit_Iv/rho_eq*10;
    nu_dl = nu/(u_eq*h0);

%                 h_stop = tau_0/rho*g*cosd(theta)/(rho/(rho-rho_f)-mu1_Iv);

    z_scale = h0;
    v_scale = u_eq;
    p_scale = crit_pb;
    t_scale = z_scale/v_scale;

    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    g_dl = g*t_scale/v_scale; 

    tau0_dl = tau0/p_scale;
%                 p_tot_grad_dl = p_tot/p_scale*z_scale;
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale; 
    rho_eq = rho_p_dl*phi_eq+rho_f_dl*(1-phi_eq);
    chi_eq = (rho_f_dl+3*rho_eq)/(4*rho_eq);
    d_dl = d/z_scale;

    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
    tol = 1e-3;
    % Solves the stepped system
    solInit1=bvpinit(xi_guess,@bvp_guess);
    opts = bvpset('RelTol',tol,'NMax',5000);
    global_max=0;
    for j=2:npts
        guess_grad = (td_guess(:,j)-td_guess(:,j-1))/(xi_guess(j)-xi_guess(j-1));
        fn_grad = (full_syst(0,td_guess(:,j))+full_syst(0,td_guess(:,j-1)))/2;
        grad_diff = guess_grad-fn_grad;
        max_diff = max(abs(grad_diff));
        global_max=max(global_max,max_diff);
    end
    [t_ode,odesoln] = ode15s(@full_syst,[0,300],td_guess(:,1));
    solN1 = bvp4c(@full_syst,@bc_vals,solInit1,opts);
    resid = solN1.stats.maxres;
    
    function guess = bvp_guess(xi)
        guessd = sum(xi_guess<xi)+1;
        if guessd == max(size(xi_guess))
            guess = td_guess(:,end);
        else
            gap = xi_guess(guessd+1)-xi_guess(guessd);
            guess = td_guess(:,guessd)*(xi-xi_guess(guessd))/gap + td_guess(:,guessd+1)*(xi_guess(guessd+1)-xi)/gap;
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
        rho_dl = (rho_p_dl*phi+rho_f_dl*(1-phi));
        p_tot_grad_dl = rho_dl*g_dl*cosd(theta);
        chi = (rho_f_dl+3*rho_dl)/(4*rho_dl);
        pb = y(7) + rho_dl*g_dl*cosd(theta)*chi.*h;
        P = (rho_dl-rho_f_dl)/rho_dl;

        zeta = 3/(2*alpha_dl*h) + P/4;
        p_p = p_tot_grad_dl.*h-pb;
        D = -2/beta_dl/h*(pb-h);
        Iv = 3*eta_f_dl*abs(u)/h/p_p;
        R_w3 = -phi*rho_f_dl/rho_dl*D;
        R_w4 = (-P/4+zeta)*D - 9/2/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));

        dhdxi = n;
        n_coeff = h/Fr^2-Q1^2/h^2;


        mu_val = p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv); 
        force_bal = h*(tand(theta)-sign(u).*mu_val)/Fr^2+tau0_dl*rho_f_dl/rho_dl+u*P*D;
        n_eq = (force_bal)./n_coeff;
        dQdxi = -P*D;

%             dndxi = 1/(2*h)*n^2 + h^(3/2)/Fr^2/nu_dl/Q1*n_coeff*(n-n_eq);
        dmdxi = rho_dl*h/lambda/rho_eq*u^(1-pres_h);

        dy6dxi = -R_w3;
        dy7dxi = R_w4/(u-u_w);
        dphidxi = (dy6dxi-phi*dQdxi)/Q1;
        dPdxi = rho_f_dl*(rho_p_dl-rho_f_dl)./rho_dl^2*dphidxi;
        dpbdxi = dy7dxi + rho_dl*g_dl*cosd(theta)*chi*dhdxi + 3/4*g_dl*cosd(theta)*(rho_p_dl-rho_f_dl)*dphidxi;
        dDdxi = -2/beta_dl*(dpbdxi/h-pb/h^2*dhdxi);
%             dndxi_old = 1/Q1/nu_dl*(n_coeff*n-force_bal+u*dQdxi)-h/Q1*(P*dDdxi+D*dPdxi);
        dndxi = 1/Q1/nu_dl*((n_coeff+nu_dl*P*(D-1))*n - force_bal+nu_dl*(h*D*dPdxi-P*2/beta_dl*dpbdxi));
        dydxi = [0,dQdxi,dhdxi,dndxi,dmdxi,dy6dxi,dy7dxi]';  
    end

    function resid = bc_vals(ya,yb)
        resid = [ya(3)-yb(3), ya(4), yb(4), ya(5), yb(5)-rf_val, ya(6)-yb(6), ya(7)-yb(7)]; 
    end
end