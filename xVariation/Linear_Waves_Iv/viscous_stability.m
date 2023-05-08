function stable = viscous_stability(theta,Fr,nu,lambda)
% Stability of the two equation model with viscosity. Can either run for a
% single wavelength or a range.
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    
    rho_p = 2500;
    rho_f = 1000;
    phi_c = 0.585;
    g=9.81;
    eta_f = 0.0010016;
    if ~exist("theta","var")
        specify_param = false;
    else
        specify_param = true;
    end
    if ~exist("lambda","var")
        single_lambda = false;
    else
        single_lambda = true;
    end
    if ~specify_param
        theta = 12;
        nu = 1.13e-4;
        Fr=0.8;
    end
    
    rho = rho_p*phi_c+(1-phi_c)*rho_f;
    crit_Iv = newt_solve_crit_Iv(theta,rho_p,rho_f);
    u_const = crit_Iv/eta_f/3*(rho_p-rho_f)*g*phi_c*cosd(theta);
    h0 = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3);  
    u_eq = u_const.*h0^2;
    nu_dl = nu/(u_eq*h0);
    
    P = (rho-rho_f)/rho;
    
    dIvdu = crit_Iv;
    dIvdh = -2*crit_Iv;

    dmudu = dmudIv_fn(crit_Iv).*dIvdu;
    dmudh = dmudIv_fn(crit_Iv).*dIvdh;
    
    if ~single_lambda
        npt = 100;
        max_imag = zeros(1,npt);
        lamda_val = linspace(8,12,npt);
        for i = 1:npt
            lambda = lamda_val(i);
            max_imag(i) = single_max_imag(lambda);
        end
        plot(lamda_val,max_imag)
        stable = [];
    else
        single_imag = single_max_imag(lambda);
        stable = single_imag<0;
    end
    
    function max_i = single_max_imag(lambda)
        k = 2*pi/lambda;
        c1 = -1;
        c2 = 2*k-1i/Fr^2*(P*dmudu)-1i*k^2*nu_dl;
        c3 = k^2*(1/Fr^2-1)+1i*k/Fr^2*(P*dmudu-P*dmudh)+1i*k^3*nu_dl;
        r = roots([c1, c2, c3]);
        max_i = max(imag(r));
    end
end