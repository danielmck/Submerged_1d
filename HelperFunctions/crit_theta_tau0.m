 function [theta, Iv] = crit_theta_tau0(h0, rho_p, rho_f, eta_f, Fr_eq, tau0, dl, rho_var, phi_param)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    
    if ~exist("dl","var")
        dl = false;
    end
    if ~exist("rho_var","var")
        rho_var = false;
    end
    if ~exist("phi_param","var")
        phi_param = rho_var;
    end
    
    tol = 1*10^(-7);
    count = 0;
    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    theta_min = atan(mu1_Iv*(rho_p*phi_c+rho_f*(1-phi_c)));
    theta_max = 30;
    Fr = -1;
    while abs(Fr-Fr_eq)>tol
        theta = (theta_max+theta_min)/2;
        [Fr, ~] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h0, tau0, dl, rho_var, phi_param);
        if Fr>Fr_eq
            theta_max = theta;
        else
            theta_min = theta;
        end
    end
    theta = (theta_max+theta_min)/2;
    [~, Iv] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h0, tau0, dl, rho_var, phi_param);
 end