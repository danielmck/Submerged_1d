function [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water()
    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    eta_f = 0.0010016; % Pa s
    rho_f = 1000; % kg/m^3
    rho_p = 2500; % kg/m^
    rho = rho_p*phi_c+rho_f*(1-phi_c);
end