[phi_c,rho_f,rho_p,~,eta_f,g] = get_params_water();
theta = 12;
tau0 = 0;
d = 1e-4;
alpha = 1e-5;

h_init = 0.0061;
h0 = 7.2822e-04;
u0 = 0.2128;
pb0 = 4.4145;
phi0 = 0.4328;

[Fr,crit_Iv] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h_init, tau0, false,true);

% twoeqn_out = two_eqn_stab(h_init,theta, rho_p, rho_f, eta_f, tau0);
% foureqn_out = single_Fr_stab(Fr,crit_Iv,theta, rho_p, rho_f, d, eta_f, alpha);
noneq_out = Fr_stab_noneq(h0,u0,phi0,pb0,theta, rho_p, rho_f, d, eta_f, alpha,tau0);
