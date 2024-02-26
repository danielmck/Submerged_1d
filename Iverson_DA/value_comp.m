h0 = 0.3; % layer height (m)
[phi_c,rho_f,rho_p,~,eta_f,g] = get_params_water();
theta = 5; % deg
alpha = 1e-5;
d=1e-5; % grain diameter (m)

Fr0 = 10;
phi0 = 0.58;
rho0 = rho_p*phi0+(1-phi0)*rho_f;
u0 = Fr0*sqrt(g*cosd(theta)*h0);
p_tot = rho0*g*cosd(theta)*h0;

Iv_init = 3*u0*eta_f/(phi0*(rho_p-rho_f)*g*cosd(theta)*h0);
Iv_eq = (phi_c-phi0)^2/phi_c^2;

kappa = ((1-phi_c).^3.*d^2)./(150*phi_c.^2);
beta = eta_f/kappa;

pp_fast_ss = 1/Iv_eq*3*u0*eta_f/h0;
pb_fast_ss = p_tot-pp_fast_ss;
ts_fast = 2*eta_f*alpha/(phi0*Iv_eq^(3/2));

pp_med_ss = g*sind(theta).*rho0.*h0./mu_Iv_fn(Iv_eq);
pb_med_ss = p_tot-pp_med_ss;
u_med_ss = pp_med_ss.*Iv_eq./3/eta_f*h0;
ts_medium = Iv_eq*rho0*h0^2/(3*eta_f*mu_Iv_fn(Iv_eq));

ts_slow = alpha*h0^2*beta;

t_scale = sqrt(h0/g);
ts_fast_dl = ts_fast/t_scale;
ts_medium_dl = ts_medium/t_scale;
ts_slow_dl = ts_slow/t_scale;

gd = 3*u0/h0;
Rey = gd*d^2*rho_f/eta_f;
Sav = gd^2*d^2*rho_p/((rho0-rho_f)*phi0*g*cosd(theta)*h0*mu_Iv_fn(Iv_eq));

beta_dim = beta*sqrt(g*h0^3)*alpha;
alpha_dim = alpha*rho0*g*cosd(theta)*h0;
pp_dim = (rho0-rho_f)*phi0*g*cosd(theta)*h0/(rho0*phi0*g*cosd(theta)*h0);
eta_dim = rho0*sqrt(g*h0)/eta_f;

t_dil_i = ts_fast/ts_medium;
t_i_diff = ts_medium/ts_slow;
t_dil_diff = ts_fast/ts_slow;

hold on
case_num = 3;
% SetPaperSize(10,10)
% scatter(case_num,ts_fast,'d','MarkerEdgeColor','#fc8d62','MarkerFaceColor','#fc8d62','DisplayName','$t_{dil}$') % 
% scatter(case_num,ts_medium,'o','MarkerEdgeColor','#8da0cb','MarkerFaceColor','#8da0cb','DisplayName','$t_{i}$')
% scatter(case_num,ts_slow,'^','MarkerEdgeColor','#abd9e9','MarkerFaceColor','#abd9e9','DisplayName','$t_{diff}$')

% scatter(case_num,ts_fast,'d','MarkerEdgeColor','#fc8d62','MarkerFaceColor','#fc8d62','HandleVisibility','off') % 
% scatter(case_num,ts_medium,'o','MarkerEdgeColor','#8da0cb','MarkerFaceColor','#8da0cb','HandleVisibility','off')
% scatter(case_num,ts_slow,'^','MarkerEdgeColor','#abd9e9','MarkerFaceColor','#abd9e9','HandleVisibility','off')