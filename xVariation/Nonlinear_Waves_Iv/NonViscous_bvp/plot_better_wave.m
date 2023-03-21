filename = "no_vis_better_master.txt"; %long_low_pe
master_file = load("Results/"+filename);
xi = master_file(1,:);
y = master_file(2:end,:);
record = readtable('Results/wave_record.csv');

in_table = strcmp(record.Name, filename);
wave_type = record.wave_type(in_table);
theta = record.theta(in_table); 
Fr = record.Fr(in_table);
tau0 = record.tau0(in_table);
if strcmp(wave_type,"full")
    full_model=true;
    d = record.d(in_table);
    alpha = record.alpha(in_table);
else
    full_model=false;
    alpha=0;
    d=0;
end

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

rho_p = 2500;
rho = rho_p*phi_c+rho_f*(1-phi_c);

chi = (rho_f+3*rho)/(4*rho);
P = (rho-rho_f)/rho;
s_c = 1-rho/(rho-rho_f)*tand(theta)/mu1_Iv;

if tau0 == 0
    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);
    pp_eq_grad = (rho_p-rho_f)*g*phi_c*cosd(theta);
    u_const = crit_Iv/eta_f/2*pp_eq_grad;
    h0 = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3);  
else
    [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0);
end
u_eq = Fr*sqrt(g*cosd(theta)*h0);
phi_eq = phi_c/(1+sqrt(crit_Iv));
p_tot = rho*g*cosd(theta);
crit_pb = rho_f*g*cosd(theta)*h0;

h_stop = tau0/(rho*g*cosd(theta))/(tand(theta)-(rho-rho_f)/rho*mu1_Iv);


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
tau0_dl = tau0/p_scale;
h_stop_dl = h_stop/z_scale;

beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);

u_w = y(1,1);
lambda = y(2,1);
h_crit_xi = y(3,1);
Q1 = y(4,:);
h = y(5,:);
u = u_w - Q1./h;
m = y(6,:);

xi = horzcat(xi(xi<1)*h_crit_xi,h_crit_xi+(xi(xi>=1)-1)*(lambda-h_crit_xi));

if full_model
    phi = y(7,:)./Q1;
    pb = y(8,:) + rho_dl*g_dl*cosd(theta)*chi.*h;
    pe = pb-h;
else
    pb=h;
end

Fr_vals = Fr.*u./sqrt(h);

h_min = roots([1,0,-u_w,-Q1(1)]);
h_static = Q1(1)/u_w;

[p_max,p_max_ind] = max(pb);
[p_min,p_min_ind] = min(pb);

[h_max,h_max_ind] = max(h);
[h_min,h_min_ind] = min(h);

h_grad = (h(2:end)-h(1:end-1))./(xi(2:end)-xi(1:end-1));

if full_model
    [phi_max,phi_max_ind] = max(phi);
    [phi_min,phi_min_ind] = min(phi);
end

pp = p_tot_grad_dl.*h-pb;
D = -2/beta_dl./h.*(pb-h);
Iv = abs(2*eta_f_dl.*u./h./pp);

if full_model
    zeta = 3./(2*alpha_dl.*h) + P/4;
    tan_psi = phi - phi_c./(1+sqrt(Iv));
    R_w3 = -phi.*rho_f_dl/rho_dl.*D;
    R_w4 = (-P.*chi+zeta).*D - 2*3/alpha_dl./h.*u.*(tan_psi);
end

Fr_equi = zeros(size(h));
Iv_equi = zeros(size(h));
for i=1:size(h,2)
    [Fr_equi(i),Iv_equi(i)] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h(i)*h0, tau0,0);
end

n_coeff = 1-Q1.^2.*Fr^2./h.^3;
Iv = 3*eta_f_dl.*abs(u)./h./pp;
mu_val = pp./(p_tot_grad_dl.*h).*mu_Iv_fn(Iv)+tau0_dl*rho_f/rho./h;
force_bal = tand(theta)-sign(u).*mu_val;
dhdxi = force_bal./n_coeff;

mu_val_min = (rho-rho_f)/rho.*mu1_Iv+tau0_dl*rho_f/rho./h;
force_bal_max = tand(theta)-sign(u).*mu_val_min;

dQdxi = -P.*D;
dmdxi = h./lambda.*u;

if full_model
    dy6dxi = -R_w3;
    dy7dxi = R_w4./(u-u_w);
    dpbdxi = dy7dxi + rho_dl*g_dl*cosd(theta)*chi.*dhdxi;

    dpbdxi_scale = dpbdxi/(p_max-p_min);
    dhdxi_scale = dhdxi/(h_max-h_min);
end

%%

C = viridis(3);
hold on
% plot(xi,2*Q1./Fr.^2.*dQdxi./(3.*h.^3), "DisplayName", "Waveform","color",C(1,:))
plot(xi(1:end-1),h_grad, "DisplayName", "Waveform","color",C(2,:))