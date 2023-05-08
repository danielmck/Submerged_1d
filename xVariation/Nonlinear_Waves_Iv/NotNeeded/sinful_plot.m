filename = "no_vis_full_master.txt"; %long_low_pe
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
s=0;


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
    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f,s);
    pp_eq_grad = (rho_p-rho_f)*g*phi_c*cosd(theta);
    u_const = crit_Iv/eta_f/2*pp_eq_grad*(1-s);
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
h_crit_pos = y(3,1);

Q1 = horzcat(y(4,:),y(9,2:end));
h = horzcat(y(5,:),y(10,2:end));
h(size(xi)) = (y(5,end)+y(10,1))/2;
u = u_w - Q1./h;
m = horzcat(y(6,:),y(11,2:end));
phi = horzcat(y(7,:),y(12,2:end))./Q1;
pb = horzcat(y(8,:),y(13,2:end)) + rho_dl*g_dl*cosd(theta)*chi.*h;

xi = horzcat(xi*h_crit_pos,h_crit_pos+xi(2:end)*(lambda-h_crit_pos));
