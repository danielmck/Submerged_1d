filename = "Fr_point6_lambda_70.txt"; %long_low_pe
master_file = load("Results/"+filename);
xi = master_file(1,:);
y = master_file(2:end,:);
record = readtable('Results/wave_record.csv');

in_table = strcmp(record.Name, filename);
wave_type = record.wave_type(in_table);
wave_type = wave_type{1};
theta = record.theta(in_table); 
lambda = record.lambda(in_table);
Fr = record.Fr(in_table);
nu = record.nu(in_table);
tau0 = record.tau0(in_table);
% if strcmp(extract(wave_type{1,1},1),"full")
if (size(wave_type,2)>5)
    full_model= all(wave_type(1:7)=='var_rho');
    var_rho = all(wave_type(1:7)=='var_rho');
else
    full_model=0;
    var_rho = 0;
end
d = record.d(in_table);
alpha = record.alpha(in_table);
% else
%     full_model=false;
%     alpha=0;
%     d=0;
% end
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

[h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0, false, var_rho);

u_eq = Fr*sqrt(g*cosd(theta)*h0);
phi_eq = phi_c/(1+sqrt(crit_Iv));
p_tot = rho*g*cosd(theta);
crit_pb = rho_f*g*cosd(theta)*h0;
nu_dl = nu/(u_eq*h0);

h_stop = tau0/(rho*g*cosd(theta))/(tand(theta)-(rho-rho_f)/rho*mu1_Iv);


z_scale = h0;
v_scale = u_eq;
p_scale = crit_pb;
t_scale = z_scale/v_scale;

eta_f_dl = eta_f/(p_scale*t_scale);
alpha_dl = alpha*p_scale;
g_dl = g*t_scale/v_scale; 

u_eq_dl = u_eq/v_scale;
rho_f_dl = rho_f*v_scale^2/p_scale;
rho_p_dl = rho_p*v_scale^2/p_scale; 
rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
d_dl = d/z_scale;
tau0_dl = tau0/p_scale;
h_stop_dl = h_stop/z_scale;

beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);

u_w = y(1,1);
Q1 = y(2,:);
h = y(3,:);
u = u_w - Q1./h;
n = y(4,:);
m = y(5,:);

if full_model
    phi = y(6,:)./Q1;
    if var_rho
        rho_dl = rho_p_dl*phi+rho_f_dl*(1-phi);
    else
        rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
        
    end
    chi = (rho_f_dl+3*rho_dl)./(4*rho_dl);
    pb = y(7,:) + rho_dl.*g_dl*cosd(theta).*chi.*h;
    pe = pb-h;
else
    pb=h;
end
if ~var_rho
    rho_dl = rho_dl*ones(size(h));
end
p_tot_grad_dl = rho_dl*g_dl*cosd(theta);
Fr_vals = Fr.*u./sqrt(h);

h_min = roots([1,0,-u_w,-Q1(1)]);
h_static = Q1(1)/u_w;

[p_max,p_max_ind] = max(pb);
[p_min,p_min_ind] = min(pb);

[h_max,h_max_ind] = max(h);
[h_min,h_min_ind] = min(h);

if full_model
    [phi_max,phi_max_ind] = max(phi);
    [phi_min,phi_min_ind] = min(phi);
end

pp = p_tot_grad_dl.*h-pb;
D = -2/beta_dl./h.*(pb-h);
Iv = 3*eta_f_dl.*abs(u)./h./pp;

if full_model
    zeta = 3./(2*alpha_dl.*h) + P/4;
    tan_psi = phi - phi_c./(1+sqrt(Iv));
    R_w3 = -phi.*rho_f_dl/rho_dl.*D;
    if var_rho
        R_w4 = (-P/4+zeta).*D - 9/2/alpha_dl./h.*u.*(tan_psi);
    else
        R_w4 = (-P.*chi+zeta).*D - 9/2/alpha_dl./h.*u.*(tan_psi);
    end
end

Fr_equi = zeros(size(h));
Iv_equi = zeros(size(h));
ave_wave = 0;
for i=1:size(h,2)
    [Fr_equi(i),Iv_equi(i)] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h(i)*h0, tau0,0);
    if i ~= 1
        ave_wave = ave_wave + (rho_dl(i)*h(i)+rho_dl(i-1)*h(i-1))/2.0*(xi(i)-xi(i-1))/lambda;
    end
end

dhdxi = n;
n_coeff = 1-Q1.^2.*Fr^2./h.^3;
mu_val = pp./(p_tot_grad_dl.*h).*mu_Iv_fn(Iv)+tau0_dl*rho_f/rho./h;
force_bal = tand(theta)-sign(u).*mu_val;

mu_val_min = (rho-rho_f)/rho.*mu1_Iv+tau0_dl*rho_f/rho./h;
force_bal_max = tand(theta)-sign(u).*mu_val_min;

n_eq = (force_bal)./n_coeff;
n_diff = n_coeff.*n - force_bal;
dQdxi = -P.*D;

dndxi = 1./h.*n.^2 + h.^2/Fr^2./nu_dl./Q1.*(n_coeff.*n-force_bal);
dn_term1 = 1./(2.*h).*n.^2;
dn_term2 = h.^(3/2)/Fr^2./nu_dl./Q1.*n_coeff.*(n-n_eq);
dmdxi = h./lambda.*u;

if full_model
    dy6dxi = -R_w3;
    dphidxi = (dy6dxi-phi.*dQdxi)./Q1;
    dPdxi = rho_f_dl*(rho_p_dl-rho_f_dl)./rho_dl.^2.*dphidxi;
    dy7dxi = R_w4./(u-u_w);
    if var_rho
        dpbdxi = dy7dxi + rho_dl.*g_dl.*cosd(theta).*chi.*dhdxi + 3/4.*g_dl.*cosd(theta).*(rho_p_dl-rho_f_dl).*dphidxi;
    else
        dpbdxi = dy7dxi+ rho_dl.*g_dl.*cosd(theta_in).*chi.*dhdxi;
    end
    dDdxi = -2/beta_dl.*(dpbdxi./h-pb./h.^2.*dhdxi);
    dndxi_add = -h.*D./Q1.*dPdxi-h.*P./Q1.*dDdxi+P.*D./Q1.*n;
    dndxi = dndxi+dndxi_add;
    dpbdxi_scale = dpbdxi/(p_max-p_min);
    dhdxi_scale = n/(h_max-h_min);
end
%%
% n_out = 5000;
% outpoints = linspace(0,lambda,n_out);
% vec_save = interp1(xi,h0*h,outpoints);
% save("time_d_load_h.txt","vec_save","-ascii");
% vec_save = interp1(xi,h0*h.*u0.*u,outpoints);
% save("time_d_load_hu.txt","vec_save","-ascii");
% vec_save = interp1(xi,h0*h.*phi,outpoints);
% save("time_d_load_hphi.txt","vec_save","-ascii");
% vec_save = interp1(xi,h0*h.*crit_pb.*y(7,:),outpoints);
% save("time_d_load_pbh.txt","vec_save","-ascii");
%%
C = viridis(4);
% SetPaperSize(7.5,7.5)
hold on
%     plot(linspace(0.5,1),get_force_bal(linspace(0.5,1)))
plot(xi,pb, "DisplayName", "$\nu = 0.1\nu_n$","color",C(1,:))%
%     plot(xi,h,"--","DisplayName","$Q_1/u_w$","color","r")
%     plot(xi,ones(size(xi))*h_stop_dl,"--","DisplayName","$h_{stop}$","color","y")
%     plot(xi(xi<5),n_eq(xi<5), "DisplayName", "Waveform","color",C(2,:))
%     plot(Iv_equi,h, "DisplayName", "Equilibrium Curve","color",C(2,:))
%     plot(xi,dn_term2, "DisplayName", "$\frac{d y_7}{d \xi}$","color",C(2,:))
%     plot(xi,h_stop_dl*ones(size(xi)), "DisplayName", "Minimum $h$","color",C(3,:))
%     plot(xi,dpbdxi(xi>20), "DisplayName", "$\frac{d p_b}{d \xi}$","color",C(3,:))

xlabel("$\xi$")
ylabel("$h$")
legend("Location","best")
% title("$\theta="+num2str(theta)+"$, $\tau_0="+num2str(tau0)+"$, No $p_e$ case")
%     plot(xi,force_bal)
%     xlim([0,0.1])
%     exp_graph(gcf,"no_pe_tau0_20_u_stop_h_Iv.pdf")

function fb = get_force_bal(h_fb)
    u_fb = u_w - Q1(1)./h_fb;
    Iv_fb = crit_Iv.*u_fb./h_fb.^2;
    mu_fb = (rho-rho_f)/rho.*mu_Iv_fn(Iv_fb)+tau0_dl*rho_f/rho./h_fb;
    fb = tand(theta)-sign(u_fb).*mu_fb;
end
