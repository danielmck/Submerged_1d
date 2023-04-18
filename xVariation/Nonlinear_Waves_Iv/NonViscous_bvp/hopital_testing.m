master_name = "no_vis_better_master.txt";
master_file = load("Results/"+master_name);
master_xi = master_file(1,:);
master_y = master_file(2:end,:);
record = readtable('Results/wave_record.csv');

in_table = strcmp(record.Name, master_name);
wave_type = record.wave_type(in_table);
theta = record.theta(in_table);
Fr = record.Fr(in_table);
tau0 = record.tau0(in_table);
alpha = record.alpha(in_table);
d = record.d(in_table);

pres_h=1;

[phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
rho = rho_p*phi_c+rho_f*(1-phi_c);
P = (rho-rho_f)/rho;
chi = (rho_f+3*rho)/(4*rho);

[h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0);
u_eq = Fr*sqrt(g*cosd(theta)*h0);
phi_eq = phi_c/(1+sqrt(crit_Iv));

crit_pb = rho_f*g*cosd(theta)*h0;

p_tot = rho*g*cosd(theta);

z_scale = h0;
v_scale = u_eq;
p_scale = crit_pb;
t_scale = z_scale/v_scale;

tau0_dl = tau0/p_scale;
eta_f_dl = eta_f/(p_scale*t_scale);
alpha_dl = alpha*p_scale;
g_dl = g*t_scale/v_scale; 

u_eq_dl = u_eq/v_scale;

rho_f_dl = rho_f*v_scale^2/p_scale;
rho_p_dl = rho_p*v_scale^2/p_scale;
d_dl = d/z_scale;

p_tot_grad_dl = p_tot/p_scale*z_scale;
rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);

extr = [diff(master_xi)>0, true];
master_y = master_y(:,extr);
master_xi = master_xi(:,extr);
crit_ind = sum(master_xi<1);
crit_ratio = master_y(3,1)/master_y(2,1);
xi = horzcat(master_xi(1:crit_ind)*crit_ratio,crit_ratio+(1-crit_ratio)*(master_xi(crit_ind+1:end)-1));
y_in = vertcat(master_y(4:8,:),master_y(1:2,:));

u_w = y_in(6,1);
lambda = y_in(7,1);

Q1 = y_in(1,:);
h = y_in(2,:);
u = (-Q1 + h.*u_w)./h;
m = y_in(3,:);
phi = y_in(4,:)./Q1;
pb = y_in(5,:) + g_dl*cosd(theta)*rho_dl*chi.*h;

zeta = 3./(2*alpha_dl*h) + P/4;
p_p = p_tot_grad_dl.*h-pb;
Iv = 3*eta_f_dl.*abs(u)./h./p_p;
D = -2/beta_dl./h.*(pb-h);

denom = (h.^3/Fr^2-Q1.^2);
fric_coeff = p_p./(p_tot_grad_dl*h).*mu_Iv_fn(Iv);
fb_val = tand(theta)-fric_coeff-tau0_dl.*rho_f./rho./h;
numer  = 1/Fr^2.*h.^3.*(fb_val+P.*D.*h);
dhdxi = numer./denom;

R_w3 = -phi.*rho_f_dl./rho_dl.*D;
R_w4 = (-P.*chi+zeta).*D - 2*3./alpha_dl./h.*u.*(phi - phi_c./(1+sqrt(Iv)));

dQdxi = -P.*D;
dmdxi = h./(lambda).*u.^(1-pres_h);
dy6dxi = -R_w3;
dy7dxi = R_w4./(u-u_w);

denom_grad_a = (3*h.^2./Fr^2.*dhdxi-2.*Q1.*dQdxi)*lambda;
denom_grad_c = (denom(2:end)-denom(1:end-1))./diff(xi);

A_pp_coeff = -mu_Iv_fn(Iv)./(p_tot_grad_dl.*h)+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./p_p;
A_pb_coeff = -P.*2./beta_dl-A_pp_coeff;
A_h_coeff = p_p./(p_tot_grad_dl.*h.^2).*mu_Iv_fn(Iv)+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*3*eta_f_dl.*abs(u)./h.^2./p_p+tau0_dl.*rho_f./rho./h.^2-P.*2./beta_dl+A_pp_coeff.*rho_dl.*g_dl.*cosd(theta);

full_h_coeff = (h.^3./Fr^2.*(g_dl*cosd(theta)*rho_dl.*chi.*A_pb_coeff+A_h_coeff-p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./u.*Q1./h.^2)+3.*h.^2./Fr^2.*(fb_val+P.*D.*h));
numer_h_term = full_h_coeff.*dhdxi;
numer_other_term = h.^3./Fr^2.*(dy7dxi.*A_pb_coeff+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./u.*dQdxi./h.^2);
numer_grad_a = (numer_h_term+numer_other_term)*lambda;
numer_grad_c = (numer(2:end)-numer(1:end-1))./diff(xi);

dhdxi_alter = dhdxi;

for i = 1:size(xi,2)
    if abs(denom(i))<1e-2
        grad_roots = roots([3*h(i).^2./Fr^2, -full_h_coeff(i)-2.*Q1(i).*dQdxi(i), -numer_other_term(i)]);
        if imag(grad_roots(1))==0
            h_grad = min(grad_roots(grad_roots>0));
            denom_grad_a(i) = (3*h(i).^2./Fr^2.*h_grad-2.*Q1(i).*dQdxi(i))*lambda;
            numer_grad_a(i) = (full_h_coeff(i)*h_grad+numer_other_term(i))*lambda;
        else
            h_grad = -1;
        end
        dhdxi_alter(i) = h_grad;
    end
end



hold on
plot(xi,denom)
% plot(xi,dhdxi)
% ylim([-1,2])