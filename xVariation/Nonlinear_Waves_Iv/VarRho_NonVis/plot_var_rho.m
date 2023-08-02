filelist = ["var_rho_master_pres_h.txt"]; %long_low_pe
n_files = size(filelist,2);

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

theta = zeros(n_files,1);
Fr = zeros(n_files,1);
tau0 = zeros(n_files,1);
full_model=true;
d = zeros(n_files,1);
alpha = zeros(n_files,1);
u_w = zeros(n_files,1);
lambda = zeros(n_files,1);
h_crit_xi = zeros(n_files,1);
crit_Iv = zeros(n_files,1);
pp_eq_grad = zeros(n_files,1);
u_const = zeros(n_files,1);
h0 = zeros(n_files,1);
u_eq = zeros(n_files,1);
phi_eq = zeros(n_files,1);
p_tot = zeros(n_files,1);
crit_pb = zeros(n_files,1);
h_stop = zeros(n_files,1);
z_scale = zeros(n_files,1);
v_scale = zeros(n_files,1);
p_scale = zeros(n_files,1);
t_scale = zeros(n_files,1);
eta_f_dl = zeros(n_files,1);
alpha_dl = zeros(n_files,1);
g_dl = zeros(n_files,1);
u_eq_dl = zeros(n_files,1);
rho_f_dl = zeros(n_files,1);
rho_p_dl = zeros(n_files,1);
rho_eq_dl = zeros(n_files,1);
d_dl = zeros(n_files,1);
tau0_dl = zeros(n_files,1);
h_stop_dl = zeros(n_files,1);
h_min = zeros(n_files,1);
beta_dl = zeros(n_files,1);

p_tot_grad_dl = cell([n_files,1]);
Q1 = cell([n_files,1]);
h = cell([n_files,1]);
u = cell([n_files,1]);
m = cell([n_files,1]);
xi = cell([n_files,1]);
phi = cell([n_files,1]);
rho = cell([n_files,1]);
P = cell([n_files,1]);
chi = cell([n_files,1]);
y5 = cell([n_files,1]);
pb = cell([n_files,1]);
pe = cell([n_files,1]);
Fr_vals = cell([n_files,1]);
h_static = cell([n_files,1]);
h_grad = cell([n_files,1]);
pp = cell([n_files,1]);
D = cell([n_files,1]);
Iv = cell([n_files,1]);
zeta = cell([n_files,1]);
phi_Iv = cell([n_files,1]);
tan_psi = cell([n_files,1]);
R_w3 = cell([n_files,1]);
diffusion = cell([n_files,1]);
dilatancy = cell([n_files,1]);
R_w4 = cell([n_files,1]);
Fr_equi = cell([n_files,1]);
Iv_equi = cell([n_files,1]);
n_coeff = cell([n_files,1]);
mu_val = cell([n_files,1]);
force_bal = cell([n_files,1]);
dhdxi = cell([n_files,1]);
mu_val_min = cell([n_files,1]);
force_bal_max = cell([n_files,1]);
dQdxi = cell([n_files,1]);
dmdxi = cell([n_files,1]);
dudxi = cell([n_files,1]);
dy4dxi = cell([n_files,1]);
dy5dxi = cell([n_files,1]);
dhy5dxi = cell([n_files,1]);
dpbdxi = cell([n_files,1]);
dpbdxi_scale = cell([n_files,1]);
dhdxi_scale = cell([n_files,1]);



% mu_val_min
% force_bal_max
% [p_max,p_max_ind] = max(pb{i,1});
% [p_min,p_min_ind] = min(pb{i,1});
% 
% [h_max,h_max_ind] = max(h{i,1});
% [h_min,h_min_ind] = min(h{i,1});
% 
% if full_model
%     [phi_max,phi_max_ind] = max(phi);
%     [phi_min,phi_min_ind] = min(phi);
% end

for i=1:n_files
    filename = filelist(i);
    master_file = load("Results/"+filename);
    xi_temp = master_file(1,:);
    y_temp = master_file(2:end,:);
    record = readtable('Results/full_record.csv');

    in_table = strcmp(record.Name, filename);
    wave_type = record.wave_type(in_table);
    theta(i) = record.theta(in_table); 
    Fr(i) = record.Fr(in_table);
    tau0(i) = record.tau0(in_table);
    full_model=true;
    d(i) = record.d(in_table);
    alpha(i) = record.alpha(in_table);
    u_w(i) = record.u_w(in_table);
    lambda(i) = record.lambda(in_table);
    h_crit_xi(i) = record.crit_xi(in_table);
    
    

% s_c = 1-rho/(rho-rho_f)*tand(theta)/mu1_Iv;
    [h0(i), crit_Iv(i)] = crit_Iv_tau0(theta(i), rho_p, rho_f, eta_f, Fr(i), tau0(i),false,true);
    u_eq(i) = Fr(i)*sqrt(g*cosd(theta(i))*h0(i));
    phi_eq(i) = phi_c/(1+sqrt(crit_Iv(i)));
%     p_tot(i) = rho*g*cosd(theta(i));
    crit_pb(i) = rho_f*g*cosd(theta(i))*h0(i);
%     h_stop(i) = tau0(i)/(rho*g*cosd(theta(i)))/(tand(theta(i))-(rho-rho_f)/rho*mu1_Iv);


    z_scale(i) = h0(i);
    v_scale(i) = u_eq(i);
    p_scale(i) = crit_pb(i);
    t_scale(i) = z_scale(i)/v_scale(i);

    eta_f_dl(i) = eta_f/(p_scale(i)*t_scale(i));
    alpha_dl(i) = alpha(i)*p_scale(i);
    g_dl(i) = g*t_scale(i)/v_scale(i); 

    u_eq_dl(i) = u_eq(i)/v_scale(i);
    rho_f_dl(i) = rho_f*v_scale(i)^2/p_scale(i);
    rho_p_dl(i) = rho_p*v_scale(i)^2/p_scale(i); 
    rho_eq_dl(i) = rho_p_dl(i)*phi_eq(i)+rho_f_dl(i)*(1-phi_eq(i));
    d_dl(i) = d(i)/z_scale(i);
    tau0_dl(i) = tau0(i)/p_scale(i);
    h_stop_dl(i) = h_stop(i)/z_scale(i);

    beta_dl(i) = 150*phi_c.^2.*eta_f_dl(i)./((1-phi_c).^3.*d_dl(i)^2);
    flux_eq = (phi_eq(i)*rho_p_dl(i)+(1-phi_eq(i))*rho_f_dl(i));
    St(i) = rho_p*d^2*u_eq(i)/3/h0(i)/eta_f;
% u_w(i) = y(1,1);
% lambda(i) = y(2,1);
% h_crit_xi(i) = y(3,1);
    Q1{i,1} = y_temp(1,:);
    h{i,1} = y_temp(2,:);
    u{i,1} = u_w(i) - Q1{i,1}./h{i,1};
    m{i,1} = y_temp(3,:);

    xi{i,1} = horzcat(xi_temp(xi_temp<1)*h_crit_xi(i),h_crit_xi(i)+(xi_temp(xi_temp>=1)-1)*(lambda(i)-h_crit_xi(i)));

    if full_model
        phi{i,1} = y_temp(4,:)./Q1{i,1};
        rho{i,1} = rho_p_dl*phi{i,1}+rho_f_dl*(1-phi{i,1});
        chi{i,1} = (rho_f_dl+3*rho{i,1})./(4*rho{i,1});
        P{i,1} = (rho{i,1}-rho_f_dl)./rho{i,1};
        y5{i,1} = y_temp(5,:);
        pb{i,1} = y_temp(5,:) + rho{i,1}.*g_dl(i)*cosd(theta(i)).*chi{i,1}.*h{i,1};
        pe{i,1} = pb{i,1}-h{i,1};
    else
        pb{i,1} = h{i,1};
    end

    Fr_vals = Fr(i).*u{i,1}./sqrt(h{i,1});

    h_min = roots([1,0,-u_w(i),-Q1{i,1}(1)]);
    h_static = Q1{i,1}(1)/u_w(i);

    [p_max,p_max_ind] = max(pb{i,1});
    [p_min,p_min_ind] = min(pb{i,1});

    [h_max,h_max_ind] = max(h{i,1});
    [h_min,h_min_ind] = min(h{i,1});

    h_grad = (h{i,1}(2:end)-h{i,1}(1:end-1))./(xi{i,1}(2:end)-xi{i,1}(1:end-1));

    if full_model
        [phi_max,phi_max_ind] = max(phi{i,1});
        [phi_min,phi_min_ind] = min(phi{i,1});
    end
    p_tot_grad_dl{i,1} = rho{i,1}*g_dl(i)*cosd(theta(i));
    pp{i,1} = p_tot_grad_dl{i,1}.*h{i,1}-pb{i,1};
    D{i,1} = -2/beta_dl(i)./h{i,1}.*(pb{i,1}-h{i,1});
    Iv{i,1} = abs(3*eta_f_dl(i).*abs(u{i,1})./h{i,1}./pp{i,1});

    if full_model
        zeta{i,1} = 3./(2*alpha_dl(i).*h{i,1}) + P{i,1}/4;
        phi_Iv{i,1} = phi_c./(1+sqrt(Iv{i,1}));
        tan_psi{i,1} = phi{i,1} - phi_Iv{i,1};
        R_w3{i,1} = -phi{i,1}.*rho_f_dl(i)./rho{i,1}.*D{i,1};
        diffusion{i,1} = (-P{i,1}/4+zeta{i,1}).*D{i,1};
        dilatancy{i,1} = 9/2/alpha_dl(i)./h{i,1}.*u{i,1}.*(tan_psi{i,1});
        R_w4{i,1} = diffusion{i,1} - dilatancy{i,1};
    end

    Fr_equi{i,1} = zeros(size(h{i,1}));
    Iv_equi{i,1} = zeros(size(h{i,1}));
    for i=1:size(h,2)
        [Fr_equi{i,1}(i),Iv_equi{i,1}(i)] = crit_Iv_tau0_h(theta(i), rho_p, rho_f, eta_f, h{i,1}(i)*h0(i), tau0(i),0);
    end

    n_coeff{i,1} = 1-Q1{i,1}.^2.*Fr(i)^2./h{i,1}.^3;
    % Iv = 3*eta_f_dl(i).*abs(u)./h./pp;
    mu_val{i,1} = pp{i,1}./(p_tot_grad_dl{i,1}.*h{i,1}).*mu_Iv_fn(Iv{i,1})+tau0_dl(i)*rho_f_dl(i)./rho{i,1}./h{i,1};
    force_bal{i,1} = tand(theta(i))-sign(u{i,1}).*mu_val{i,1};
    dhdxi{i,1} = force_bal{i,1}./n_coeff{i,1};

    mu_val_min{i,1} = P{i,1}.*mu1_Iv+tau0_dl(i)*rho_f_dl./rho{i,1}./h{i,1};
    force_bal_max{i,1} = tand(theta(i))-sign(u{i,1}).*mu_val_min{i,1};

    dQdxi{i,1} = -P{i,1}.*D{i,1};
    dudxi{i,1} = -dQdxi{i,1}./h{i,1}+dhdxi{i,1}.*Q1{i,1}./h{i,1}.^2;
    dmdxi{i,1} = h{i,1}./lambda(i).*u{i,1};

    if full_model
        dy4dxi{i,1} = -R_w3{i,1};
        dy5dxi{i,1} = R_w4{i,1}./(u{i,1}-u_w(i));
        dhy5dxi{i,1} = R_w4{i,1}.*h{i,1}+(pb{i,1}-rho{i,1}.*g_dl(i)*cosd(theta(i)).*chi{i,1}.*h{i,1}).*dQdxi{i,1};
        dpbdxi{i,1} = dy5dxi{i,1} + rho{i,1}.*g_dl(i)*cosd(theta(i)).*chi{i,1}.*dhdxi{i,1};

        dpbdxi_scale{i,1} = dpbdxi{i,1}/(p_max(i)-p_min(i));
        dhdxi_scale{i,1} = dhdxi{i,1}/(h_max(i)-h_min(i));
    end
    
    wave_flux = 0;
    for k=2:size(h{i,1},2)
        wave_flux = wave_flux + (rho{i,1}(k-1)*u{i,1}(k-1)*h{i,1}(k-1)+rho{i,1}(k)*u{i,1}(k)*h{i,1}(k))/2*(xi{i,1}(k)-xi{i,1}(k-1))/lambda(i); 
    end
end



% Iv_coeff = 1/(u_w-u_e)/eta_f_dl/alpha_dl;
% app_rate = phi_e/(1+sqrt(Iv_e))*sqrt(Iv_e^3).*Iv_coeff;
% xi_app = linspace(lambda-0.10,lambda,100);
% tan_psi_eq = alpha_dl(i)*h{i,1}/9*2./u{i,1}.*(diffusion{i,1}+P{i,1}.*D{i,1});
% tan_psi_app = tan_psi_eq(end)+(tan_psi_e-tan_psi_eq(end))*exp(app_rate*(xi_app-lambda));

pb_approx = zeros(size(xi{i,1}));
pb2_approx = zeros(size(xi{i,1}));
pb2_resid_vec = zeros(size(xi{i,1}));
for j=1:size(xi{i,1},2)
    tan_psi_e = tan_psi{i,1}(j);
    rho_e = rho{i,1}(j);
    chi_e = chi{i,1}(j);
    phi_e = phi{i,1}(j);
    u_e = u{i,1}(j);
    h_e = h{i,1}(j);
    Iv_e = Iv{i,1}(j);
    Q1_e = Q1{i,1}(j);
    p_tot_e = p_tot_grad_dl{i,1}(j).*h_e;
    A = -phi_c(i);
    B = 3*eta_f_dl(i)*u{i,1}(j)./h{i,1}(j);
    C = p_tot_grad_dl{i,1}(j)*h{i,1}(j);
    D1 = alpha_dl(i)*h_e/9*2/u_e.*(-P{i,1}(j).*chi{i,1}(j)+zeta{i,1}(j)).*(-2/beta_dl(i)./h_e);
    E = h_e*D1+phi_e;

    x_sol = roots([D1, sqrt(B)*D1, (A-D1*C+E), -(D1*C*(sqrt(B))-E*sqrt(B))]);
    x_sol = x_sol(x_sol>0);
    p_sol = -x_sol.^2+C;
    pb_approx(1,j) = p_sol;
    
    D_pb1 = -2/beta_dl(i)./h_e.*(p_sol-h_e);
    
    Iv_pb1 = 3.*u_e.*eta_f_dl(i)/h_e./(p_tot_e-p_sol);
    mu_val_pb1 = (p_tot_e-p_sol)./(p_tot_e).*mu_Iv_fn(Iv_pb1)+tau0_dl(i)*rho_f_dl(i)./rho_e./h_e;
    force_bal_pb1 = tand(theta(i))-sign(u_e).*mu_val_pb1;
    dhdxi_pb1 = min(force_bal_pb1./n_coeff{i,1}(j),1);
    
    pb2_coeff = -u_e/Q1_e/(p_tot_e-p_sol).*(-9/2/alpha_dl(i).*u{i,1}(j)*(sqrt(Iv_pb1)./(p_tot_e-p_sol)./(1+sqrt(Iv_pb1)).^2)+(-P{i,1}(j)/4+zeta{i,1}(j)).*-2/beta_dl(i));
    pb2_resid = -u_e/Q1_e/(p_tot_e-p_sol).*(p_sol-rho_e*chi_e*g_dl(i).*cosd(theta(i))*h_e)*P{i,1}(j)*D_pb1+(Q1_e/h_e^2-u_e/h_e-u_e/(p_tot_e-p_sol)/h_e*g_dl(i).*cosd(theta(i))*(rho_e*h_e-rho_e*h_e*chi_e))*dhdxi_pb1;
    pb2_resid_vec(1,j) = pb2_resid;
    pb2 = -pb2_resid./pb2_coeff;
    pb2_approx(1,j) = pb2;
end

% pp_root = B/(phi_c(i)/phi_e-1)^2;
% x_ex = sqrt(pp_root);
% root_ex = x_ex^3*D+sqrt(B)*D*x_ex^2+(A-D*C+E)*x_ex-(D*C*(sqrt(B))-E*sqrt(B));
% plot_x = linspace(0.5,1,100);
% val_x = plot_x.^3*D+sqrt(B)*D*plot_x.^2+(A-D*C+E)*plot_x-(D*C*(sqrt(B))-E*sqrt(B));
% plot(plot_x,val_x)

h_min = h{i,1}(1);
h_pl = h{i,1}(end);
u_min = u{i,1}(1);
u_pl = u{i,1}(end);
Q1_shock = Q1{i,1}(end);

dppdxi = 1./Q1{i,1}.*dhdxi{i,1}.*(u_w-u{i,1}).*g_dl(i).*cosd(theta(i)).*(rho{i,1}.*h{i,1}-chi{i,1}.*rho{i,1}.*h{i,1})-(h{i,1}./4-pp{i,1})./Q1{i,1}.*dQdxi{i,1}+dhy5dxi{i,1}./Q1{i,1};
dIvdxi = 3*eta_f_dl(i)./h{i,1}./pp{i,1}.*(-dQdxi{i,1}./h{i,1}+(Q1{i,1}./h{i,1}.^2-u{i,1}./h{i,1}).*dhdxi{i,1}-u{i,1}./pp{i,1}.*dppdxi);

R4_term = -u{i,1}./pp{i,1}./Q1{i,1}.*dhy5dxi{i,1};
R1_term = (-1./h{i,1}+u{i,1}./pp{i,1}./Q1{i,1}.*(h{i,1}/4-pp{i,1})).*dQdxi{i,1};
dh_term = (Q1{i,1}./h{i,1}.^2-u{i,1}./h{i,1}-u{i,1}./pp{i,1}./Q1{i,1}.*(u_w-u{i,1}).*g_dl(i).*cosd(theta(i)).*(rho{i,1}.*h{i,1}-chi{i,1}.*rho{i,1}.*h{i,1})).*dhdxi{i,1};
Iv_deriv = 3*eta_f_dl(i)./h{i,1}./pp{i,1}.*(dh_term+R1_term+R4_term);

dil_cont = u{i,1}./pp{i,1}./Q1{i,1}.*h{i,1}.*dilatancy{i,1};
%%
SetPaperSize(8,8)
C = viridis(3);
hold on
% plot(xi,2*Q1./Fr(i).^2.*dQdxi./(3.*h.^3), "DisplayName", "Waveform","color",C(1,:))
for i=1:n_files
%     plot(xi{i,1}(1:end),dh_term+R1_term+R4_term, "DisplayName", "Waveform")
    plot(xi{i,1},rho{i,1}.*u{i,1}.*h{i,1}, "color",C(1,:), "DisplayName", "Waveform")
    plot(xi{i,1},wave_flux*ones(size(xi{i,1})),"--", "color",C(1,:), "DisplayName", "Waveform average")
    plot(xi{i,1},rho_eq_dl(i)*ones(size(xi{i,1})),"--", "color",C(2,:), "DisplayName", "Uniform equilibrium")
%     plot(xi{i,1}(1:end-1),diff(Iv{i,1})./diff(xi{i,1}), "DisplayName", "Waveform")
%     plot(xi{i,1}(1:end),pb{i,1}-pb_approx-pb2_approx, "DisplayName", "Waveform") %"color",C(2,:)
%     plot(xi{i,1}(1:end),tan_psi_eq, "DisplayName", "Waveform")
%     plot(xi_app,tan_psi_app, "DisplayName", "Waveform")
end
% ylim([-0.001,0.001])
% xlim([lambda-0.5,lambda])
xlabel("$\xi$")
ylabel("Flux")
legend("Location","best")
% exp_graph(gcf,"flux_wave_uniform_comp.pdf")