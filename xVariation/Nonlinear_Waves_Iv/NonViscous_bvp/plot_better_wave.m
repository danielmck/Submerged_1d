filelist = ["no_vis_better_pres_h.txt"]; %long_low_pe
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
rho = rho_p*phi_c+rho_f*(1-phi_c);

chi = (rho_f+3*rho)/(4*rho);
P = (rho-rho_f)/rho;

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
p_tot_grad_dl = zeros(n_files,1);
rho_f_dl = zeros(n_files,1);
rho_p_dl = zeros(n_files,1);
rho_dl = zeros(n_files,1);
d_dl = zeros(n_files,1);
tau0_dl = zeros(n_files,1);
h_stop_dl = zeros(n_files,1);
h_min = zeros(n_files,1);
beta_dl = zeros(n_files,1);

Q1 = cell([n_files,1]);
h = cell([n_files,1]);
u = cell([n_files,1]);
m = cell([n_files,1]);
xi = cell([n_files,1]);
phi = cell([n_files,1]);
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
dy4dxi = cell([n_files,1]);
dy5dxi = cell([n_files,1]);
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

    if tau0(i) == 0
        crit_Iv(i) = newt_solve_crit_Iv(theta(i), rho_p, rho_f);
        pp_eq_grad(i) = (rho_p-rho_f)*g*phi_c*cosd(theta(i));
        u_const(i) = crit_Iv(i)/eta_f/3*pp_eq_grad(i);
        h0(i) = ((Fr(i)*sqrt(g*cosd(theta(i))))./u_const(i))^(2/3);  
    else
        [h0(i), crit_Iv(i)] = crit_Iv_tau0(theta(i), rho_p, rho_f, eta_f, Fr(i), tau0(i));
    end
    u_eq(i) = Fr(i)*sqrt(g*cosd(theta(i))*h0(i));
    phi_eq(i) = phi_c/(1+sqrt(crit_Iv(i)));
    p_tot(i) = rho*g*cosd(theta(i));
    crit_pb(i) = rho_f*g*cosd(theta(i))*h0(i);

    h_stop(i) = tau0(i)/(rho*g*cosd(theta(i)))/(tand(theta(i))-(rho-rho_f)/rho*mu1_Iv);


    z_scale(i) = h0(i);
    v_scale(i) = u_eq(i);
    p_scale(i) = crit_pb(i);
    t_scale(i) = z_scale(i)/v_scale(i);

    eta_f_dl(i) = eta_f/(p_scale(i)*t_scale(i));
    alpha_dl(i) = alpha(i)*p_scale(i);
    g_dl(i) = g*t_scale(i)/v_scale(i); 

    u_eq_dl(i) = u_eq(i)/v_scale(i);
    p_tot_grad_dl(i) = p_tot(i)/p_scale(i)*z_scale(i);
    rho_f_dl(i) = rho_f*v_scale(i)^2/p_scale(i);
    rho_p_dl(i) = rho_p*v_scale(i)^2/p_scale(i); 
    rho_dl(i) = rho_p_dl(i)*phi_c+rho_f_dl(i)*(1-phi_c);
    d_dl(i) = d(i)/z_scale(i);
    tau0_dl(i) = tau0(i)/p_scale(i);
    h_stop_dl(i) = h_stop(i)/z_scale(i);

    beta_dl(i) = 150*phi_c.^2.*eta_f_dl(i)./((1-phi_c).^3.*d_dl(i)^2);

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
        y5{i,1} = y_temp(5,:);
        pb{i,1} = y_temp(5,:) + rho_dl(i)*g_dl(i)*cosd(theta(i))*chi.*h{i,1};
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

    pp{i,1} = p_tot_grad_dl(i).*h{i,1}-pb{i,1};
    D{i,1} = -2/beta_dl(i)./h{i,1}.*(pb{i,1}-h{i,1});
    Iv{i,1} = abs(3*eta_f_dl(i).*abs(u{i,1})./h{i,1}./pp{i,1});

    if full_model
        zeta{i,1} = 3./(2*alpha_dl(i).*h{i,1}) + P/4;
        phi_Iv{i,1} = phi_c./(1+sqrt(Iv{i,1}));
        tan_psi{i,1} = phi{i,1} - phi_Iv{i,1};
        R_w3{i,1} = -phi{i,1}.*rho_f_dl(i)./rho_dl(i).*D{i,1};
        diffusion{i,1} = (-P.*chi+zeta{i,1}).*D{i,1};
        dilatancy{i,1} = 2*3/alpha_dl(i)./h{i,1}.*u{i,1}.*(tan_psi{i,1});
        R_w4{i,1} = diffusion{i,1} - dilatancy{i,1};
    end

    Fr_equi{i,1} = zeros(size(h{i,1}));
    Iv_equi{i,1} = zeros(size(h{i,1}));
    for i=1:size(h,2)
        [Fr_equi{i,1}(i),Iv_equi{i,1}(i)] = crit_Iv_tau0_h(theta(i), rho_p, rho_f, eta_f, h{i,1}(i)*h0(i), tau0(i),0);
    end

    n_coeff{i,1} = 1-Q1{i,1}.^2.*Fr(i)^2./h{i,1}.^3;
    % Iv = 3*eta_f_dl(i).*abs(u)./h./pp;
    mu_val{i,1} = pp{i,1}./(p_tot_grad_dl(i).*h{i,1}).*mu_Iv_fn(Iv{i,1})+tau0_dl(i)*rho_f/rho./h{i,1};
    force_bal{i,1} = tand(theta(i))-sign(u{i,1}).*mu_val{i,1};
    dhdxi{i,1} = force_bal{i,1}./n_coeff{i,1};

    mu_val_min{i,1} = P.*mu1_Iv+tau0_dl(i)*rho_f/rho./h{i,1};
    force_bal_max{i,1} = tand(theta(i))-sign(u{i,1}).*mu_val_min{i,1};

    dQdxi{i,1} = -P.*D{i,1};
    dmdxi{i,1} = h{i,1}./lambda(i).*u{i,1};

    if full_model
        dy4dxi{i,1} = -R_w3{i,1};
        dy5dxi{i,1} = R_w4{i,1}./(u{i,1}-u_w(i));
        dpbdxi{i,1} = dy5dxi{i,1} + rho_dl(i)*g_dl(i)*cosd(theta(i))*chi.*dhdxi{i,1};

        dpbdxi_scale{i,1} = dpbdxi{i,1}/(p_max(i)-p_min(i));
        dhdxi_scale{i,1} = dhdxi{i,1}/(h_max(i)-h_min(i));
    end
end

tan_psi_e = tan_psi{i,1}(end);
phi_e = phi{i,1}(end);
u_e = u{i,1}(end);
h_e = h{i,1}(end);
Iv_e = Iv{i,1}(end);
pp_e = pp{i,1}(end);

app_rate = phi_e/(1+sqrt(Iv_e))*sqrt(Iv_e^3)/(u_w-u_e)/eta_f_dl/alpha_dl;
xi_app = linspace(19.95,20,100);
tan_psi_app = tan_psi_e*exp(app_rate*(xi_app-lambda));

%%

% C = viridis(3);
hold on
% plot(xi,2*Q1./Fr(i).^2.*dQdxi./(3.*h.^3), "DisplayName", "Waveform","color",C(1,:))
for i=1:n_files
    plot(xi{i,1}(1:end),tan_psi{i,1}, "DisplayName", "Waveform"), %"color",C(2,:)
    plot(xi_app,tan_psi_app, "DisplayName", "Waveform")
end
xlim([19.95,20])