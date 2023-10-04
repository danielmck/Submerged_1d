filelist = ["report_demo_wave_v2.txt"]; %long_low_pe"",no_pe_static_show.txt
n_files = size(filelist,2);

[phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();

chi = (rho_f+3*rho)/(4*rho);
P = (rho-rho_f)/rho;

theta = zeros(n_files,1);
Fr = zeros(n_files,1);
tau0 = zeros(n_files,1);
full_model=true;
d = zeros(n_files,1);
alpha = zeros(n_files,1);
u_w = zeros(n_files,1);
Q1 = zeros(n_files,1);
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

h = cell([n_files,1]);
u = cell([n_files,1]);
m = cell([n_files,1]);
xi = cell([n_files,1]);
phi = cell([n_files,1]);
pb = cell([n_files,1]);
pe = cell([n_files,1]);
Fr_vals = cell([n_files,1]);
h_static = cell([n_files,1]);
h_grad = cell([n_files,1]);
pp = cell([n_files,1]);
Iv = cell([n_files,1]);
phi_Iv = cell([n_files,1]);
tan_psi = cell([n_files,1]);
Fr_equi = cell([n_files,1]);
Iv_equi = cell([n_files,1]);
n_coeff = cell([n_files,1]);
mu_val = cell([n_files,1]);
force_bal = cell([n_files,1]);
dhdxi = cell([n_files,1]);
dmdxi = cell([n_files,1]);



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
    record = readtable('Results/wave_record.csv');

    in_table = strcmp(record.Name, filename);
    wave_type = record.wave_type(in_table);
    theta(i) = record.theta(in_table); 
    Fr(i) = record.Fr(in_table);
    tau0(i) = record.tau0(in_table);

% s_c = 1-rho/(rho-rho_f)*tand(theta)/mu1_Iv;

    [h0(i), crit_Iv(i)] = crit_Iv_tau0(theta(i), rho_p, rho_f, eta_f, Fr(i), tau0(i));
    u_eq(i) = Fr(i)*sqrt(g*cosd(theta(i))*h0(i));
    phi_eq(i) = phi_c/(1+sqrt(crit_Iv(i)));
    p_tot(i) = rho*g*cosd(theta(i));
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
    p_tot_grad_dl(i) = p_tot(i)/p_scale(i)*z_scale(i);
    rho_f_dl(i) = rho_f*v_scale(i)^2/p_scale(i);
    rho_p_dl(i) = rho_p*v_scale(i)^2/p_scale(i); 
    rho_dl(i) = rho_p_dl(i)*phi_c+rho_f_dl(i)*(1-phi_c);
    d_dl(i) = d(i)/z_scale(i);
    tau0_dl(i) = tau0(i)/p_scale(i);
%     h_stop_dl(i) = h_stop(i)/z_scale(i);

    beta_dl(i) = 150*phi_c.^2.*eta_f_dl(i)./((1-phi_c).^3.*d_dl(i)^2);

% u_w(i) = y(1,1);
% lambda(i) = y(2,1);
% h_crit_xi(i) = y(3,1);
    Q1(i) = y_temp(2,1);
    u_w(i) = y_temp(1,1);
    lambda(i) = y_temp(3,1);
    h{i,1} = y_temp(4,:);
    u{i,1} = u_w(i) - Q1(i)./h{i,1};
    m{i,1} = y_temp(5,:);
    
    if size(y_temp,1) == 6
        stat_len = y_temp(6,1);
    else
        stat_len = 0;
    end

    xi{i,1} = stat_len+xi_temp*lambda(i);
    
    if stat_len > 0
%         Q1 = horzcat(Q1(1),Q1);
%         u_w = horzcat(u_w(1),u_w);
        xi{i,1} = horzcat([0,xi{i,1}(1)*0.99],xi{i,1});
        h{i,1} = horzcat([Q1(i)/u_w(i),Q1(i)/u_w(i)],h{i,1});
        u{i,1} = horzcat([0,0],u{i,1});
        m{i,1} = horzcat([0,0],m{i,1});
    end

    pb{i,1} = h{i,1};

    Fr_vals = Fr(i).*u{i,1}./sqrt(h{i,1});

    h_min = roots([1,0,-u_w(i),-Q1(1)]);
    h_static = Q1(i)/u_w(i);

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
    Iv{i,1} = abs(3*eta_f_dl(i).*abs(u{i,1})./h{i,1}./pp{i,1});

    h_alt = linspace(h_static*h0(i)+1e-4,1.5*h0(i),100);
    Fr_equi{i,1} = zeros(size(h_min));
    Iv_equi{i,1} = zeros(size(h_min));
    for j=1:size(h_alt,2)
        [Fr_equi{i,1}(j),Iv_equi{i,1}(j)] = crit_Iv_tau0_h(theta(i), rho_p, rho_f, eta_f, h_alt(j), tau0(i),0);
    end

    n_coeff{i,1} = 1-Q1(i).^2.*Fr(i)^2./h{i,1}.^3;
    % Iv = 3*eta_f_dl(i).*abs(u)./h./pp;
    mu_val{i,1} = pp{i,1}./(p_tot_grad_dl(i).*h{i,1}).*mu_Iv_fn(Iv{i,1})+tau0_dl(i)*rho_f/rho./h{i,1};
    force_bal{i,1} = tand(theta(i))-sign(u{i,1}).*mu_val{i,1};
    dhdxi{i,1} = force_bal{i,1}./n_coeff{i,1};

    dmdxi{i,1} = h{i,1}./lambda(i).*u{i,1};
    C = viridis(3);
    SetPaperSize(8,8)
    hold on
% plot(xi,2*Q1./Fr(i).^2.*dQdxi./(3.*h.^3), "DisplayName", "Waveform","color",C(1,:))
    h_alt_dl = h_alt/h0(i);
    uval = u_w(i) - Q1(i)./h_alt_dl;
    Iv_in = crit_Iv(i).*uval./h_alt_dl.^2;
    mu = P.*mu_Iv_fn(Iv_in)+tau0_dl(i)*rho_f/rho./h_alt_dl;
    fb = tand(theta(i))-sign(uval).*mu;
%     plot(h_alt_dl,fb, "DisplayName", "$\tau_0 = "+num2str(tau0(i))+"Pa$", "color",C(i,:))
%     if (i==2)
%         plot(Q1(i)/u_w(i),fb(1), 'rx','HandleVisibility','on',"DisplayName","$\frac{Q}{u_w}$")
%     else
%         plot(Q1(i)/u_w(i),fb(1), 'rx','HandleVisibility','off',"DisplayName","$\frac{Q}{u_w}$")
%     end
end


%%

% C = viridis(3);
% SetPaperSize(8,8)
% hold on
plot(xi{i,1}*h0,h{i,1}*h0, "DisplayName", "Waveform","color",C(1,:))
for i=1:n_files
%     plot(h_alt,force_bal{i,1}, "DisplayName", "Waveform relation", "color",C(1,:))
%     plot(h_alt,max(u_eq.*(u_w(i) - Q1(i)./h_alt*h0),0), "DisplayName", "Waveform relation", "color",C(1,:))
%     plot(h_alt, Fr_equi{i,1}.*sqrt(h_alt*g*cosd(theta)), "DisplayName", "Equilibrium Curve", "color",C(2,:))
end
ax=gca;
ax.YAxis.Exponent = 0;
xlabel("Dimensionless Flow height") %$p_b$ ($Pa$)
ylabel("Flow height $h$ ($m$)")
% ylabel("Velocity $u$ ($ms^{-1}$)")
% ylabel("Dimensionless force balance")
xlabel("$\xi$ ($m$)")
% xlabel("Froude number")
% xlim([0,])
% ylim([0,5])
% legend("Location","best")
title("$\theta = "+num2str(theta(1))+"^{\circ}$, Base flow $Fr = "+num2str(Fr(1))+"$, $\tau_0 = "+num2str(tau0)+"Pa$"); %$
exp_graph(gcf,"two_eqn_master_h.pdf")