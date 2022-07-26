mu1_Iv = 0.32;
mu2_Iv = 0.7;
Iv_0 = 0.005;

reg_param = 1*10^7;

rho_p = 2500;

%     rho_f = 1;
%     eta_f = 1.18e-5;

rho_f = 1000;
eta_f = 0.0010016;

d = 1e-3;
alpha = 1e-4;

phi_c=0.585; % Volume fraction

g=9.81; % m/s^2
rho = phi_c*rho_p + (1-phi_c)*rho_f;

theta = [8.6,15];
n_cases = size(theta,2);

k = 100;
Fr = 0.2;
h = 1e-2;
A_mat = cell(2,1);
for i = 1:n_cases
    crit_Iv(i) = newt_solve_crit_Iv(theta(i), rho_p, rho_f);
    crit_phi(i) = phi_c./(1+sqrt(crit_Iv(i)));
    u_const(i) = crit_Iv(i)/eta_f/2*(rho_p-rho_f)*g*phi_c*cosd(theta(i));

%         h(i) = ((Fr*sqrt(g*cosd(theta(i))))./u_const(i))^(2/3); % layer height (m)
    h = h.*ones(1,2);

    p_p(i) = (rho_p-rho_f)*g*phi_c*cosd(theta(i))*h(i);
    crit_pb(i) = rho_f*g*cosd(theta(i))*h(i);
    crit_u(i) = crit_Iv(i)/eta_f/2*p_p(i)*h(i);
    Fr(i) = crit_u(i)/sqrt(g*cosd(theta(i))*h(i));

    v_scale(i) = crit_u(i);
    p_scale(i) = crit_pb(i);
    z_scale(i) = h(i);
    t_scale(i) = z_scale(i)/v_scale(i);

    eta_f_dl(i) = eta_f/(p_scale(i)*t_scale(i));
    alpha_dl(i) = alpha*p_scale(i);
    g_dl(i) = g*t_scale(i)/v_scale(i);
    d_dl(i) = d/z_scale(i);

    rho_f_dl(i) = rho_f*v_scale(i)^2/p_scale(i);
    rho_p_dl(i) = rho_p*v_scale(i)^2/p_scale(i);
    rho_dl(i) = rho_p_dl(i)*phi_c+rho_f_dl(i)*(1-phi_c);
    p_p_dl(i) = p_p/p_scale;

    chi_dl(i) = (rho_f_dl(i)+3*rho_dl(i))/(4*rho_dl(i));
    P(i) = (rho_dl(i)-rho_f_dl(i))/rho_dl(i);
    zeta(i) = 3/(2*alpha_dl(i)) + P(i)/4;
    beta_dl(i) = 150*phi_c.^2.*eta_f_dl(i)./((1-phi_c).^3.*d_dl(i)^2);

    dmudu(i) = dmudIv_fn(crit_Iv(i)).*crit_Iv(i);
    dmudp(i) = dmudIv_fn(crit_Iv(i)).*-crit_Iv(i)/p_p_dl(i);
    dmudh(i) = dmudIv_fn(crit_Iv(i)).*-2*crit_Iv(i);

    A_mat{i,1}(1,1) = k; %+2*1i*P(i)/beta_dl;
    A_mat{i,1}(1,2) = k;
    A_mat{i,1}(1,3) = -2*P(i)*1i/beta_dl(i);

    A_mat{i,1}(2,1) = g_dl(i)*cosd(theta(i))*k-P(i)*g_dl(i)*cosd(theta(i))*dmudh(i)*1i;
    A_mat{i,1}(2,2) = k - 1i*P(i)*g_dl(i)*cosd(theta(i))*dmudu(i); 
    A_mat{i,1}(2,3) = 1i*tand(theta(i))/(rho_dl(i)-rho_f_dl(i)) - 1i*P(i)*g_dl(i)*cosd(theta(i))*dmudp(i);

    A_mat{i,1}(3,1) =  3*sqrt(crit_Iv(i))*1i/alpha_dl(i)/(1+sqrt(crit_Iv(i)))^2*phi_c; %- 1i*zeta*2/beta_dl-2/beta_dl(i)*P(i)*1i;
    A_mat{i,1}(3,2) = (chi_dl(i)*rho_dl(i)-rho_f_dl(i))*g_dl(i)*cosd(theta(i))*k - 3/2*sqrt(crit_Iv(i))/(1+sqrt(crit_Iv(i)))^2*phi_c/alpha_dl(i)*1i;
    A_mat{i,1}(3,3) = k + 1i*2/beta_dl(i)*(P(i)-zeta(i)) - 3/2*sqrt(crit_Iv(i))/(1+sqrt(crit_Iv(i)))^2*phi_c/alpha_dl(i)*1i/P(i);
    A_mat{i,1}(3,4) = -3*1i*crit_phi(i)/alpha_dl(i);

    A_mat{i,1}(4,1) = 1i * -(P(i)-rho_f_dl(i)/rho_dl(i)).*2./beta_dl(i);
    A_mat{i,1}(4,3) = (P(i)+rho_f_dl(i)/rho_dl(i))*2/beta_dl(i)*1i;
    A_mat{i,1}(4,4) = k;

    A_eig(i,:) = eig(A_mat{i,1});
end

    function mu_val = mu_Iv_fn(Iv)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    phi_c = 0.585;

    reg_param = 1*10^7;
    mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
end

function dmudIv = dmudIv_fn(Iv)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    phi_c = 0.585;

    reg_param = 1*10^7;
    dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
end
