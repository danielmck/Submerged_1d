mu1_Iv = 0.32;
mu2_Iv = 0.7;
Iv_0 = 0.005;

reg_param = 1*10^7;

rho_p = 2500;

%     rho_f = 1;
%     eta_f = 1.18e-5;

rho_f = 1000;
eta_f = 0.0010016;

d = 1e-4;
alpha = [1e-6 1e-6];

phi_c=0.585; % Volume fraction

g=9.81; % m/s^2
rho = phi_c*rho_p + (1-phi_c)*rho_f;

theta = [18,25];
n_cases = size(theta,2);
phase = zeros(4,2);

k = [1 1];
Fr = [0.05 0.05];
% h = 1e-2;
A_mat = cell(2,1);
for i = 1:n_cases
    crit_Iv(i) = newt_solve_crit_Iv(theta(i), rho_p, rho_f);
    crit_phi(i) = phi_c./(1+sqrt(crit_Iv(i)));
    u_const(i) = crit_Iv(i)/eta_f/2*(rho_p-rho_f)*g*phi_c*cosd(theta(i));

    h(i) = ((Fr(i)*sqrt(g*cosd(theta(i))))./u_const(i))^(2/3); % layer height (m)
%     h = h.*ones(1,2);

    p_p(i) = (rho_p-rho_f)*g*phi_c*cosd(theta(i))*h(i);
    crit_pb(i) = rho_f*g*cosd(theta(i))*h(i);
    crit_u(i) = crit_Iv(i)/eta_f/2*p_p(i)*h(i);
    Fr(i) = crit_u(i)/sqrt(g*cosd(theta(i))*h(i));

    v_scale(i) = crit_u(i);
    p_scale(i) = crit_pb(i);
    z_scale(i) = h(i);
    t_scale(i) = z_scale(i)/v_scale(i);

    eta_f_dl(i) = eta_f/(p_scale(i)*t_scale(i));
    alpha_dl(i) = alpha(i)*p_scale(i);
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

    dIvdu(i) = crit_Iv(i);
    dIvdh(i) = -2*crit_Iv(i);
    dIvdp(i) = crit_Iv(i)/p_p_dl(i);

    dmudu(i) = dmudIv_fn(crit_Iv(i)).*dIvdu(i);
    dmudp(i) = dmudIv_fn(crit_Iv(i)).*dIvdp(i);
    dmudh(i) = dmudIv_fn(crit_Iv(i)).*dIvdh(i);

    dpsidIv(i) = phi_c/2/(1+sqrt(crit_Iv(i)))^2/sqrt(crit_Iv(i));
    
    A_mat{i,1}(1,1) = k(i); %+2*1i*P/beta_dl
    A_mat{i,1}(1,2) = k(i);
    A_mat{i,1}(1,3) = -2*P(i)*1i/beta_dl(i);
    
    A_mat{i,1}(2,1) = g_dl(i)*cosd(theta(i))*k(i)-P(i)*g_dl(i)*cosd(theta(i))*dmudh(i)*1i;
    A_mat{i,1}(2,2) = k(i) - 1i*P(i)*g_dl(i)*cosd(theta(i))*dmudu(i); 
    A_mat{i,1}(2,3) = 1i*tand(theta(i))/(rho_dl(i)-rho_f_dl(i)) - 1i*P(i)*g_dl(i)*cosd(theta(i))*dmudp(i);
    
    A_mat{i,1}(3,1) = -2*3*1i/alpha_dl(i)*dpsidIv(i)*dIvdh(i); %  - 1i*zeta(i)*2/beta_dl(i)-2/beta_dl(i)*P(i)*1i;
    A_mat{i,1}(3,2) = (chi_dl(i)*rho_dl(i)-rho_f_dl(i))*g_dl(i)*cosd(theta(i))*k(i) - 2*3*1i/alpha_dl(i)*dpsidIv(i)*dIvdu(i);
    A_mat{i,1}(3,3) = k(i) -2*3*1i/alpha_dl(i)*dpsidIv(i)*dIvdp(i)  + 1i*2/beta_dl(i)*(P(i)-zeta(i));
    A_mat{i,1}(3,4) = -2*3*1i*crit_phi(i)/alpha_dl(i);
    
    A_mat{i,1}(4,3) = (P(i)+rho_f_dl(i)/rho_dl(i))*2/beta_dl(i)*1i;
    A_mat{i,1}(4,4) = k(i);

    t = 1;
    ind = 4;
%     
    [A_vec A_val] = eig(A_mat{i,1});
    A_eig(i,:) = diag(A_val);
    [~, ind] = max(imag(A_eig(i,:)));
    x_val = linspace(0,2*pi/k(i),100);
    wave = real(exp(A_eig(i,4)*t).*A_vec(:,ind).*exp(1i.*k(i).*x_val));
    peak = max(wave,[],2);
    wave = wave./peak;
    
    for l = 1:4
        phase(l,i) = get_phase(wave(l,:),x_val,k(i));
    end
    plot(x_val,real(wave(:,:)))
    legend('Height','Velocity','Pressure','Volume Fraction')
end

function phase = get_phase(wave, x_val,k)
    found_rt = 0;
    if (size(wave,2) ~= size(x_val,2))
        "Vectors need to be same size"
    else
        for i = 1:(size(wave,2)-1)
            if (wave(i)*wave(i+1)<0)
                grad = (wave(i+1)-wave(i))/(x_val(i+1)-x_val(i));
                root = x_val(i) + wave(i)/grad;
                phase = mod(root*k + (grad > 0)*pi,2*pi);
                found_rt = 1;
                break
            end
        end
    end
    if (~found_rt)
        phase = NaN;
        "No root found"
    end
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
