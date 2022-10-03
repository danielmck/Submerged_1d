function int_to_wave_no_pe    
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    theta = 25;
    
    alpha = 1e-5;
    rho_p = 2500;
    d = 1e-4;

%     rho_f = 1000;
%     eta_f = 0.0010016; % Pa s

    rho_f = 1; % kg/m^3 
    eta_f = 1.18e-5; % Pa s
    
    rho = phi_c*rho_p + (1-phi_c)*rho_f;
    
    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);
    crit_phi = phi_c./(1+sqrt(crit_Iv));
    u_const = crit_Iv/eta_f/2*(rho_p-rho_f)*g*phi_c*cosd(theta);
    
    crit_Fr = get_critical_Fr(theta, rho_p, rho_f, d, eta_f, alpha);
    Fr = 0.75;
    h0 = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3);
    
    crit_u = u_const*h0^2;
    p_tot = rho*g*cosd(theta);
    crit_pb = rho_f*g*cosd(theta)*h0;
    nu = eta_f/rho_f;
    
    z_scale = h0;
    v_scale = crit_u;
    p_scale = crit_pb;
    t_scale = z_scale/v_scale;

    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    g_dl = g*t_scale/v_scale; 

    crit_u_dl = crit_u/v_scale;
    p_tot_grad_dl = p_tot/p_scale*z_scale;
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale; 
    rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
    d_dl = d/z_scale;
    nu_dl = nu/z_scale/v_scale;

    P = (rho_dl-rho_f_dl)/rho_dl;
    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
    chi = (rho_f+3*rho)/(4*rho);
    
%     u_w = 1+sqrt(g_dl*cosd(theta));
    u_w = 1+1/Fr-1e-3;
    Q1 = 1-u_w;
    eq_heights = roots([1 0 -u_w (u_w-1)]);
    
%     init_h = 0.78;
%     init_u = u_w+(1-u_w)/init_h;
%     init_pb = init_h;
%     init_phi = crit_phi;
%     init_n = get_n_eq(init_h,init_u);
%     init_w1 = init_h*(init_u-u_w);
%     init_w3 = (1-u_w)/sqrt(init_h)*init_n;
% %     init_vals = [init_w1, init_h, init_w3];
%     init_vals = [init_w1, init_h, init_n];
    
%     f=figure;
%     width = 10;
%     height = 10;
%     set(f, 'PaperUnits', 'centimeters');
%     set(f, 'PaperSize', [width height]);
    hold on
    n_vals = 2;
    init_h_list = linspace(0.999,1.001,n_vals);
    init_n_list = linspace(-0.001,0.001,n_vals);
%     for i = 1:n_vals
%         init_n = init_n_list(i);
%         for j = 1:n_vals
%             init_h = init_h_list(j);
%             init_u = u_w+(1-u_w)/init_h;
%             init_w1 = init_h*(init_u-u_w);
%             init_vals = [init_w1, init_h, init_n];
%             [xi_vals_f,out_vals_f]=ode15s(@full_system,[50, 100],init_vals);
% %             if out_vals_f(end,2)<2
%             plot(out_vals_f(:,2),out_vals_f(:,3),'color','k')
% %             end
% %             xlim([0.75,1.5])
%             ylim([-2,2])
% %             [xi_vals_b,out_vals_b]=ode15s(@full_system,[50, 0],init_vals);
% %             plot(out_vals_b(:,2),out_vals_b(:,3),'color',colours(j))
%         end
%     end
    [xi_vals_f,out_vals_f]=ode15s(@full_system,[0, 100],[Q1,1.001,0]);
    plot(xi_vals_f,out_vals_f(:,2),'color','k')
    ylabel("$\frac{\partial h}{\partial \eta}$")
    xlabel("$h$")
%     legend("$h_0 = $"+num2str(init_h_list(1)),"$h_0 = $"+num2str(init_h_list(2)),"$h_0 = $"+num2str(init_h_list(3)),"$h_0 = $"+num2str(init_h_list(4)),"$h_0 = $"+num2str(init_h_list(5)))
    title("Forwards Integration")
%     exp_graph(f,"ViscousWave_h_n_forward_int.png")
    opts = odeset('Events',@WaveEndFn);
%     [xi_wave,y_wave,lambda,x_out,ie]=ode15s(@full_system,[0, 10],final_vals,opts);
%     [~,min_ind] = min(y_wave(:,2));
%     y_wave = vertcat(y_wave(min_ind:end,:),y_wave(1:min_ind-1,:));
%     xi_wave = vertcat(xi_wave(min_ind:end,:),xi_wave(1:min_ind-1,:));
%     xi_wave = mod(xi_wave-xi_wave(1),lambda);
%     h1 = y_wave(1,2);
    
%     [xi_no_diff, h_no_diff, u_no_diff, n_no_diff] = single_roll_wave(init_h,theta,rho_dl,rho_f_dl,g_dl,eta_f_dl,p_tot_grad_dl-1,crit_Iv);
%     h_wave = y_wave(:,2);
%     u_wave = (y_wave(:,1)+u_w.*h_wave)./h_wave;
%     n_wave = y_wave(:,3);
%     pb_vals = out_vals(:,5) + rho_dl*g_dl*cosd(theta)*chi.*h_vals;
%     D_vals = -2/beta_dl./h_vals.*(pb_vals-rho_f_dl*g_dl*cosd(theta).*h_vals);
%     n_vals = (-out_vals(:,3)./h_vals)./(u_vals-u_w);
%     plot(xi_vals,y_wave);
    hold on
%     plot(h_no_diff,n_no_diff);
    
    function [position,isterminal,direction] = WaveEndFn(t,y)
        persistent grad_inv;
        if (t==0)
            grad_inv = 0;
        end
        init_grad = 2*(final_vals(3)>0)-1;
        if (grad_inv == 0)
            if (init_grad*y(3)<0)
                grad_inv = 1;
            end
        end
        if (~grad_inv || t<1e-6)
            position = 1*sign(y(2)-final_vals(2)+y(3)*1e-10);
        else
            position = y(2)-final_vals(2);
        end
        % The value that we want to be zero
        isterminal = 1;  % Halt integration 
        direction = init_grad;   % The zero can be approached from either direction
    end

    function dvecdx = full_system(x,y)
        h = y(2);
        u = (y(1)+u_w.*h)./h;
%         u = u_w+(1-u_w)/h;
        phi = crit_phi;
        pb = h;
        
        p_p = p_tot_grad_dl*h-pb;
        
        Iv = abs(2*eta_f_dl*u/h/p_p);
        R_w1 = 0;
        R_w3 = g_dl*cosd(theta)*h*force_bal(h,u);
        
%         n = (-y(3)/h)/(u-u_w);
%         n=y(3)*sqrt(h)/(1-u_w);

        n = y(3);
        n_coeff = 1-Q1^2.*Fr^2/h^3;
        n_eq = force_bal(h,u)/n_coeff;
        
        dy1dx = R_w1;
%         dy2dx = n;
%         dy3dx = (R_w3 - (g_dl*cosd(theta)*h-(u-u_w)^2)*n)/nu_dl;
        dy2dx = n;
%         dy2dx = R_w3 - (g_dl*cosd(theta)*h-(u-u_w)^2)*n/nu_dl*rho_dl;
        dy3dx = 1/(2*h)*n^2 + h^(3/2)/Fr^2/nu_dl/Q1*n_coeff*(n-n_eq);
%         dy3dx = ((g_dl*h*cosd(theta)-(1-u_w)^2/h^2)-R_w3)/nu_dl;
%         dy3dx = -((1-u_w)^2/h^2*n-g_dl*h*cosd(theta)*n+R_w3)/nu_dl;
        dvecdx = [dy1dx, dy2dx, dy3dx]';
    end

%     function n_eq = get_n_eq(h,u)
%         n_eq = g_dl.*h.^3.*cosd(theta).*force_bal(h,u)./(h.^3.*-Q1.^2);
%     end

    function num_val = force_bal(h,u)
        Iv = 2.*abs(u).*eta_f_dl./((p_tot_grad_dl-1).*h.^2);
        num_val = tand(theta)-sign(u)*mu_Iv_fn(Iv).*(rho_dl-rho_f_dl)/rho_dl;
    end

    function num_val = force_bal_h(h)
        u = (-(1-u_w) + h.*u_w)./h;
        Iv = 2.*u.*eta_f_dl./((p_tot_grad_dl-1).*h.^2);
        num_val = tand(theta)-mu_Iv_fn(Iv).*(rho_dl-rho_f_dl)/rho_dl;
    end

    function dhdxi = h_deriv_fn(xi,h)
        if (abs(h-1)>1e-6)
            dhdxi = g_dl.*cosd(theta).*h.^3.*force_bal_h(h)./(g_dl.*h.^3.*cosd(theta)-(1-u_w).^2);
        else
            dhdxi = (g_dl.*cosd(theta).*(3*h.^2.*force_bal(h)+h.^3.*(force_bal(h+del_h)-force_bal(h))./del_h))./(g_dl.*3.*h.^2.*cosd(theta));
        end
    end
    
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
    
    function dmudIv = dmudIv_fn(Iv)
        dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
    end
end