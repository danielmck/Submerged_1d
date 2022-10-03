function nonlin_ode15s
    h_crit = 0.8; % layer height (m)
    d=1.43e-5; % grain diameter (m)
    d_dl = d/h_crit;

    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    eta_f = 0.0010016; % Pa s
    g=9.81; % m/s^2
    rho_f = 1000; % kg/m^3
    rho_p = 2500; % kg/m^3
    theta = 15; % deg
    alpha = 1e-4; % 1/Pa

    Fr_min = get_critical_Fr(theta, rho_p, rho_f, d_dl, eta_f, alpha);
    rho = phi_c*rho_p+(1-phi_c)*rho_f;

    p_tot_grad = rho*g*cosd(theta);
    
    v_scale = sqrt(g.*h_crit);
    p_scale = rho_f.*g.*h_crit;
    t_scale = sqrt(h_crit./g);
    z_scale = h_crit;
    
    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    g_dl = g*t_scale/v_scale; 

    crit_u_dl = crit_u/v_scale;
    p_tot_grad_dl = p_tot/p_scale*z_scale;
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale; 
    rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
    d_dl = d/z_scale;

    P = (rho_dl-rho_f_dl)/rho_dl;
    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
    
    chi = (rho_f+3*rho)/(4*rho);
    P = (rho-rho_f)/rho;
    zeta = 3/(2*alpha_dl) +g_dl*cosd(theta)*rho_f_dl*P/4;

    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);
    h_crit_min = (Fr_min*sqrt(g*cosd(theta))/(crit_Iv*pb_grad))^(2/3);
    if (h_crit < h_crit_min)
        error("Flow not deep enough for travelling waves to occur, must be above "+num2str(h_crit_min))  
    else
        crit_u = crit_Iv.*h_crit^2.*pb_grad./eta_f./2;
%         Q1 = h_crit^(3/2).*sqrt(g.*cosd(theta));
        u_w = crit_u + sqrt(h_crit.*g.*cosd(theta));

%         determ = -4*u_w^3+27*Q1^2*crit_Iv*pb_grad./eta_f./2;
%         eq_heights = roots([crit_Iv*pb_grad./eta_f./2 0 -u_w Q1]);
%         h_min = max(eq_heights(eq_heights<(h_crit-1e-6) & eq_heights>0));
%         h_max = get_h_plus(h_min);
        hold on
        n_line = 10;
        h_low = 5e-2;
        h_high = 0.95*h_crit;
        p = h_crit*crit_u;
        fluxes = zeros(n_line,1);
%         SetPaperSize(12,10)
        for i = 1:n_line
            h1 = h_low+(i-1)/n_line*(h_high-h_low);
            h2 = get_h_plus(h1);
            [xi_vals,h_vals]=ode15s(@h_deriv_fn,linspace(0,500,10001),h1);
            u_vals = get_u(h_vals);
            n_point = sum(h_vals<h2);
            plot(xi_vals(h_vals<h2),h_vals(h_vals<h2),'r')
            flux = (h_vals(1).*u_vals(1)+h_vals(n_point).*u_vals(n_point) + 2*sum(h_vals(2:n_point-1).*u_vals(2:n_point-1)))./(2*n_point);
            fluxes(i) = flux;
        end
        xlabel('$\xi(m)$')
        ylabel('h(m)')
        title('Roll waves on a 15 degree slope')
%         PrintFig('Ive_15deg_roll_wave')
    end

    function h_plus = get_h_plus(h_minus)
        h_plus = (-h_minus + sqrt(h_minus.^2+8*Q1^2./(g.*h_minus.*cosd(theta))))/2;
    end
    
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end

    function dvecdx = full_system(x,y)
        w1 = y(1);
        w2 = y(2);
        w3 = y(3);
        w4 = y(4);
        
        persistent h_min;
        if isempty(h_min)
            h_min = 0;
        end
        
        h = zeros(1,size(w1,2));
        for i = 1:size(w1,2)
            cube_roots = roots([0.5*g_dl*cosd(theta) 0 (-w2(i)+u_w*w1(i)) w1(i)^2]);
            real_roots = real(cube_roots(abs(imag(cube_roots))<1e-12));
            
            
            if (size(real_roots,1) > 1)
%                 "Multiple real roots, need to check we have the right one"
                u_roots = (w1+u_w.*real_roots)./real_roots;
                pb_roots = w4 + rho_dl*g_dl*cosd(theta)*chi.*real_roots;
                pos = (real_roots>0) & (u_roots>0) & (pb_roots>0);
                if sum(pos)>1
                    warning("More than 1 viable root")
                    h(i) = min(real_roots(real_roots>h_min));
                elseif (sum(pos)<1)
                    error("No viable roots")
                else
                    h(i) = real_roots(pos);
                end
            else
                if (real_roots(1) > 0)
                    h(i) = real_roots(1);
                else
                    error("No viable roots")
                end
            end
        end
        u = (w1+u_w.*h)./h;
        phi = w3./w1;
        pb = w4 + rho_dl*g_dl*cosd(theta)*chi.*h;
        
        
        p_p = p_tot_grad_dl*h-pb;
        D = -2/beta_dl/h*(pb-rho_f_dl*g_dl*cosd(theta)*h);
        Iv = 2*eta_f_dl*u/h/p_p;
        R_w1 = P*D;
        R_w2 = g_dl*sind(theta)*h - mu_Iv_fn(Iv)/rho_dl*p_p;
        R_w3 = -phi*rho_f_dl/rho_dl*D;
        R_w4 = (-P*rho_dl*g_dl*cosd(theta)*chi+zeta)*D + 2*3/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));
        
        dw1dx = R_w1;
        dw2dx = R_w2;
        dw3dx = R_w3;
        dw4dx = R_w4/(u-u_w);
        dvecdx = [dw1dx dw2dx dw3dx dw4dx];
    end

    function u = get_u(h)
        u = (-Q1 + h.*u_w)./h;
    end

    function num_val = force_bal(h)
        u = get_u(h);
        Iv = 2.*u.*eta_f./(pb_grad.*h.^2);
        num_val = tand(theta)-mu_Iv_fn(Iv).*(rho-rho_f)/rho;
    end
end