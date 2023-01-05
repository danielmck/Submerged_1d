function nonlinear_syst
% Something else that does not work, trying to find solutions to the full
% system with the bvp method
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    theta = 10;
    
    alpha = 1e-5;
    rho_p = 2500;
    d = 1e-4;

    rho_f = 1000;
    eta_f = 0.0010016; % Pa s

%     rho_f = 1; % kg/m^3
%     eta_f = 1.18e-5; % Pa s
    
    rho = phi_c*rho_p + (1-phi_c)*rho_f;
    
    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);
    crit_phi = phi_c./(1+sqrt(crit_Iv));
    u_const = crit_Iv/eta_f/3*(rho_p-rho_f)*g*phi_c*cosd(theta);
    
    crit_Fr = get_critical_Fr(theta, rho_p, rho_f, d, eta_f, alpha);
    Fr = 0.6;
    h0 = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3);
    
    crit_u = u_const*h0^2;
    p_tot = rho*g*cosd(theta);
    crit_pb = rho_f*g*cosd(theta)*h0;
    
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

    P = (rho_dl-rho_f_dl)/rho_dl;
    
    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
    
    chi = (rho_f+3*rho)/(4*rho);
%     zeta = 3/(2*alpha_dl) +g_dl*cosd(theta)*rho_f_dl*P/4;
    
    k = 50;
    
    a=0;
    b=2*pi/k+a;
    N=100;
    
    [h_guess, u_guess] = roll_waves(b-a,theta,rho_f_dl,rho_p_dl,eta_f_dl,crit_Fr);
    guess_pts = size(h_guess,1);
    
    A_mat = make_A_mat(k,rho_p,rho_f,theta,eta_f,d,alpha,Fr,crit_Iv);
    A_eig = eigs(A_mat);
    [~, idx] = sort(imag(A_eig),'descend');
    A_eig = A_eig(idx);
    
    if (imag(A_eig(1))<0)
        error("Wave is not unstable, try a different value")
    else
        u_w = real(A_eig(1))/k;
        solInit1=bvpinit(linspace(a,b,N),@RollWaveGuess);
        solN1 = bvp4c(@full_system,@bcfcn,solInit1);

%         solN1.parameters;
        x = solN1.x;
        w1_sol = solN1.y(1,:);
        w2_sol = solN1.y(2,:);
        w3_sol = solN1.y(3,:);
        w4_sol = solN1.y(4,:);
        phi_sol = w3_sol./w1_sol;
        h_sol = w2_sol;
        u_sol = (w1_sol+u_w.*h_sol)./h_sol;
        pb_sol = w4_sol + rho_dl*g_dl*cosd(theta)*chi.*h_sol;
        plot(x,h_sol)
    end
    
    function dvecdx = full_system(x,y)
        w1 = y(1);
        w2 = y(2);
        w3 = y(3);
        w4 = y(4);
        
        h = w2;
%         h = zeros(1,size(w1,2));
%         for i = 1:size(w1,2)
%             cube_roots = roots([0.5*g_dl*cosd(theta) 0 (-w2(i)+u_w*w1(i)) w1(i)^2]);
%             real_roots = real(cube_roots(abs(imag(cube_roots))<1e-12));
%             
%             
%             if (size(real_roots,1) > 1)
% %                 "Multiple real roots, need to check we have the right one"
%                 u_roots = (w1+u_w.*real_roots)./real_roots;
%                 pb_roots = w4 + rho_dl*g_dl*cosd(theta)*chi.*real_roots;
%                 pos = (real_roots>0) & (u_roots>0) & (pb_roots>0);
%                 if sum(pos)>1
%                     error("More than 1 viable root")
%                 elseif (sum(pos)<1)
%                     error("No viable roots")
%                 else
%                     h(i) = real_roots(pos);
%                 end
%             else    
%                 h(i) = real_roots(1);
%             end
%         end
        u = (w1+u_w.*h)./h;
        phi = w3./w1;
        pb = w4 + rho_dl*g_dl*cosd(theta)*chi.*h;
        
%         pb = h;
        zeta = 3/(2*alpha_dl*h) + P/4;
        p_p = p_tot_grad_dl*h-pb;
        D = -2/beta_dl/h*(pb-rho_f_dl*g_dl*cosd(theta)*h);
        Iv = abs(2*eta_f_dl*u/h/p_p);
        R_w1 = P*D;
        R_w2 = g_dl*sind(theta)*h - sign(u)*mu_Iv_fn(Iv)/rho_dl*p_p;
        R_w3 = -phi*rho_f_dl/rho_dl*D;
        R_w4 = (-P*rho_dl*g_dl*cosd(theta)*chi+zeta)*D - 2*3/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));
        
        dw1dx = R_w1;
        dw2dx = (R_w2-(2*u-u_w)*P*D)./(g_dl*cosd(theta)*h-(u-u_w)^2);
        dw3dx = R_w3;
        dw4dx = R_w4/(u-u_w);
        dvecdx = [dw1dx dw2dx dw3dx dw4dx];
%         dvecdx = [dw1dx dw2dx];
    end

    function res = bcfcn(ya,yb)
        ua = (ya(1)+u_w.*ya(2))./ya(2);
        ub = (yb(1)+u_w.*yb(2))./yb(2);
        fa = ya(2).*ua.*(ua-u_w)+0.5*g_dl*cosd(theta)*ya(2)^2;
        fb = yb(2).*ub.*(ub-u_w)+0.5*g_dl*cosd(theta)*yb(2)^2;
        res = [ya(1)-yb(1) fa-fb ya(3)-yb(3) ya(4)-yb(4)];
%         res = [ya(1)-yb(1) fa-fb];
    end

    function v=evp1Guess1(x)
        x_fn = sin(k/2*x-pi/2); %(x/(b-a)-0.5)*2;
        h = 1+0.65*x_fn;
        u = crit_u_dl+ 0.4*x_fn;
        phi = crit_phi+0.02*x_fn;
        pb = 1+0.7*x_fn;
        w1 = h.*(u-u_w);
        w3 = phi.*h.*(u-u_w);
%         v=[h.*(u-u_w), h.*u.*(u-u_w)+0.5*g_dl*cosd(theta)*h^2, phi.*h.*(u-u_w), pb-rho_dl*g_dl*cosd(theta)*chi.*h];
        v=[w1, h, w3, pb-rho_dl*g_dl*cosd(theta)*chi.*h];
%         v=[w1, h];
    end

    function v=RollWaveGuess(x)
        pt_ind = 1+(x-a)*(guess_pts-1)/(b-a);
        pt_frac = pt_ind-floor(pt_ind);
        if (pt_ind<guess_pts)
            h = h_guess(floor(pt_ind))*(pt_frac) + h_guess(floor(pt_ind)+1)*(1-pt_frac);
            u = u_guess(floor(pt_ind))*(pt_frac) + u_guess(floor(pt_ind)+1)*(1-pt_frac);
        else
            h = h_guess(end);
            u = u_guess(end);
        end
        phi = crit_phi+0*h;
        Iv_val = (phi_c-phi).^2./phi.^2;
        pp = u*eta_f_dl*2./h./Iv_val;
        pb = p_tot_grad_dl*h-pp;
%         pb = h;
        w1 = h.*(u-u_w);
        w3 = phi.*h.*(u-u_w);
%         v=[h.*(u-u_w), h.*u.*(u-u_w)+0.5*g_dl*cosd(theta)*h^2, phi.*h.*(u-u_w), pb-rho_dl*g_dl*cosd(theta)*chi.*h];
        v=[w1, h, w3, pb-rho_dl*g_dl*cosd(theta)*chi.*h];
    end

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
    
    function dmudIv = dmudIv_fn(Iv)
        dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
    end
end