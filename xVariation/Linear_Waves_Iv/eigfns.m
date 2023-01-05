function eigfns
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    theta = 15;
    
    alpha = 1e-4;
    rho_p = 2500;
    d = 1e-4;

    rho_f = 1000;
    eta_f = 0.0010016; % Pa s

    %     rho_f = 1; % kg/m^3
    %     eta_f = 1.18e-5; % Pa s
    
    h0 = 0.5;
 
    rho = phi_c*rho_p + (1-phi_c)*rho_f;
    ptot = rho*g*phi_c*cosd(theta)*h0;

    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);
    phi0 = phi_c./(1+sqrt(crit_Iv));
    u_const = crit_Iv/eta_f/3*(rho_p-rho_f)*g*phi_c*cosd(theta);
    u0 = u_const*h0^2;
    pb0 = rho_f*g*cosd(theta)*h0;
    
    chi = (rho_f+3*rho)/(4*rho);
    P = (rho-rho_f)/rho;
    zeta = 3/(2*alpha*h0) +g*cosd(theta)*rho_f*P/4;
    
    dIvdu = crit_Iv/u0;
    dIvdh = -2*crit_Iv/h0;
    dIvdp = crit_Iv/(ptot-pb0);
    
    dmudu = dmudIv_fn(crit_Iv).*dIvdu;
    dmudp = dmudIv_fn(crit_Iv).*dIvdp;
    dmudh = dmudIv_fn(crit_Iv).*dIvdh;
    
    dpsidIv = phi_c/2/(1+sqrt(crit_Iv))^2/sqrt(crit_Iv);
    
    beta = 150*phi_c.^2.*eta_f./((1-phi_c).^3.*d^2);
    
    a=0; b=1; N=100; lambda=1;
    solInit1=bvpinit(linspace(a,b,N),@evp1Guess1,lambda);
    solN1 = bvp4c(@lin_system,@bcfcn,solInit1);
    solN1.parameters
    x = solN1.x;
    yh = solN1.y(1,:);
    yu = solN1.y(2,:);
    yphi = solN1.y(3,:);
    ypb = solN1.y(4,:);
    plot(x,yu)
    
    function dvecdx = lin_system(x,y,lambda)
        h1 = y(1);
        u1 = y(2);
        phi1 = y(3);
        pb1 = y(4);
        
        D = -2/beta/h0*pb1;
        R_h = P*D;
        R_u = tand(theta)/(rho-rho_f)/h0*pb1 - P*g*cosd(theta)*(dmudu*u1+dmudp*pb1+dmudh*h1);
        R_phi = -phi0*(P + rho_f/rho)*D;
        R_pb = (-rho_f*g*cosd(theta)+zeta)*D + 2*3/alpha/h0*u0*(phi1 + dpsidIv*(dIvdu*u1+dIvdh*h1+dIvdp*pb1));
        
        dh1dx = (R_u - lambda*(h0*u1 - u0*h1)-2*u0*R_h)/(h0*g*cosd(theta)-u0^2);
        du1dx = (u0*(R_u - lambda*(h0*u1+u0*h1))-(g*cosd(theta)+u0^2)*(R_h - lambda * h1))/(h0*u0^2-h0^2*(g*cosd(theta)));
        dphi1dx = (R_phi - h0*lambda*phi1)/h0/u0;
        dpb1dx = (R_pb - g*cosd(theta)*(chi*rho-rho_f)*du1dx - lambda*pb1)/u0;
        dvecdx = [dh1dx du1dx dphi1dx dpb1dx];
    end

    function res = bcfcn(ya,yb,lambda)
        res = [ya(1)-yb(1) ya(2)-yb(2) ya(3)-yb(3) ya(4)-yb(4) ya(1)];
    end

    function v=evp1Guess1(x)
        v=[0.01*sin(2*pi/b*x) 0.1*sin(2*pi/b*x) 0.01*sin(2*pi/b*x) 0.1*sin(2*pi/b*x)]; 
    end

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
    
    function dmudIv = dmudIv_fn(Iv)
        dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
    end
end
