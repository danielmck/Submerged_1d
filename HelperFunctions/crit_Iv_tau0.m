 function [crit_h, Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr_eq, tau0, dl, rho_var, phi_param)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    
    if ~exist("dl","var")
        dl = false;
    end
    if ~exist("rho_var","var")
        rho_var = false;
    end
    if ~exist("phi_param","var")
        phi_param = rho_var;
    end
    
    reg_param = 1*10^7;
    count = 0;
    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    if ~rho_var
        rho = phi_c*rho_p + (1-phi_c)*rho_f;
    end
    if tau0>0
        max_tol = 1e-7;
        resid = 1;   
        del_Iv = 1e-8;
        max_Iv = 1e-2;
        min_Iv = 1e-8;
        Iv = newt_solve_crit_Iv(theta,rho_p,rho_f,rho_var,phi_param);
        while (abs(resid)>max_tol && max_Iv>1.1e-8)
            resid = force_bal_fn(Iv);
            if Iv > 1e-7
                force_deriv = force_bal_deriv(Iv);
            else
                force_deriv = (force_bal_fn(Iv+del_Iv)-resid)/del_Iv;
            end
            if resid < 0
                min_Iv = Iv;
            else
                max_Iv = Iv;
            end
            if ((Iv-resid/force_deriv<min_Iv) || (Iv-resid/force_deriv>min_Iv))
                Iv = (max_Iv+min_Iv)/2;
            else
                Iv = Iv-resid/force_deriv;
            end
            crit_h = (3*Fr_eq*eta_f/(Iv*(rho-rho_f)*sqrt(g*cosd(theta))))^(2/3);
            count=count+1;
        end
    else
        Iv = newt_solve_crit_Iv(theta,rho_p,rho_f,rho_var,phi_param);
        if rho_var
            phi = phi_c/(1+sqrt(Iv));
            rho = phi_param*(phi*rho_p + (1-phi)*rho_f)+(1-phi_param)*(phi_c*rho_p + (1-phi_c)*rho_f);
        end
        crit_h = (3*Fr_eq*eta_f/(Iv*(rho-rho_f)*sqrt(g*cosd(theta))))^(2/3);
    end
    
    function fb = force_bal_fn(Iv)
        if rho_var
            phi = phi_c/(1+sqrt(Iv));
            rho = phi_param*(phi*rho_p + (1-phi)*rho_f)+(1-phi_param)*(phi_c*rho_p + (1-phi_c)*rho_f);
        end
        h=(3*Fr_eq*eta_f/(Iv*(rho-rho_f)*sqrt(g*cosd(theta))))^(2/3);
        pp = (rho-rho_f)*g*cosd(theta)*(h*(1-dl)+dl/(rho_f*g*cosd(theta)));
        crit_mu = rho/(rho-rho_f)*tand(theta)-tau0/pp;
        fb = mu_Iv_fn(Iv)-crit_mu;
    end
    
    function h=get_h(Iv)
        h = (3*Fr_eq*eta_f/(Iv*(rho-rho_f)*sqrt(g*cosd(theta))))^(2/3);
    end

    function hd=get_h_deriv(Iv)
        hd = (3*Fr_eq*eta_f/((rho-rho_f)*sqrt(g*cosd(theta))))^(2/3)*(-2/3)*Iv^(-5/3);
    end

    function rho_deriv=get_rho_deriv(Iv)
        phi_deriv = -phi_c/2/(1+sqrt(Iv))^2/sqrt(Iv);
        rho_deriv = (rho_p-rho_f)*phi_deriv;
    end

    function deriv = force_bal_deriv(Iv)
        h=get_h(Iv);
        if rho_var
            phi = phi_c/(1+sqrt(Iv));
            rho=phi_param*(phi*rho_p + (1-phi)*rho_f)+(1-phi_param)*(phi_c*rho_p + (1-phi_c)*rho_f);
            rho_deriv = phi_param*get_rho_deriv(Iv);
        else
            rho_deriv = 0;
        end
        pp = (rho-rho_f)*g*cosd(theta)*(h*(1-dl)+dl/(rho_f*g*cosd(theta)));
        hd=get_h_deriv(Iv);
%         pp_d = (rho-rho_f)*g*cosd(theta)*hd*(1-dl);
        deriv = dmudIv_fn(Iv)-tau0/pp^2*(pp/h*hd+pp/(rho-rho_f)*rho_deriv)+rho_f/(rho-rho_f)^2*tand(theta)*rho_deriv;
    end
end