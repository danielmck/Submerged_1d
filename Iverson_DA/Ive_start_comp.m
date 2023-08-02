function Ive_start_comp


phi_c=0.585; % Volume fraction
eta_f = 0.0010016; % Pa s
g=9.81; % m/s^2
rho_f = 1000; % kg/m^3
rho_p = 2500; % kg/m^3
theta = 5; % deg
theta0 = 13;
% fric_ang = 0.65;

% h0 = 4e-2; % layer height (m)
d=1e-4; % grain diameter (m)
Fr = 5;
[h0, crit_Iv] = crit_Iv_tau0(theta0, rho_p, rho_f, eta_f, Fr, 0);
phi0 = phi_c/(1+sqrt(crit_Iv));
rho_eq = phi0*rho_p+(1-phi0)*rho_f;
u0 = crit_Iv*h0*((rho_eq-rho_f)*g*cosd(theta)*h0)/3/eta_f;
pb0 = rho_f.*g.*cosd(theta).*h0;
mass0 = rho_eq*h0;

alpha = 1e-4;
kappa = ((1-phi_c).^3.*d^2)./(150*phi_c.^2);

buoyancy = -(rho_p-rho_f)*g*phi_c*cosd(theta);

v_scale = u0;
p_scale = pb0;
z_scale = h0;
t_scale = z_scale/v_scale;
% density_ratio = rho_p./rho_f;
% p_tot_grad_dl = p_tot/p_scale*z_scale;
rho_f_dl = rho_f*v_scale^2/p_scale;
rho_p_dl = rho_p*v_scale^2/p_scale; 

% rho = phi_c*density_ratio+1-phi_c;
g_dl = g*t_scale/v_scale; 
d_dl = d/z_scale;
eta_f_dl = eta_f/(p_scale*t_scale);
alpha_dl = alpha*p_scale;
kappa_dl = kappa/(z_scale)^2;
beta_dl = eta_f_dl/kappa_dl;
mass0_dl = mass0/z_scale*v_scale^2/p_scale;

run_Ive_da_sim()


% movefile(fname,'Results/');
% EOS_write_record(fname,N,h,d,reg_param,density_ratio,phi_c,theta,eta_f_dl,a_dl,phi_rcp,phi_rlp,t_step,2,shear_lim_dl);

    function dvecdt=Ive_depth_ave(t,vec)
        h = vec(1);
        phi = vec(2);
        u = vec(3);
        pb = vec(4);
 
        rho_dl = rho_p_dl*phi+rho_f_dl*(1-phi);
        P = (rho_dl-rho_f_dl)/rho_dl;
        pe = (pb-rho_f_dl*g_dl*cosd(theta)*h);
        pp = rho_dl*h*g_dl*cosd(theta)-pb;
        Iv = 3*u*eta_f_dl/(h*(pp));
        tan_psi = phi-phi_c/(1+sqrt(Iv));
        tau_zx = pp*mu_Iv_fn(Iv);
        D = -2/(beta_dl*h)*pe;
        zeta = 3/(2*alpha_dl*h) + P/4;
        
        dhdt = (rho_dl-1)/rho_dl*D;
%         dhdt=0;
        dphidt = -phi*D/h;
        dudt = g_dl*sind(theta)*h-tau_zx/rho_dl+u*P*D;
        dpbdt = zeta*D-9/2*u/(h*alpha_dl)*tan_psi;
%         dIvdt = dudt.*3.*eta_f_dl./pp+1./pp.*Iv.*(-3*u/(h*alpha_dl)*(tan_psi));
%         tan_psi_eq = -(sind(theta)-tand(theta0)/cosd(theta)).*3.*eta_f_dl./((rho-1).*cosd(theta0))/(-init_Iv^2./eta_f_dl/((h*alpha_dl)));
        
        dvecdt = [dhdt,dphidt,dudt,dpbdt]';
    end

    function success = run_Ive_da_sim()        
        % No initial pressure of phihat and initial values of u_p and u_f
        % defined above
        time_vals = linspace(0,30000,1500);
%         opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');

%         [~,vec]=ode15s(@Ive_depth_ave,time_vals,[1,depth_phi_orig,depth_u,pb_init],opts);

        n_ind = 20;
        phi_val = phi0; %linspace(phi0-0.03,phi0+0.005,n_ind);
        u_val = 1; %linspace(0.8,2,n_ind);
        p_val = 1; %linspace(0.8,1.2,n_ind);
        
        opts = odeset('AbsTol',1e-6,'RelTol',1e-6,'Events',@IntEnd);
        for phi = phi_val
            rho_in = rho_p_dl*phi + (1-phi)*rho_f_dl;
            h_in = mass0_dl/rho_in;
            for u = u_val
                for pb = p_val
%                     ind = n_ind^2*(i-1)+n_ind*(j-1)+k;
                    
                    [tvals,vec]=ode15s(@Ive_depth_ave,[0,10000],[h_in,phi,u,pb]);
                    hold on
                    plot(vec(:,2),vec(:,3)*v_scale) %,vec(:,4)
                end
            end
        end
        ylabel("$u$")
        xlabel("$\phi$")
        zlabel("$p_b$")
%         [tvals,vec]=ode15s(@Ive_depth_ave,[0,1000],[1,init_phi,init_u,cosd(theta0)]);
%         vec = horzcat(tvals,vec);
%         size(vec)
%         save(fname, 'vec','-ascii')
%         success=1;
    end

    function [position,isterminal,direction] = IntEnd(t,y)
        position = abs(y(4)-rho_f_dl*g_dl*cosd(theta)*y(1))>1e-2 || (t<10);
        isterminal = 1;
        direction = 0;
    end
        
end