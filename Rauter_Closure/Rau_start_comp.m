function Rau_start_comp

phi_c=0.585; % Volume fraction
eta_f = 0.0010016; % Pa s
g=9.81; % m/s^2
rho_f = 1000; % kg/m^3
rho_p = 2500; % kg/m^3
theta = 10; % deg
a = 130;
phi_rlp = 0.53;
phi_rcp = 0.63;
% theta0 = 13;
% fric_ang = 0.65;

% h0 = 4e-2; % layer height (m)
d=1e-4; % grain diameter (m)
Fr = 1.5;
[h0, crit_Iv] = crit_Iv_rauter(theta, rho_p, rho_f, eta_f, Fr, a, phi_rlp, phi_rcp, 0);
u0 = Fr*sqrt(g*cosd(theta)*h0);
pp0 = 3*eta_f*u0/h0/crit_Iv; %(rho_eq-rho_f).*g.*cosd(theta).*h0-a*(phi0-phi_rlp)/(phi_rcp-phi0);
rho_eq = pp0/h0/g/cosd(theta)+rho_f;
phi0 = (rho_eq-rho_f)/(rho_p-rho_f);
% rho_eq = phi0*rho_p+(1-phi0)*rho_f;
pb0 = rho_f.*g.*cosd(theta).*h0;
pp_shear0 = 3*eta_f*u0/(phi_c/phi0-1)^2/h0;
pp_contact0 = a*(phi0-phi_rlp)/(phi_rcp-phi0);

% u0 = crit_Iv*h0*pp_shear0/3/eta_f;
mass0 = rho_eq*h0;

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
kappa_dl = kappa/(z_scale)^2;
beta_dl = eta_f_dl/kappa_dl;
mass0_dl = mass0/z_scale*v_scale^2/p_scale;
a_dl = a/p_scale;

run_Ive_da_sim()


% movefile(fname,'Results/');
% EOS_write_record(fname,N,h,d,reg_param,density_ratio,phi_c,theta,eta_f_dl,a_dl,phi_rcp,phi_rlp,t_step,2,shear_lim_dl);

    function dvecdt=Ive_depth_ave(t,vec)
        h = vec(1);
        phi = vec(2);
        u = vec(3);
 
        rho_dl = rho_p_dl*phi+rho_f_dl*(1-phi);
        P = (rho_dl-rho_f_dl)/rho_dl;
        pp_contact = a_dl*(phi-phi_rlp)/(phi_rcp-phi);
        Iv_phi = (phi_c/phi-1)^2;
        pp_shear = 3*u*eta_f_dl/h/Iv_phi;
        pp = pp_contact+pp_shear;
        pb = rho_dl*h*g_dl*cosd(theta)-pp;
        pe = (pb-rho_f_dl*g_dl*cosd(theta)*h);
        
        Iv = 3*u*eta_f_dl/(h*(pp));
        tau_zx = pp*mu_Iv_fn(Iv);
        D = -2/(beta_dl*h)*pe;
        
        dhdt = (rho_dl-1)/rho_dl*D;
%         dhdt=0;
        dphidt = -phi*D/h;
        dudt = g_dl*sind(theta)*h-tau_zx/rho_dl+u*P*D;
%         dIvdt = dudt.*3.*eta_f_dl./pp+1./pp.*Iv.*(-3*u/(h*alpha_dl)*(tan_psi));
%         tan_psi_eq = -(sind(theta)-tand(theta0)/cosd(theta)).*3.*eta_f_dl./((rho-1).*cosd(theta0))/(-init_Iv^2./eta_f_dl/((h*alpha_dl)));
        
        dvecdt = [dhdt,dphidt,dudt]';
    end

    function success = run_Ive_da_sim()        
        % No initial pressure of phihat and initial values of u_p and u_f
        % defined above
        time_vals = linspace(0,30000,1500);
%         opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');

%         [~,vec]=ode15s(@Ive_depth_ave,time_vals,[1,depth_phi_orig,depth_u,pb_init],opts);

        phi_val = [phi0-0.02,phi0-0.01,phi0-0.005,phi0];
        u_val = [0.8,0.9,1,1.1];
        n_ind = 4;
        opts = odeset('AbsTol',1e-6,'RelTol',1e-6,'Events',@IntEnd);
        for i = 1:n_ind
            rho_in = rho_p_dl*phi_val(i) + (1-phi_val(i))*rho_f_dl;
            h_in = mass0_dl/rho_in;
            for j = 1:n_ind
%                 ind = n_ind^2*(i-1)+n_ind*(j-1);

                [tvals,vec]=ode15s(@Ive_depth_ave,[0,10000],[h_in,phi_val(i),u_val(j)],opts);
                hold on
                plot(vec(:,2),vec(:,3))
            end
        end
        ylabel("$u$")
        xlabel("$\phi$")
%         [tvals,vec]=ode15s(@Ive_depth_ave,[0,1000],[1,init_phi,init_u,cosd(theta0)]);
%         vec = horzcat(tvals,vec);
%         size(vec)
%         save(fname, 'vec','-ascii')
%         success=1;
    end

    function [position,isterminal,direction] = IntEnd(t,y)
        h = y(1);
        u = y(3);
        phi = y(2);
        pp_contact = a_dl*(phi-phi_rlp)/(phi_rcp-phi);
        rho_dl = rho_p_dl*phi+rho_f_dl*(1-phi);
        Iv_phi = (phi_c/phi-1)^2;
        pp_shear = 3*u*eta_f_dl/h/Iv_phi;
        pp = pp_contact+pp_shear;
        pb = rho_dl*h*g_dl*cosd(theta)-pp;
        
        position = abs(pb-rho_f_dl*g_dl*cosd(theta)*h)>1e-3 || (t<10);
        isterminal = 1;
        direction = 0;
    end
        
end