function u_sinusoid
h0 = 1e-1; % layer height (m)
d=1e-4; % grain diameter (m)

phi_c=0.585; % Volume fraction
eta_f = 0.0010016; % Pa s
g=9.81; % m/s^2
rho_f = 1000; % kg/m^3
rho_p = 2500; % kg/m^3
theta = 0; % deg
theta0 = 10; % deg

% fric_ang = 0.65;

alpha = 1e-4;
kappa = ((1-phi_c).^3.*d^2)./(150*phi_c.^2);

buoyancy = -(rho_p-rho_f)*g*phi_c*cosd(theta);

v_scale = sqrt(g.*h0);
p_scale = rho_f.*g.*h0;
t_scale = sqrt(h0./g);
z_scale = h0;
density_ratio = rho_p./rho_f;
rho = phi_c*density_ratio+1-phi_c;

d_dl = d/z_scale;
eta_f_dl = eta_f/(p_scale*t_scale);
alpha_dl = alpha*p_scale;
kappa_dl = kappa/(z_scale)^2;
beta_dl = eta_f_dl/kappa_dl;
init_Iv = newt_solve_crit_Iv(theta0,rho_p,rho_f,false,true);
init_phi = phi_c/(1+sqrt(init_Iv));
init_rho = init_phi*density_ratio+1-init_phi;
init_u = init_Iv/eta_f_dl/3.*(rho-1).*cosd(theta0);

%     movefile(fname,'Results/');
% end
phase = 0.01;
u_peak = 1;
run_Ive_da_sim()
% EOS_write_record(fname,N,h,d,reg_param,density_ratio,phi_c,theta,eta_f_dl,a_dl,phi_rcp,phi_rlp,t_step,2,shear_lim_dl);

    function dvecdt=Ive_depth_ave(t,vec)
        h = vec(1);
        phi = vec(2)/h;
        pb = vec(3);
        u = u_peak*(sin(t*phase)+sin(2*t*phase));
 
        rho = density_ratio*phi+(1-phi);
        P = (rho-1)./rho;
        pe = (pb-cosd(theta)*h);
        pp = rho*h*cosd(theta)-pb;
        Iv = 3*abs(u)*eta_f_dl/(h*(pp));
        tan_psi = phi-phi_c/(1+sqrt(Iv));
        D = -2/(beta_dl*h)*pe;
        zeta = 3/(2*alpha_dl*h) + cosd(theta)*P/4;
        dhdt = P*D;
%         dhdt=0;
        dphidt = -phi*D/rho;
        dpbdt = zeta*D-9/2*abs(u)/(h*alpha_dl)*(tan_psi);
%         dIvdt = dudt.*3.*eta_f_dl./pp+1./pp.*Iv.*(-3*u/(h*alpha_dl)*(tan_psi));
%         tan_psi_eq = -(sind(theta)-tand(theta0)/cosd(theta)).*3.*eta_f_dl./((rho-1).*cosd(theta0))/(-init_Iv^2./eta_f_dl/((h*alpha_dl)));
        
        dvecdt = [dhdt,dphidt,dpbdt]';
    end

    function success = run_Ive_da_sim()        
        % No initial pressure of phihat and initial values of u_p and u_f
        % defined above
        time_vals = [0,5000];
%         opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');

%         [~,vec]=ode15s(@Ive_depth_ave,time_vals,[1,depth_phi_orig,depth_u,pb_init],opts);
        [tvals,vec]=ode15s(@Ive_depth_ave,time_vals,[1,0.575,1]);
        vec = horzcat(tvals,vec(:,1:2),u_peak*(sin(tvals*phase)+sin(2*tvals*phase)),vec(:,3));
        fname = "sin_u_evo_play.txt";
        save("Results/"+fname, 'vec','-ascii')
        success=1;
    end    
end