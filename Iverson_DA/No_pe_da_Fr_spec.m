function No_pe_da_Fr_spec
h0 = 0.5; % layer height (m)
d=5e-3; % grain diameter (m)
tlen = 100;

[phi_c,rho_f,rho_p,~,eta_f,g] = get_params_water();
theta = 5; % deg
Fr0 = 5;
u0 = Fr0*sqrt(g*cosd(theta)*h0);
change_t = 0;
% fric_ang = 0.65;
rho = rho_p*phi_c+(1-phi_c)*rho_f;

alpha = 1e-5;
kappa = ((1-phi_c).^3.*d^2)./(150*phi_c.^2);

buoyancy = -(rho_p-rho_f)*g*phi_c*cosd(theta);

init_u =  u0; %init_Iv/eta_f/3.*(init_rho-rho_f).*cosd(theta0).*g*h0^2;

v_scale = sqrt(g.*h0);
p_scale = rho_f.*g.*h0;
t_scale = sqrt(h0./g);
z_scale = h0;
density_ratio = rho_p./rho_f;
rho_dl = phi_c*density_ratio+1-phi_c;

d_dl = d/z_scale;
eta_f_dl = eta_f/(p_scale*t_scale);
alpha_dl = alpha*p_scale;
kappa_dl = kappa/(z_scale)^2;
beta_dl = eta_f_dl/kappa_dl;
% init_Iv = newt_solve_crit_Iv(theta0,rho_p,rho_f);
change_t_dl = change_t/t_scale;
init_u_dl = init_u/v_scale;
% flux_con = load('../Iverson_Closure/Results/flux_conditions.txt');
% for l=1:50
%     flux=flux_con(l,1);
%     theta_init=flux_con(l,2);
%     init_phi=flux_con(l,3);
%     unit = floor(flux);pp_approx_new
%     tenth = mod(l,10);
%     if (tenth == 0)
%         fname = "Ive_DA_4_5_deep_"+num2str(unit)+"_flux.txt";
%     else
%         fname = "Ive_DA_4_5_deep_"+num2str(unit)+"_"+num2str(tenth)+"_flux.txt";
%     end
%     run_Ive_da_sim()
%     movefile(fname,'Results/');
% end
fname = "No_pe_da_"+num2str(theta)+"deg_Fr_"+num2str(Fr0)+"_gravel.txt";
run_Ive_da_sim()
movefile(fname,'Results/');
% EOS_write_record(fname,N,h,d,reg_param,density_ratio,phi_c,theta,eta_f_dl,a_dl,phi_rcp,phi_rlp,t_step,2,shear_lim_dl);

    function dvecdt=Ive_depth_ave(t,vec)
        h = vec(1);
        u = vec(2);
        
        if (change_t_dl > 0)
            if t < change_t_dl
                theta_in = theta0+(theta-theta0)/change_t_dl*t;
            else
                theta_in = theta;
            end
        else
            theta_in = theta;
        end
 
        pp = (rho_dl-1)*cosd(theta)*h;
        Iv = 3*abs(u)*eta_f_dl/(h*pp);
        tau_zx = pp*sign(u)*mu_Iv_fn(Iv);
        
        dhdt = 0;
%         dhdt=0;
        dudt = sind(theta_in)*h-tau_zx/rho_dl;
        
        dvecdt = [dhdt,dudt]';
    end

    function success = run_Ive_da_sim()
        time_vals = linspace(0,tlen,1500);
%         opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');

%         [~,vec]=ode15s(@Ive_depth_ave,time_vals,[1,depth_phi_orig,depth_u,pb_init],opts);
        [tvals,vec]=ode15s(@Ive_depth_ave,time_vals,[1,init_u_dl]);
        vec = horzcat(tvals,vec);
        size(vec);
%         vec = vec.*[t_scale,z_scale,1,v_scale,p_scale];
        save(fname, 'vec','-ascii')
        success=1;
    end
end