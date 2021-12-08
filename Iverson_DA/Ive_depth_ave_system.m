function Ive_depth_ave_system
N = 200; % number of discretisation points in z for each of tau, u
h0 = 4e-2; % layer height (m)
d=1.43e-5; % grain diameter (m)

mu1_Iv = 0.32;
mu2_Iv = 0.7;
Iv_0 = 0.005;

reg_param = 10^8;

phi_c=0.585; % Volume fraction
eta_f = 0.0010016; % Pa s
g=9.81; % m/s^2
rho_f = 1000; % kg/m^3
rho_p = 2500; % kg/m^3
theta = 5; % deg
fric_ang = 0.65;

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

t_step = 10;
% for l=1:70
%     q=fix(l/10);
%     r=mod(l,10);
%     fname = "Ive_da_"+num2str(q);
%     if (r==0)
%         fname = append(fname,"_deep_v3.txt");
%     else
%         fname = append(fname,"_"+num2str(r)+"_deep_v3.txt");
%     end
%     theta = 0.1*l;
%     run_Ive_da_sim()
%     movefile(fname,'DA_Results/');
% end
fname = "Ive_da_5_deep_25_start.txt";
run_Ive_da_sim()
movefile(fname,'DA_Results/');
% EOS_write_record(fname,N,h,d,reg_param,density_ratio,phi_c,theta,eta_f_dl,a_dl,phi_rcp,phi_rlp,t_step,2,shear_lim_dl);

    function dvecdt=Ive_depth_ave(t,vec)
        h = vec(1);
        phi = vec(2);
        u = vec(3);
        pb = vec(4);
 
        rho = density_ratio*phi+(1-phi);
        pe = (pb-cosd(theta)*h);
        pp = rho*h*cosd(theta)-pb;
        Iv = 3*u*eta_f_dl/(h*(pp));
        Iv_base = 3*u*eta_f_dl/(h*(pp));
        tan_psi = phi-phi_c/(1+sqrt(Iv));
        tau_zx = pp*mu_Iv_fn(Iv_base)+(1-phi)*eta_f_dl*2*u/h;
        D = -2*kappa_dl/(eta_f_dl*h)*pe;
        
        dhdt = (rho-1)/rho*D;
%         dhdt=0;
        dphidt = -phi*D/h;
        dudt = sind(theta)*h-tau_zx/rho;
        dpbdt = -3*kappa_dl/(alpha_dl*eta_f_dl*h^2)*pe+cosd(theta)*dhdt/4-3*u/(h*alpha_dl)*(tan_psi);
        dvecdt = [dhdt,dphidt,dudt,dpbdt]';
    end

    function success = run_Ive_da_sim()
%         cd Results
%         init_data = load('../Iverson_Closure/Results/dil_5deg_deep_comp_new_phic.txt');
% %         steady_state = init_data(5001,:);
% %         cd ../
%         vec = init_data(1,:);
%         depth_u = depth_average(vec(601:800)',200,1);
%         depth_phi_orig = phi_c+depth_average(vec(201:400)',200,1);
%         pb_init=cosd(theta)+vec(1);
%         Iv_init = 3*depth_u*eta_f_dl/((2.5*0.6+0.4)*cosd(theta)-pb_init);
%         depth_phi = phi_c/(1+sqrt(Iv_init));
        da_init_data = load('DA_Results/Ive_da_25_deep.txt');
        init_vec = da_init_data(end-1,2:end);
        
        % No initial pressure of phihat and initial values of u_p and u_f
        % defined above
        time_vals = [0,5000*t_step];
        opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');

%         [~,vec]=ode15s(@Ive_depth_ave,time_vals,[1,depth_phi_orig,depth_u,pb_init],opts);
        [tvals,vec]=ode15s(@Ive_depth_ave,time_vals,init_vec,opts);
        vec = horzcat(tvals,vec);
        size(vec)
        save(fname, 'vec','-ascii')
        success=1;
    end

    function beta_val=beta_fn(phihat)
        beta_val = 150*phi_c.^2./((1-phi_c).^3.*d_dl^2);
    end
   
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
end