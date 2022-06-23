function DA_altered_sim
% All parameters are defined dimensionally and then non-dimensionalised
N = 200; % number of discretisation points in z for each of tau, u
h0 = 4e-2; % layer height (m)
d=1.43e-5; % grain diameter (m)

mu1_Iv = 0.32; % mu(Iv) rheology parameters
mu2_Iv = 0.7;
Iv_0 = 0.005;

reg_param = 10^8;

phi_c=0.585; % Volume fraction
eta_f = 0.0010016; % Fluid Viscosity (Pa s)
g=9.81; % Gravity (m/s^2)
rho_f = 1000; % Fluid Density (kg/m^3)
rho_p = 2500; % Particle Density (kg/m^3)
theta = 5; % Slope Angle (deg)

alpha = 1e-4; % Particle Compressibility
kappa = ((1-phi_c).^3.*d^2)./(150*phi_c.^2); % Matrix Permeability

v_scale = sqrt(g.*h0); % Typical Scales
p_scale = rho_f.*g.*h0;
t_scale = sqrt(h0./g);
z_scale = h0;
density_ratio = rho_p./rho_f;
rho = phi_c*density_ratio+1-phi_c;

d_dl = d/z_scale; % Dimensionless Quantities
eta_f_dl = eta_f/(p_scale*t_scale);
alpha_dl = alpha*p_scale;
kappa_dl = kappa/(z_scale)^2;

%     t_step = 1; % Time between points when specifying time values
%     num_points = 5000; % Num time values
%     time_vals = (0:num_points)*t_step; % time values of interest
    
time_vals = [0,5000]; % Start and end time

fname = "Ive_da_theta_4_deep_9_2_start.txt";
run_Ive_da_theta_sim()
movefile(fname,'DA_Results/');
% EOS_write_record(fname,N,h,d,reg_param,density_ratio,phi_c,theta,eta_f_dl,a_dl,phi_rcp,phi_rlp,t_step,2,shear_lim_dl);

    % This function is the same as the full depth averaged sytem apart from
    % the fact that the material is always at the dilatant equilibrium.
    % This means that there is no evolution equation for the basal
    % pressure.
    function dvecdt=Ive_depth_ave_red(t,vec)
        h = vec(1);
        phi = vec(2);
        u = vec(3);
 
        rho = density_ratio*phi+(1-phi);
        pp = (phi/(phi-phi_c))^2*3*u*eta_f_dl/h;
        pe = h*(rho-1)*cosd(theta)-pp;
        Iv_base = 3*u*eta_f_dl/(h*(pp));
        tau_zx = pp*mu_Iv_fn(Iv_base)+(1-phi)*eta_f_dl*2*u/h;
        D = -2*kappa_dl/(eta_f_dl*h)*pe;
        
        dhdt = (rho-1)/rho*D;
%         dhdt=0;
        dphidt = -phi*D/h;
        dudt = sind(theta)*h-tau_zx/rho;
        dvecdt = [dhdt,dphidt,dudt]';
    end

    % The depth averaged model with a slope angle that varies linearly from
    % theta_start to theta_stop over the first theta_delta_t of the
    % evolution.
    function dvecdt=Ive_depth_ave_theta(t,vec)
        h = vec(1);
        phi = vec(2);
        u = vec(3);
        pb = vec(4);
 
        theta_start = 9.2;
        theta_stop = 4;
        theta_delta_t = 10;
        
        theta = max(theta_stop + (theta_start-theta_stop)*(theta_delta_t-t)/theta_delta_t,theta_stop);
        rho = density_ratio*phi+(1-phi);
        pe = (pb-cosd(theta)*h);
        pp = rho*h*cosd(theta)-pb;
        Iv = 3*u*eta_f_dl/(h*(pp));
        Iv_base = 3*u*eta_f_dl/(h*(pp));
        tan_psi = phi-phi_c/(1+sqrt(Iv));
        tau_zx = pp*mu_Iv_fn(Iv_base)+(1-phi)*eta_f_dl*2*u/h;
        D = -2*kappa_dl/(eta_f_dl*h)*pe;
        
        dhdt = (rho-1)/rho*D;
        dphidt = -phi*D/h;
        dudt = sind(theta)*h-tau_zx/rho;
        dpbdt = -3*kappa_dl/(alpha_dl*eta_f_dl*h^2)*pe+cosd(theta)*dhdt/4-3*u/(h*alpha_dl)*(tan_psi);
        dvecdt = [dhdt,dphidt,dudt,dpbdt]';
    end
    
    % Runs the reduced simulation
    function success = run_Ive_da_red_sim()
        cd Results
        da_init_data = load('DA_Results/Ive_da_9_2_deep_13_start.txt');
        init_vec = [1 da_init_data(end,2:4)];
        cd ../
        % Uses the final state of a full depth system but has no need for
        % the pressure value.

        opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');
        [tvals,vec]=ode15s(@Ive_depth_ave_red,time_vals,init_vec,opts);
        vec = horzcat(tvals,vec);
        size(vec)
        save(fname, 'vec','-ascii')
        success=1;
    end

    % Runs the simulation with variable theta
    function success = run_Ive_da_theta_sim()
        da_init_data = load('DA_Results/Ive_da_9_2_deep_13_start.txt');
        init_vec = [1 da_init_data(end-1,2:end)];
        % Uses the final state of a depth averaged system simulation as the
        % ICs.
        
        opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');
        [tvals,vec]=ode15s(@Ive_depth_ave_theta,time_vals,init_vec,opts);
        vec = horzcat(tvals,vec);
        size(vec)
        save(fname, 'vec','-ascii')
        success=1;
    end

    function beta_val=beta_fn()
        beta_val = 150*phi_c.^2./((1-phi_c).^3.*d_dl^2);
    end
   
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
end