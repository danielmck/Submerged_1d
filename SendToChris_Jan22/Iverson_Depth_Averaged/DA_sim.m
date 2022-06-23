function DA_sim
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
theta = 13; % Slope Angle (deg)

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

fname = "Ive_da_5_deep.txt"; % file name
run_Ive_da_sim()
movefile(fname,'Results/');
Ive_da_write_record(fname,h0,d,reg_param,density_ratio,alpha_dl,phi_c,theta,eta_f_dl);

    % Depth sveraged system
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
        % Iv and Iv_base were for when we assumed that Iv was not constant
        % across the domain, I have left them here in case we try and alter
        % the system to allow for that again.
        
        tan_psi = phi-phi_c/(1+sqrt(Iv)); % Dilatancy is const across domain
        tau_zx = pp*mu_Iv_fn(Iv_base)+(1-phi)*eta_f_dl*2*u/h; % Basal friction
        D = -2*kappa_dl/(eta_f_dl*h)*pe; %Rate of dilation
        
        dhdt = (rho-1)/rho*D;
        dphidt = -phi*D/h;
        dudt = sind(theta)*h-tau_zx/rho;
        dpbdt = -3*kappa_dl/(alpha_dl*eta_f_dl*h^2)*pe+cosd(theta)*dhdt/4-3*u/(h*alpha_dl)*(tan_psi);
        dvecdt = [dhdt,dphidt,dudt,dpbdt]';
    end

    function success = run_Ive_da_sim()

%         init_data = load('../Iverson_1D/Results/dil_5deg_deep_comp_new_phic.txt');
%         vec = init_data(1,:);
%         depth_u = depth_average(vec(601:800)',200,1);
%         depth_phi_orig = phi_c+depth_average(vec(201:400)',200,1);
%         pb_init=cosd(theta)+vec(1);
%         init_vec = [1,depth_phi_orig,depth_u,pb_init];
        % Depth averages ICs of 1D sim to be used as ICs here

        cd Results
        da_init_data = load('Ive_da_13_deep.txt');
        init_vec = [1 da_init_data(end,3:end)];
        cd ../        
        % Uses final values from DA sim as ICs
        
        opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');
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