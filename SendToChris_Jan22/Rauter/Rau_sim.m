function Rau_sim
N = 200; % number of discretisation points in z for each of tau, u
h = 4e-2; % layer height (m)
d=1.43e-5; % grain diameter (m)

mu1_I=0.342; 
mu2_I=0.557;  
I_0 = 0.069;

mu1_Iv = 0.32;
mu2_Iv = 0.7;
Iv_0 = 0.005;

reg_param = 10^8;

phi_c=0.585; % Volume fraction
eta_f = 0.0010016; % Fluid Viscosity (Pa s)
g=9.81; % Gravity (m/s^2)
rho_f = 1000; % Fluid Density (kg/m^3)
rho_p = 2500; % Particle Density (kg/m^3)
theta = 13; % Slope Angle (deg)

shear_lim = 5; % Shear Limit (1/s) - Max shear rate for material denser than phi_c
phi_rlp = 0.53; % Random loose packing
phi_rcp = 0.63; % Random close packing
a = 20; % Inverse of Particle Compressibility (Pa)

dz = h/(N-0.5); % z point spacing (m)
z_pe = linspace(dz/2,h,N);
z_u = linspace(0,h-dz/2,N);

p_b = (rho_p-rho_f)*g*phi_c*cosd(theta)*(h-z_pe); % Total Pressure (Pa)
buoyancy = -(rho_p-rho_f)*g*phi_c*cosd(theta);

v_scale = sqrt(g.*h); % Typical Scales
p_scale = rho_f.*g.*h;
t_scale = sqrt(h./g);
z_scale = h;

density_ratio = rho_p./rho_f;
rho = phi_c*density_ratio+1-phi_c;

d_dl = d/z_scale; % Dimensionless Quantities
dz_dl = dz/z_scale;
a_dl = a./p_scale;
shear_lim_dl = shear_lim*t_scale;
eta_f_dl = eta_f/(p_scale*t_scale);
z_pe_dl = z_pe/z_scale;
z_u_dl = z_u/z_scale;
p_b_dl = p_b/p_scale;
buoyancy_dl = buoyancy/p_scale*z_scale;

% t_step = 1; % Time between points when specifying time values
% num_points = 5000; % Num time values
% time_vals = (0:num_points)*t_step; % time values of interest
    
time_vals = [0,5000]; % Start and end time

fname = "Rauter_13_deep.txt";
run_rauter_sim()
movefile(fname,'Results/');
Rau_write_record(fname,N,h,d,reg_param,density_ratio,phi_c,theta,eta_f_dl,a_dl,phi_rcp,phi_rlp,shear_lim_dl);

    % The function that defines the evolution of p_e, phihat, u_p and u_f
    % with the Rauter closure eqution. There is no evolution eq
    function dvecdt=rauter_closure(t,vec)
        len = length(vec);
        phi_hat = vec(1:len/3)';
        u_f = vec(len/3+1:2*len/3)';
        u_p = vec(2*len/3+1:end)';
        
        dufdz=diff(u_f)./dz_dl;
        dupdz=diff(u_p)./dz_dl;
        
        p_c = a_dl*(phi_c+phi_hat-phi_rlp)./(phi_rcp-phi_c-phi_hat);
        p_c(p_c<0) = 0; % Ensures the contact pressure is positive
        phi_m = phi_c+(phi_rcp-phi_c).*(abs(dupdz)<shear_lim_dl).*(shear_lim_dl-abs(dupdz)).^2/shear_lim_dl.^2;
        p_i = eta_f_dl.*dupdz./((phi_m./(phi_c+phi_hat(1:end-1))-1).^2);
        
        p_p = [p_c(1:end-1) + p_i 0];
        p_e = p_b_dl-p_p;
        dpdz = [0 diff(p_e)./dz_dl];

        Iv = eta_f_dl.*abs(dupdz)./(p_p(1:end-1));
        
        beta_pe = beta_fn(phi_hat);
        beta_u = interp1(z_pe,beta_pe,z_u,'linear','extrap');
        
        dphidt = -phi_c.*diff((1./beta_u.*dpdz))./dz_dl;
        dphidt = interp1(z_pe(1:end-1),dphidt,z_pe,'linear','extrap');
       
        tau_f = [eta_f_dl.* dufdz 0];
        tau_p = [mu_Iv_fn(Iv).*sign(dupdz).*p_p(1:end-1) 0];

        drag_force = (1-phi_c)^2.*beta_u.*(u_f-u_p);
        dufdt = [0 (1/(1-phi_c).*diff(tau_f)./dz_dl + sind(theta)-drag_force(2:end)./(1-phi_c))];
        dupdt = [0 (1/(density_ratio*phi_c).*diff(tau_p)./dz_dl + sind(theta)+drag_force(2:end)./(density_ratio*phi_c))];
        dvecdt = [dphidt dufdt dupdt]';
    end

    function success = run_rauter_sim()
        cd Results
        init_data = load('Rauter_5_deep.txt');
        cd ../
        vec = init_data(end,:);
        % Takes the final results of a previous simulation as the ICs
        
        opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');
        [~,vec]=ode15s(@rauter_closure,time_vals,vec,opts);
        size(vec)
        save(fname, 'vec','-ascii')
        success=1;
    end

    function beta_val=beta_fn(phihat)
        beta_val = 150*(phi_c + phihat).^2.*eta_f_dl./((1-phi_c-phihat).^3.*d_dl^2);
    end
   
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
end