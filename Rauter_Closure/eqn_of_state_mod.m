function eqn_of_state_mod
N = 200; % number of discretisation points in z for each of tau, u
h = 4e-2; % layer height (m)
d=1.43e-5; % grain diameter (m)

mu1_I=0.342; % \frac{u+u_d}{\rho_f \phi_c}
mu2_I=0.557; % 
I_0 = 0.069;

mu1_Iv = 0.32;
mu2_Iv = 0.7;
Iv_0 = 0.005;

reg_param = 10^8;

phi_c=0.585; % Volume fraction
eta_f = 0.0010016; % Pa s
g=9.81; % m/s^2
rho_f = 1000; % kg/m^3
rho_p = 2500; % kg/m^3
theta = 8.5; % deg

shear_lim = 5;
phi_rlp = 0.53;
phi_rcp = 0.63;
a = 20;

dz = h/(N-0.5); % z spacing
z_pe = linspace(dz/2,h,N);
z_u = linspace(0,h-dz/2,N);
p_b = (rho_p-rho_f)*g*phi_c*cosd(theta)*(h-z_pe);
buoyancy = -(rho_p-rho_f)*g*phi_c*cosd(theta);

v_scale = sqrt(g.*h);
p_scale = rho_f.*g.*h;
t_scale = sqrt(h./g);
z_scale = h;
density_ratio = rho_p./rho_f;
rho = phi_c*density_ratio+1-phi_c;

d_dl = d/z_scale;
dz_dl = dz/z_scale;
eta_f_dl = eta_f/(p_scale*t_scale);
a_dl = a/p_scale;
shear_lim_dl = shear_lim*t_scale;
z_pe_dl = z_pe/z_scale;
z_u_dl = z_u/z_scale;
p_b_dl = p_b/p_scale;
buoyancy_dl = buoyancy/p_scale*z_scale;

t_step = 0.1;
fname = "rau_8_5_deep_no_pe.txt";
run_no_pe()

    function dvecdt=no_excess_pressure(t,vec,phi_hat_const)
        len = length(vec);
        u_f = vec(1:len/2)';
        u_p = vec(len/2+1:end)';
        
        dufdz=diff(u_f)./dz_dl;
        dupdz=diff(u_p)./dz_dl;
        
        p_p = p_b_dl;

        Iv = eta_f_dl.*abs(dupdz)./(p_p(1:end-1));
%         I_u = interp1(z_pe,I,z_u,'linear','extrap');
        
        beta_pe = beta_fn(phi_hat_const);
        beta_u = interp1(z_pe,beta_pe,z_u,'linear','extrap');
       
        tau_f = [eta_f_dl.* dufdz 0];
        tau_p = [mu_Iv_fn(Iv).*sign(dupdz).*p_p(1:end-1) 0];
%         diff_tau_p = (mu_I_fn(I).*[0 diff(sign(dupdz).*p_p)]+(abs(I_u)>10/reg_param).*sign(dupdz).*p_p.*[0 diff(mu_I_fn(I))])./dz_dl;
        drag_force = (1-phi_c)^2.*beta_u.*(u_f-u_p);
        dufdt = [0 (1/(1-phi_c).*diff(tau_f)./dz_dl + sind(theta)-drag_force(2:end)./(1-phi_c))];
        dupdt = [0 (1/(density_ratio*phi_c).*diff(tau_p)./dz_dl + sind(theta)+drag_force(2:end)./(density_ratio*phi_c))];
        dvecdt = [dufdt dupdt]';
    end

    function dvecdt=const_up(t,vec,const_up)
        phi_hat = vec(1:end)';       
        
        const_dupdz=diff(const_up)./dz_dl;
        
        p_c = a_dl*(phi_c+phi_hat-phi_rlp)./(phi_rcp-phi_c-phi_hat);
        p_c(p_c<0) = 0;
        phi_m = phi_c+(phi_rcp-phi_c).*(abs(const_dupdz)<shear_lim_dl).*(shear_lim_dl-abs(const_dupdz)).^2/shear_lim_dl.^2;
        p_i = eta_f_dl.*const_dupdz./((phi_m./(phi_c+phi_hat(1:end-1))-1).^2);
        
        p_p = [p_c(1:end-1) + p_i 0];
        p_e = p_b_dl-p_p;
        dpdz = [0 diff(p_e)./dz_dl];

%         I_u = interp1(z_pe,I,z_u,'linear','extrap');
        
        beta_pe = beta_fn(phi_hat);
        beta_u = interp1(z_pe,beta_pe,z_u,'linear','extrap');
        
        dphidt = -phi_c.*diff((1./beta_u.*dpdz))./dz_dl;
        dphidt = interp1(z_pe(1:end-1),dphidt,z_pe,'linear','extrap');
       
        dvecdt = dphidt';
    end

    function success = run_no_pe()
        cd EqnOfState_Results
        init_data = load('rau_13_deep.txt');
%         steady_state = init_data(5001,:);
        cd ../
%         vec=zeros(N,1);
%         whole_data = ;
%         for j=1:N
%             vec(j,1) = (phi_rcp*p_b_dl(j)+a_dl*phi_rlp)/(a_dl+p_b_dl(j))-phi_c;
%             vec(N+j,1) = 0;
%             vec(2*N+j,1) = 0;
%         end
        vec = init_data(5000,N+1:end);
        final_phi = init_data(5000,1:N);
        % No initial pressure of phihat and initial values of u_p and u_f
        % defined above
        time_vals = (0:5000)*t_step;
        opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');

        [~,vec]=ode15s(@(t,vec) no_excess_pressure(t,vec,final_phi),time_vals,vec,opts);
        size(vec)
        save(fname, 'vec','-ascii')
        success=1;
    end
    
    function success = run_const_up()
        cd EqnOfState_Results
        init_data = load('rau_8_5_deep.txt');
%         steady_state = init_data(5001,:);
        cd ../
%         vec=zeros(2*N,1);
%         whole_data = ;
%         for j=1:N
%             vec(j,1) = (phi_rcp*p_b_dl(j)+a_dl*phi_rlp)/(a_dl+p_b_dl(j))-phi_c;
%             vec(N+j,1) = 0;
%             vec(2*N+j,1) = 0;
%         end
        vec = init_data(5,1:N);
        % No initial pressure of phihat and initial values of u_p and u_f
        % defined above
        up_max = init_data(5,3*N);
        const_u_p = (1-(1-z_u_dl).^(3/2))*up_max;
        time_vals = (0:5000)*t_step;
        opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');

        [~,vec]=ode15s(@(t,vec) const_up(t,vec,const_u_p),time_vals,vec,opts);
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