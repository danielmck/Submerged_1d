function eqn_of_state_system
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
theta = 5; % deg

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

t_step = 10;

% for l=1:4
%     q=fix(l/2);
%     r=mod(l,2);
    fname = "Rauter_6_1_deep.txt";
%     if (r==0)
%         fname = append(fname,"_deep.txt");
%     else
%         fname = append(fname,"_"+num2str(r)+"_deep.txt");
%     end
    theta = 6.1;
    run_rauter_sim()
    movefile(fname,'EqnOfState_Results/');
    EOS_write_record(fname,N,h,d,reg_param,density_ratio,phi_c,theta,eta_f_dl,a_dl,phi_rcp,phi_rlp,t_step,2,shear_lim_dl);
    % The function that defines the evolution of p_e, phihat, u_p and u_f
    % including dilatancy.
    function dvecdt=rauter_closure(t,vec)
        len = length(vec);
        phi_hat = vec(1:len/3)';
        u_f = vec(len/3+1:2*len/3)';
        u_p = vec(2*len/3+1:end)';
        
        dufdz=diff(u_f)./dz_dl;
        dupdz=diff(u_p)./dz_dl;
        
        p_c = a_dl*(phi_c+phi_hat-phi_rlp)./(phi_rcp-phi_c-phi_hat);
        p_c(p_c<0) = 0;
        phi_m = phi_c+(phi_rcp-phi_c).*(abs(dupdz)<shear_lim_dl).*(shear_lim_dl-abs(dupdz)).^2/shear_lim_dl.^2;
        p_i = eta_f_dl.*dupdz./((phi_m./(phi_c+phi_hat(1:end-1))-1).^2);
        
        p_p = [p_c(1:end-1) + p_i 0];
        p_e = p_b_dl-p_p;
        dpdz = [0 diff(p_e)./dz_dl];

        Iv = eta_f_dl.*abs(dupdz)./(p_p(1:end-1));
%         I_u = interp1(z_pe,I,z_u,'linear','extrap');
        
        beta_pe = beta_fn(phi_hat);
        beta_u = interp1(z_pe,beta_pe,z_u,'linear','extrap');
        
        dphidt = -phi_c.*diff((1./beta_u.*dpdz))./dz_dl;
        dphidt = interp1(z_pe(1:end-1),dphidt,z_pe,'linear','extrap');
       
        tau_f = [eta_f_dl.* dufdz 0];
        tau_p = [mu_Iv_fn(Iv).*sign(dupdz).*p_p(1:end-1) 0];
%         diff_tau_p = (mu_I_fn(I).*[0 diff(sign(dupdz).*p_p)]+(abs(I_u)>10/reg_param).*sign(dupdz).*p_p.*[0 diff(mu_I_fn(I))])./dz_dl;
        drag_force = (1-phi_c)^2.*beta_u.*(u_f-u_p);
        dufdt = [0 (1/(1-phi_c).*diff(tau_f)./dz_dl + sind(theta)-drag_force(2:end)./(1-phi_c))];
        dupdt = [0 (1/(density_ratio*phi_c).*diff(tau_p)./dz_dl + sind(theta)+drag_force(2:end)./(density_ratio*phi_c))];
        dvecdt = [dphidt dufdt dupdt]';
    end

    function success = run_rauter_sim()
        cd EqnOfState_Results
        init_data = load('rau_13_deep.txt');
%         steady_state = init_data(5001,:);
        cd ../
        vec=zeros(3*N,1);
%         whole_data = ;
%         for j=1:N
%             vec(j,1) = (phi_rcp*p_b_dl(j)+a_dl*phi_rlp)/(a_dl+p_b_dl(j))-phi_c;
%             vec(N+j,1) = 0;
%             vec(2*N+j,1) = 0;
%         end
        vec = init_data(5000,:);
        % No initial pressure of phihat and initial values of u_p and u_f
        % defined above
        time_vals = (0:5000)*t_step;
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