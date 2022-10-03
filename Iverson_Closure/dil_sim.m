function dil_sim
    N = 200; % number of discretisation points in z for each of tau, u
    h = 4e-2; % layer height (m)
    d= 1e-5; % grain diameter (m)

    mu1_I=0.342; % \frac{u+u_d}{\rho_f \phi_c}
    mu2_I=0.557; % 
    I_0 = 0.069;
    
%     mu1_I=tand(20); 
%     mu2_I=tand(33);  
%     I_0 = 0.3;
%     del_phi = 0.2;
    
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    
    eta_f = 0.0010016; % Pa s
    rho_f = 1000; % kg/m^3
    
%     eta_f = 1.18e-5; % Pa s
%     rho_f = 1; % kg/m^3
    
    g=9.81; % m/s^2
    rho_p = 2500; % kg/m^3
    theta = 5; % deg
    theta0 = 13;

%     init_file_name = "dil_13deg_deep_comp_new_phic.txt";
    alpha = 1e-4; % 1/Pa
    dz = h/(N-0.5); % z spacing
    z_pe = linspace(dz/2,h,N);
    z_u = linspace(0,h-dz/2,N);

    v_scale = sqrt(g.*h);
    p_scale = rho_f.*g.*h;
    t_scale = sqrt(h./g);
    z_scale = h;
    density_ratio = rho_p./rho_f;
    rho = phi_c*density_ratio+1-phi_c;
%     crit_grad = -(density_ratio*phi_c+(1-phi_c))*sind(theta)/mu1_I;
    
    p_b = (rho_p-rho_f)*g*phi_c*cosd(theta)*(h-z_pe);
    buoyancy = -(rho_p-rho_f)*g*phi_c*cosd(theta);
        
    p_b_dl = p_b/p_scale;
    buoyancy_dl = buoyancy/p_scale*z_scale;
    d_dl = d/z_scale;
    dz_dl = dz/z_scale;
    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    z_pe_dl = z_pe/z_scale;
    z_u_dl = z_u/z_scale;
    
    init_Iv = newt_solve_crit_Iv(theta0,rho_p,rho_f);
    init_umax = init_Iv/eta_f_dl*(density_ratio-1)*phi_c*cosd(theta0)/2;
%     pe_crit = p_b_dl+crit_grad*(1-z_pe);
    t_step = 1;
    % for l=1:90
%     q= fix(l/10);
%     r=mod(l,10);
%     fname = "Ive_comp_"+num2str(q);
%     if (r==0)
%         fname = append(fname,"_deep.txt");
%     else
%         fname = append(fname,"_"+num2str(r)+"_deep.txt");
%     end
%     fname = "Ive_comp_4_deep_v2.txt";
    flux_con = load('Results/flux_conditions.txt');
%     for l=10:10
%     l = 20;
%         flux=flux_con(l,1);
%         theta_init=flux_con(l,2);
%         init_phi=flux_con(l,3);
%         unit = floor(flux);
%         tenth = mod(l,10);
%         if (tenth == 0)
% %             init_file_name = "Ive_comp_"+num2str(unit)+"_deep.txt";
%             fname = "Ive_comp_4_deep_"+num2str(unit)+"_flux_big_part.txt";
% %             fname = "Ive_comp_"+num2str(unit)+"_deep.txt";
%         else
% %             init_file_name = "Ive_comp_"+num2str(unit)+"_"+num2str(tenth)+"_deep.txt";
%             fname = "Ive_comp_4_deep_"+num2str(unit)+"_"+num2str(tenth)+"_flux_big_part.txt";
% %             fname = "Ive_comp_"+num2str(unit)+"_"+num2str(tenth)+"_deep.txt";
%         end
%         init_file_name = fname;
%     init_file_name = "Ive_comp_9_2_deep_long_run.txt";
        fname = "Ive_"+num2str(theta)+"deg_"+num2str(theta0)+"init_short_ts.txt";
        run_dil_sim()
        movefile(fname,'Results/');
        write_record(fname,'dil',N,h,d,reg_param,density_ratio,alpha_dl,phi_c,theta,eta_f_dl,t_step);
%     end

    % The function that defines the evolution of p_e, phihat, u_p and u_f
    % including dilatancy.
    function dvecdt=dilatancy_derivs(t,vec)
        % Code for having a value of theta that changes linearly with time
        % at the start of the simulation before becoming constant
%         theta_start = 13;
%         theta_stop = 5;
%         theta_delta_t = 10;
%         
%         theta = max(theta_stop + (theta_start-theta_stop)*(theta_delta_t-t)/theta_delta_t,theta_stop);
        
        len = length(vec);
        p_e = vec(1:len/4)';
        phi = vec(len/4+1:len/2)';
        u_f = vec(len/2+1:3*len/4)';
        u_p = vec(3*len/4+1:end)';
        
        dpdz=[0 diff(p_e)./dz_dl]; % transposes u and then finds differences to approx deriv
        dufdz=diff(u_f)./dz_dl;
        dupdz=diff(u_p)./dz_dl;
        d2updz2=[0 diff(dupdz)./dz_dl];
        
        p_p = p_b_dl - p_e;

        Iv = eta_f_dl.*abs(dupdz)./(p_p(1:end-1));
%         Iv(end) = eta_f_dl.*d2updz2(end)./(buoyancy_dl-dpdz(end));
%         I = 2.*d_dl.*abs(dupdz).*sqrt(density_ratio)./sqrt(abs(p_p(1:end-1)));
%         I_u = interp1(z_pe,I,z_u,'linear','extrap');
%         K = sqrt(I.^2+2*Iv);
        
        beta_pe = beta_fn(phi);
        beta_u = interp1(z_pe,beta_pe,z_u,'linear','extrap');

        dilatancy = 1./alpha_dl.*(dupdz>1e-3).*abs(dupdz).*(phi(1:end-1)- phi_c./(1+sqrt(abs(Iv))));
%         dilatancy = 1./alpha_dl.*(dupdz>5e-4).*abs(dupdz).*(phi_hat(1:end-1) + (del_phi*I));
        
        dpdt = [1./alpha_dl.*diff((1./beta_u.*dpdz))./dz_dl-dilatancy 0];
        dphidt = -phi_c.*diff((1./beta_u.*dpdz))./dz_dl;
        dphidt = interp1(z_pe(1:end-1),dphidt,z_pe,'linear','extrap');
       
        tau_f = [eta_f_dl.* dufdz 0];
        tau_p = [mu_Iv_fn(abs(Iv)).*sign(dupdz).*p_p(1:end-1) 0];
%         diff_tau_p = (mu_I_fn(I).*[0 diff(sign(dupdz).*p_p)]+(abs(I_u)>10/reg_param).*sign(dupdz).*p_p.*[0 diff(mu_I_fn(I))])./dz_dl;
        drag_force = (1-phi_c)^2.*beta_u.*(u_f-u_p);
        dufdt = [0 (1/(1-phi_c).*diff(tau_f)./dz_dl + sind(theta)-drag_force(2:end)./(1-phi_c))];
        dupdt = [0 (1/(density_ratio*phi_c).*diff(tau_p)./dz_dl + sind(theta)+drag_force(2:end)./(density_ratio*phi_c))];
        dvecdt = [dpdt dphidt dufdt dupdt]';
    end

    % Runs the full dilatant system
    function success = run_dil_sim()
        % Start by running simulation for initial velocity profile. Want
%         % the steady state velocity at zero pressure and zero phihat
%         init_vec=zeros(2*N,1); 
%         % particles are initially stationary
% 
%         init_time_vals = [0,0.01,0.05,0.1,0.5,1.0,5.0,100.0];
%         opts=odeset('AbsTol',1e-10,'RelTol',1e-10,'Stats','on');
%         
%         [~,init_vec]=ode15s(@p_const_derivs,init_time_vals,init_vec,opts);
%         save('initial_velocity.txt', 'init_vec','-ascii')
%         % Saves initial velocity in a file for reference
%         init_uf=init_vec(7,1:N);
%         init_up=init_vec(7,N+1:end);
%         cd Results
%         init_file_name = 'Ive_comp_10_deep.txt';
%         init_data = load(init_file_name);
%         steady_state = init_data(5001,:);
%         cd ../
%         vec=init_data(end,end-799:end);
        vec = zeros(4*N,1);
%         whole_data = ;
        for j=1:N
%             vec(j,1) = p_b_dl(j)*-0.4;
            vec(j,1) = (cosd(theta0)-cosd(theta))*(1-z_u_dl(j));
            vec(N+j,1) = phi_c/(1+sqrt(init_Iv));
            vec(2*N+j,1) = (1-(1-z_u_dl(j))^2)*init_umax; %(1-(1-z_u_dl(j))^2)*3/2*flux+(1-phi_c).*sind(theta_init)./beta_fn(init_phi)./(1-phi_c)^2; %init_uf(j);
            vec(3*N+j,1) = (1-(1-z_u_dl(j))^2)*init_umax; %(1-(1-z_u_dl(j))^2).*3/2*flux; %init_up(j);
        end

        % No initial pressure of phihat and initial values of u_p and u_f
        % defined above
        time_vals = linspace(0,0.5);
        opts=odeset('AbsTol',1e-8,'RelTol',1e-8,'Stats','on');

        [tvals,vec]=ode15s(@dilatancy_derivs,time_vals,vec,opts);
        vec = horzcat(tvals,vec);
        size(vec)
        save(fname, 'vec','-ascii')
        success=1;
    end

    function beta_val=beta_fn(phi)
        beta_val = 150*phi.^2.*eta_f_dl./((1-phi).^3.*d_dl^2);
    end

    function mu_val = mu_I_fn(I)
        mu_val = tanh(reg_param*I).*(mu1_I+(mu2_I-mu1_I)./(1+I_0./abs(I)));
    end

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
end