function simple_submerged_system
    N = 200; % number of discretisation points in z for each of tau, u
    h = 4e-2; % layer height (m)
    d=1.43e-5; % grain diameter (m)

    mu1_I=0.342; % \frac{u+u_d}{\rho_f \phi_c}
    mu2_I=0.557; % 
    I_0 = 0.069;

    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    mu1_K = 0.2764;
    mu2_K = 0.8797;
    K_0 = 0.1931;
    a_K = 0.7071;
    reg_param = 10^8;

    phi_c=0.585; % Volume fraction
    eta_f = 0.0010016; % Pa s
    g=9.81; % m/s^2
    rho_f = 1000; % kg/m^3
    rho_p = 2500; % kg/m^3
    theta = 4; % deg
% for l=1:9
%     q= 9; %fix(l/10);
%     r=l; %mod(l,10);
%     fname = "Ive_comp_4";
%     if (r==0)
%         fname = append(fname,+"_deep_9"+"_start.txt");
%     else
%         fname = append(fname,"_deep_9_"+num2str(r)+"_start.txt");
%     end
%     fname = "Ive_comp_4_deep_v2.txt";
%     theta = 9+0.1*r;
    init_file_name = "Ive_comp_9_2_deep.txt";
    alpha = 0.0001; % 1/Pa
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
    alpha_dl = alpha*p_scale;
    z_pe_dl = z_pe/z_scale;
    z_u_dl = z_u/z_scale;
    p_b_dl = p_b/p_scale;
    buoyancy_dl = buoyancy/p_scale*z_scale;
    t_step = 1;
    fname = "Ive_comp_4_deep_9_2_start_early.txt";
    run_dil_sim()
    movefile(fname,'Results/');
    write_record(fname,'dil',N,h,d,reg_param,density_ratio,alpha_dl,phi_c,theta,eta_f_dl,t_step);
% end


% fname = "dil_8_5deg_deep_comp_new_phic.txt";
% run_dil_sim()
% movefile(fname,'Results/');


    % The function that defines the evolution of p_e, u_p and u_f for the case
    % where the evolution of p_e is purely diffusive
    function dvecdt=p_e_derivs(t,vec)
        len = length(vec);
        p_e = vec(1:len/3)';
        u_f = vec(len/3+1:2*len/3)';
        u_p = vec(2*len/3+1:end)';

        dpdz=[0 diff(p_e)./dz_dl]; % transposes u and then finds differences to approx deriv
        dufdz=diff(u_f)./dz_dl;
        dupdz=diff(u_p)./dz_dl;
        
        p_p = p_b_dl - p_e;
        phi_hat = p_p.*alpha_dl;

        I = 2.*d_dl.*dupdz.*sqrt(density_ratio)./sqrt(p_p(1:end-1));
        
        beta_pe = beta_fn(phi_hat);
        beta_u = interp1(z_pe_dl,beta_pe,z_u_dl,'linear','extrap');
        dpdt = [phi_c./alpha_dl.*diff((1./beta_u.*dpdz))./dz_dl 0];
       
        tau_f = [eta_f_dl.* dufdz 0];
        tau_p = [mu_I_fn(I).*p_p(1:end-1) 0];
        drag_force = (1-phi_c)^2.*beta_u.*(u_f-u_p);
        dufdt = [0 (1/(1-phi_c).*diff(tau_f)./dz_dl+sind(theta)-drag_force(2:end)./(1-phi_c))];
        dupdt = [0 (1/(density_ratio*phi_c).*diff(tau_p)./dz_dl + sind(theta)+drag_force(2:end)./(density_ratio*phi_c))];
        dvecdt = [dpdt dufdt dupdt]';
    end

    % The function that defines the evolution of p_e, phihat, u_p and u_f
    % including dilatancy.
    function dvecdt=dilatancy_derivs(t,vec)
        len = length(vec);
        p_e = vec(1:len/4)';
        phi_hat = vec(len/4+1:len/2)';
        u_f = vec(len/2+1:3*len/4)';
        u_p = vec(3*len/4+1:end)';
        
        dpdz=[0 diff(p_e)./dz_dl]; % transposes u and then finds differences to approx deriv
        dufdz=diff(u_f)./dz_dl;
        dupdz=diff(u_p)./dz_dl;
        d2updz2=[0 diff(dupdz)./dz_dl];
        
        p_p = p_b_dl - p_e;

        Iv = eta_f_dl.*abs(dupdz)./(p_p(1:end-1));
%         Iv(end) = eta_f_dl.*d2updz2(end)./(buoyancy_dl-dpdz(end));
        I = 2.*d_dl.*abs(dupdz).*sqrt(density_ratio)./sqrt(abs(p_p(1:end-1)));
%         I_u = interp1(z_pe,I,z_u,'linear','extrap');
        K = sqrt(I.^2+2*Iv);
        
        beta_pe = beta_fn(phi_hat);
        beta_u = interp1(z_pe,beta_pe,z_u,'linear','extrap');

        dilatancy = 1./alpha_dl.*(dupdz>5e-4).*abs(dupdz).*(phi_hat(1:end-1)+ (sqrt(abs(Iv)).*phi_c./(1+sqrt(abs(Iv)))));
        dpdt = [1./alpha_dl.*diff((1./beta_u.*dpdz))./dz_dl-dilatancy 0];
        dphidt = -phi_c.*diff((1./beta_u.*dpdz))./dz_dl;
        dphidt = interp1(z_pe(1:end-1),dphidt,z_pe,'linear','extrap');
       
        tau_f = [eta_f_dl.* dufdz 0];
        tau_p = [mu_Iv_fn(Iv).*sign(dupdz).*p_p(1:end-1) 0];
%         diff_tau_p = (mu_I_fn(I).*[0 diff(sign(dupdz).*p_p)]+(abs(I_u)>10/reg_param).*sign(dupdz).*p_p.*[0 diff(mu_I_fn(I))])./dz_dl;
        drag_force = (1-phi_c)^2.*beta_u.*(u_f-u_p);
        dufdt = [0 (1/(1-phi_c).*diff(tau_f)./dz_dl + sind(theta)-drag_force(2:end)./(1-phi_c))];
        dupdt = [0 (1/(density_ratio*phi_c).*diff(tau_p)./dz_dl + sind(theta)+drag_force(2:end)./(density_ratio*phi_c))];
        dvecdt = [dpdt dphidt dufdt dupdt]';
    end

    function dvecdt=single_phase_derivs(t,vec)
        len = length(vec);
        p_e = vec(1:len/3)';
        phi_hat = vec(len/3+1:2*len/3)';
        u = vec(2*len/3+1:end)';
        
        dpdz=[0 diff(p_e)./dz_dl]; % transposes u and then finds differences to approx deriv
        dudz=[diff(u)./dz_dl 0];
        d2udz2=[0 diff(dudz)./dz_dl];
        
        p_p = p_b_dl - p_e;

        Iv = eta_f_dl.*abs(dudz)./(p_p+1e-8);
        Iv(end) = eta_f_dl.*d2udz2(end)./(buoyancy_dl-dpdz(end));
        I = 2.*d_dl.*abs(dudz).*sqrt(density_ratio)./sqrt(p_p+1e-8);
%         I_u = interp1(z_pe,I,z_u,'linear','extrap');
        
        beta_pe = beta_fn(phi_hat);
        beta_u = interp1(z_pe,beta_pe,z_u,'linear','extrap');

        dilatancy = 1./alpha_dl.*abs(dudz).*(phi_hat+ (sqrt(abs(Iv)).*phi_c./(1+sqrt(abs(Iv)))));
        dpdt = [1./alpha_dl.*diff((1./beta_u.*dpdz))./dz_dl-dilatancy(1:end-1) 0];
        dphidt = -phi_c.*diff((1./beta_u.*dpdz))./dz_dl;
        dphidt = interp1(z_pe(1:end-1),dphidt,z_pe,'linear','extrap');
       
        tau_f = eta_f_dl.* dudz;
        tau_p = mu_I_fn(I).*sign(dudz).*p_p;
        dudt = [0 (1/(rho).*(diff(tau_f)./dz_dl + diff(tau_p)./dz_dl)+ sind(theta))];
        dvecdt = [dpdt dphidt dudt]';
    end

    % The function that defines evolution of p_e and phihat with constant
    % velocities in the dilatant case
    function dvecdt=p_phi_derivs(t,vec)
        len = length(vec);
        p_e = vec(1:len/2)';
        phi_hat = vec(len/2+1:end)';

        dpdz=[0 diff(p_e)./dz_dl]; % transposes u and then finds differences to approx deriv
        
        p_p = p_b_dl - p_e;

        Iv = eta_f_dl.*abs(init_dupdz)./(p_p+1e-8);
        
        beta_pe = beta_fn(phi_hat);
        beta_u = interp1(z_pe_dl,beta_pe,z_u_dl,'linear','extrap');

        dilatancy = 1./alpha_dl.*init_dupdz.*(phi_hat+sqrt(abs(Iv)).*phi_c./(1+sqrt(abs(Iv))));
        dpdt = [1./alpha_dl.*diff((1./beta_u.*dpdz))./dz_dl-dilatancy(1:end-1) 0];
        dphidt = [-phi_c.*diff((1./beta_u.*dpdz))./dz_dl -phi_c*dilatancy(end)];
        dvecdt = [dpdt dphidt]';
    end
    
    % The function that defines evolution of u_p and u_f with a constant
    % pressure profile in the diffusion only case.
    function dvecdt=p_const_derivs(t,vec)
        len = length(vec);
        u_f = vec(1:len/2)';
        u_p = vec(len/2+1:end)';

        dufdz=[diff(u_f)./dz_dl 0];
        dupdz=[diff(u_p)./dz_dl 0];
        
        p_p = 0.6*p_b_dl;
        phi_hat = alpha_dl*p_p;
        I = 2.*d_dl.*dupdz.*sqrt(density_ratio)./sqrt(p_p+0.000001);
        beta_pe = beta_fn(phi_hat);
        beta_u = interp1(z_pe_dl,beta_pe,z_u_dl,'linear','extrap');
        
        tau_f = eta_f.* dufdz;
        tau_p = mu_I_fn(I).*p_p;
        drag_force = (1-phi_c)^2.*beta_u(2:end).*(u_f(2:end)-u_p(2:end));
        dufdt = [0 1/(1-phi_c).*diff(tau_f)./dz_dl+sind(theta)-drag_force./(1-phi_c)];
        dupdt = [0 1/(density_ratio*phi_c).*diff(tau_p)./dz_dl + sind(theta)+drag_force./(density_ratio*phi_c)];
        dvecdt = [dufdt dupdt]';
    end

    % Runs the simulation for the constant pressure case
    function success = run_const_p_sim()
        vec=zeros(2*N,1); 
        for j=1:N
            vec(j,1) = 0;
            vec(N+j,1) = 0;
        end
        % particles are initially stationary
        
        time_vals = (0:100);
        % time values of interest
        opts=odeset('AbsTol',1e-10,'RelTol',1e-10,'Stats','on');
        
        [~,vec]=ode15s(@p_const_derivs,time_vals,vec,opts);
        save(fname, 'vec','-ascii')
        success=1;
    end

    % Runs the simulation where only p_e and phi evolve
    function success = run_p_phi_sim()
        vec=zeros(2*N,1); 
        for j=1:N
            vec(j,1) = 0;
            vec(N+j,1) = 0;
        end
        % Sets initial pressure and phihat values
        time_vals = (0:10)*100;

        opts=odeset('AbsTol',1e-10,'RelTol',1e-10,'Stats','on');

        [~,vec]=ode15s(@p_phi_derivs,time_vals,vec,opts);
        size(vec)
        save(fname, 'vec','-ascii')
        success=1;
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
        cd Results
%         init_file_name = 'Ive_comp_13_deep_custom_time.txt';
        init_data = load(init_file_name);
%         steady_state = init_data(5001,:);
        cd ../
        vec=init_data(5001,:);
%         vec = zeros(4*N,1);
%         whole_data = ;
%         for j=1:N
% %             vec(1,N+j)=vec(1,N+j)+0.015;
%             vec(j,1) = p_b_dl(j)*0.4;
%             vec(N+j,1) = p_b_dl(j)*0.6*alpha_dl;
%             vec(2*N+j,1) = 0; %(1-(1-z_u_dl(j))^(3/2)); %init_uf(j);
%             vec(3*N+j,1) = 0; %(1-(1-z_u_dl(j))^(3/2)); %init_up(j);
%         end

        % No initial pressure of phihat and initial values of u_p and u_f
        % defined above
        time_vals = [0 250];
        opts=odeset('AbsTol',1e-8,'RelTol',1e-8,'Stats','on');

        [tvals,vec]=ode15s(@dilatancy_derivs,time_vals,vec,opts);
        vec = horzcat(tvals,vec);
        size(vec)
        save(fname, 'vec','-ascii')
        success=1;
    end

    % Runs the single phase simulation
    function success = run_single_phase_sim()
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
        cd Results
        init_data = load("dil15deg_acc.txt");
        steady_state = init_data(6000,:);
        cd ../
        vec=zeros(3*N,1);
%         whole_data = load('dil9_5deg_small.txt');
        for j=1:N
            vec(j,1) = 0;%0.6.*sind(theta)*phi_c*(1-z_pe_dl(j));
            vec(N+j,1) = steady_state(200+j);
            vec(2*N+j,1) = steady_state(400+j); %init_uf(j);
        end

        % No initial pressure of phihat and initial values of u_p and u_f
        % defined above
        time_vals = (0:3000)*t_step;
        opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');

        [~,vec]=ode15s(@single_phase_derivs,time_vals,vec,opts);
        size(vec)
        save(fname, 'vec','-ascii')
        success=1;
    end

    % Runs the diffusion only system
    function success = run_p_e_sim()
%         init_vec=zeros(2*N,1); 
        % particles are initially stationary

%         init_time_vals = [0,5.0,5000.0];
%         opts=odeset('AbsTol',1e-10,'RelTol',1e-10,'Stats','on');
%         
%         [~,init_vec]=ode15s(@p_const_derivs,init_time_vals,init_vec,opts);
%         init_uf=init_vec(3,1:N);
%         init_up=init_vec(3,N+1:end);
        % Sets pressure and velocity initial conditions
        vec=zeros(3*N,1);
        for j=1:N
            vec(j,1) = 0.6*p_b_dl(j);
            vec(N+j,1) = 0;%init_uf(j);
            vec(2*N+j,1) = 0;%init_up(j);
        end
        time_vals = (0:400)*t_step;

        opts=odeset('AbsTol',1e-10,'RelTol',1e-10,'Stats','on');

        [~,vec]=ode15s(@p_e_derivs,time_vals,vec,opts);
        size(vec)
        save(fname, 'vec','-ascii')
        success=1;
    end

    function beta_val=beta_fn(phihat)
        beta_val = 150*(phi_c + phihat).^2.*eta_f_dl./((1-phi_c-phihat).^3.*d_dl^2);
    end

    function beta_val=dimless_beta_fn(phihat)
        beta_val = (phi_c + phihat).^2./((1-phi_c-phihat).^3);
    end
    
    function deriv_val=beta_deriv(phihat)
        deriv_val = 150*eta_f/d^2.*(phi_c + phihat).*(2+(phi_c + phihat))./(1-(phi_c + phihat)).^4;
    end

    function mu_val = mu_I_fn(I)
        mu_val = tanh(reg_param*I).*(mu1_I+(mu2_I-mu1_I)./(1+I_0./abs(I)));
    end

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end

    function mu_val = mu_K_fn(K,Iv)
        mu_val = tanh(reg_param*K).*(mu1_K+(mu2_K-mu1_K)./(1+K_0./abs(K))+5/2*phi_c/(a_K)*(Iv./(K+1e-8)));
    end
    
end