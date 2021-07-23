function simple_submerged_system
opengl software
N = 200; % number of discretisation points in z for each of tau, u
h = 4e-3; % layer height (m)
d=1.43e-4; % grain diameter (m)

mu1_I=0.342; % \frac{u+u_d}{\rho_f \phi_c}
mu2_I=0.557; % 
I_0 = 0.069;

phi_c=0.6; % Volume fraction
eta_f = 0.0010016; % Pa s
g=9.81; % m/s^2
rho_f = 1000; % kg/m^3
rho_p = 2500; % kg/m^3
theta = 9.5; % deg
alpha = 0.001; % 1/Pa

buoyancy = -(rho_p-rho_f)*g*phi_c*cosd(theta);

dz = h/(N-0.5); % z spacing
z_pe = linspace(dz/2,h,N);
z_u = linspace(0,h-dz/2,N);

p_b = (rho_p-rho_f)*g*phi_c*cosd(theta)*(h-z_pe);


fname = "pconst_5deg_params1.txt";
run_const_p_sim()

    % The function that defines the evolution of p_e, u_p and u_f for the case
    % where the evolution of p_e is purely diffusive
    function dvecdt=p_e_derivs(t,vec)
        len = length(vec);
        p_e = vec(1:len/3)';
        u_f = vec(len/3+1:2*len/3)';
        u_p = vec(2*len/3+1:end)';

        dpdz=[0 diff(p_e)./dz]; % transposes u and then finds differences to approx deriv
        dufdz=[diff(u_f)./dz 0];
        dupdz=[diff(u_p)./dz 0];
        
        p_p = p_b - p_e;
        phi_hat = p_p.*alpha;

         
        I = 2.*d.*dupdz.*sqrt(rho_p)./sqrt(p_p+0.000001);
        
        beta_pe = beta_fn(phi_hat);
        beta_u = interp1(z_pe,beta_pe,z_u,'linear','extrap');
        dpdt = [phi_c./alpha.*diff((1./beta_u.*dpdz))./dz 0];
       
        tau_f = eta_f.* dufdz;
        tau_p = mu_I_fn(I).*p_p;
        drag_force = (1-phi_c)^2.*beta_u.*(u_f-u_p);
        dufdt = [0 1/(rho_f*(1-phi_c)).*diff(tau_f)./dz+g*sind(theta)-drag_force(2:end)./(rho_f*(1-phi_c))];
        dupdt = [0 1/(rho_p*phi_c).*diff(tau_p)./dz + g*sind(theta)+drag_force(2:end)./(rho_p*phi_c)];
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

%         dphidz_u = diff(phi_hat)./dz; % N-1 elements as do not care about value at base
%         dphidz_h = interp1(z_u(2:end),dphidz_u,h,'linear','extrap');
        
        dpdz=[0 diff(p_e)./dz]; % transposes u and then finds differences to approx deriv
        dufdz=[diff(u_f)./dz 0];
        dupdz=[diff(u_p)./dz 0];
        d2updz2=[0 diff(dupdz)./dz];

        
        p_p = p_b - p_e;

        Iv = eta_f.*abs(dupdz)./(p_p+1e-8);
        Iv(end) = eta_f.*d2updz2(end)./(buoyancy-dpdz(end));
        I = 2.*d.*abs(dupdz).*sqrt(rho_p)./sqrt(p_p+1e-8);
        
        beta_pe = beta_fn(phi_hat);
        beta_u = interp1(z_pe,beta_pe,z_u,'linear','extrap');

        dilatancy = 1./alpha.*abs(dupdz).*(phi_hat+ (sqrt(abs(Iv)).*phi_c./(1+sqrt(abs(Iv)))));
        dpdt = [1./alpha.*diff((1./beta_u.*dpdz))./dz-dilatancy(1:end-1) 0];
        dphidt = -phi_c.*diff((1./beta_u.*dpdz))./dz;
        dphidt = interp1(z_pe(1:end-1),dphidt,z_pe,'linear','extrap');
       
        tau_f = eta_f.* dufdz;
        tau_p = sign(dupdz).*mu_I_fn(I).*p_p;
        drag_force = (1-phi_c)^2.*beta_u.*(u_f-u_p);
        dufdt = [0 1/(rho_f*(1-phi_c)).*diff(tau_f)./dz] + g*sind(theta)-drag_force./(rho_f*(1-phi_c));
        dupdt = [0 1/(rho_p*phi_c).*diff(tau_p)./dz + g*sind(theta)+drag_force(2:end)./(rho_p*phi_c)];
        dvecdt = [dpdt dphidt dufdt dupdt]';
    end

    % The function that defines evolution of p_e and phihat with constant
    % velocities in the dilatant case
    function dvecdt=p_phi_derivs(t,vec)
        len = length(vec);
        p_e = vec(1:len/2)';
        phi_hat = vec(len/2+1:end)';

        dpdz=[0 diff(p_e)./dz]; % transposes u and then finds differences to approx deriv
        
        p_p = p_b - p_e;

        Iv = eta_f.*abs(init_dupdz)./(p_p+1e-8);
        
        beta_pe = beta_fn(phi_hat);
        beta_u = interp1(z_pe,beta_pe,z_u,'linear','extrap');

        dilatancy = 1./alpha.*init_dupdz.*(phi_hat+sqrt(abs(Iv)).*phi_c./(1+sqrt(abs(Iv))));
        dpdt = [1./alpha.*diff((1./beta_u.*dpdz))./dz-dilatancy(1:end-1) 0];
        dphidt = [-phi_c.*diff((1./beta_u.*dpdz))./dz -phi_c*dilatancy(end)];
        dvecdt = [dpdt dphidt]';
    end
    
    % The function that defines evolution of u_p and u_f with a constant
    % pressure profile in the diffusion only case.
    function dvecdt=p_const_derivs(t,vec)
        len = length(vec);
        u_f = vec(1:len/2)';
        u_p = vec(len/2+1:end)';

        dufdz=[diff(u_f)./dz 0];
        dupdz=[diff(u_p)./dz 0];
        
        p_p = p_b;
        phi_hat = 0.4*p_p;
        I = 2.*d.*dupdz.*sqrt(rho_p)./sqrt(p_p+0.000001);
        beta_pe = beta_fn(phi_hat);
        beta_u = interp1(z_pe,beta_pe,z_u,'linear','extrap');
        
        tau_f = eta_f.* dufdz;
        tau_p = mu_I_fn(I).*p_p;
        drag_force = (1-phi_c)^2.*beta_u(2:end).*(u_f(2:end)-u_p(2:end));
        dufdt = [0 1/(rho_f*(1-phi_c)).*diff(tau_f)./dz+g*sind(theta)-drag_force./(rho_f*(1-phi_c))];
        dupdt = [0 1/(rho_p*phi_c).*diff(tau_p)./dz + g*sind(theta)+drag_force./(rho_p*phi_c)];
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
        
        time_vals = [0,0.01,0.05,0.1,0.5,1.0,5.0,100];
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
        time_vals = (0:400);

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
        init_vec=zeros(2*N,1); 
        % particles are initially stationary

        init_time_vals = [0,0.01,0.05,0.1,0.5,1.0,5.0,100.0];
        opts=odeset('AbsTol',1e-10,'RelTol',1e-10,'Stats','on');
        
        [~,init_vec]=ode15s(@p_const_derivs,init_time_vals,init_vec,opts);
        save('initial_velocity.txt', 'init_vec','-ascii')
        % Saves initial velocity in a file for reference
        init_uf=init_vec(7,1:N);
        init_up=init_vec(7,N+1:end);
        vec=zeros(4*N,1);
%         whole_data = load('dil9_5deg_small.txt');
        for j=1:N
            vec(j,1) = 0; %-mu1_I.*g*sind(theta)*(rho_p*phi_c+(1-phi_c)*rho_f)*(h-z_pe(j));
            vec(N+j,1) = 0;
            vec(2*N+j,1) = init_uf(j);
            vec(3*N+j,1) = init_up(j);
        end

        % No initial pressure of phihat and initial values of u_p and u_f
        % defined above
        time_vals = (0:4000);
        opts=odeset('AbsTol',1e-10,'RelTol',1e-10,'Stats','on');

        [~,vec]=ode15s(@dilatancy_derivs,time_vals,vec,opts);
        size(vec)
        save(fname, 'vec','-ascii')
        success=1;
    end

    % Runs the diffusion only system
    function success = run_p_e_sim()
        % Sets pressure and velocity initial conditions
        vec=zeros(3*N,1);
        for j=1:N
            vec(j,1) = 0.6*abs(p_b(j));
            vec(N+j,1) = 0;
            vec(2*N+j,1) = 0;
        end
        time_vals = (0:100);

        opts=odeset('AbsTol',1e-10,'RelTol',1e-10,'Stats','on');

        [~,vec]=ode15s(@p_e_derivs,time_vals,vec,opts);
        size(vec)
        save(fname, 'vec','-ascii')
        success=1;
    end

    function beta_val=beta_fn(phihat)
        beta_val = 150*(phi_c + phihat).^2.*eta_f./((1-phi_c-phihat).^3.*d^2);
    end

    function beta_val=dimless_beta_fn(phihat)
        beta_val = (phi_c + phihat).^2./((1-phi_c-phihat).^3);
    end
    
    function deriv_val=beta_deriv(phihat)
        deriv_val = 150*eta_f/d^2.*(phi_c + phihat).*(2+(phi_c + phihat))./(1-(phi_c + phihat)).^4;
    end

    function mu_val = mu_I_fn(I)
            mu_val = tanh(1e8*I).*(mu1_I+(mu2_I-mu1_I)./(1+I_0./abs(I)));
    end
end