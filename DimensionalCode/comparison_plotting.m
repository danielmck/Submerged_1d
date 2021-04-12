opengl software
sim1 = load("../Results/slower_p_evo_data.txt");
% sim2 = load("dil9_5_Iv_large.txt");
% sim3 = load("sim3_out.txt");
sim_list = {sim1};
% Loads the data to compare and puts them in a list of simulations

sim_num = size(sim_list,2);
n_times = zeros(1,sim_num);
for j=1:sim_num
    n_times(1,j) = size(sim_list{1,j},1);
end

% Need to specify the type of simulation:
% dilatancy - "dil"
% diffusion only - "pdriv"
% constant pressure profile - "pcon"
% constant velocity profile - "ucon"
% Need to create lists storing the parameters in positions that corrrelate to
% the data in sim_list
sim_type = ["pdriv","dil","pdriv"];

N=[200,200,200];
h = [4e-3,4e-3,4e-3];
d=[1.43e-5, 1.43e-4,1.43e-4];
dz = 1./(N-0.5); % z spacing
phi_c=[0.6,0.6,0.6]; % Volume fraction
eta_f = [0.0010016,0.0010016,0.0010016]; % Pa s
g=[9.81,9.81,9.81]; % m/s^2
rho_f = [1000,1000,1000]; % kg/m^3
rho_p = [2500,2500,2500]; % kg/m^3
theta = [9.5,9.5,5.0]; % deg
alpha = [0.001,0.001,0.001]; % 1/Pa
s_frac = [0.0,0,0.7];

v_scale = sqrt(g.*h);
p_scale = rho_f.*g.*h;
t_scale = sqrt(h./g);
z_scale = h;
density_ratio = rho_p./rho_f;

% Creates cells to store the matrices of values

buoyancy = zeros(sim_num);
z_pe = cell(sim_num,1);
z_u = cell(sim_num,1);
p_b = cell(sim_num,1);
p_e = cell(sim_num,1);
p_p = cell(sim_num,1);
phi = cell(sim_num,1);
u_f = cell(sim_num,1);
u_p = cell(sim_num,1);
I = cell(sim_num,1);
Iv = cell(sim_num,1);
mu = cell(sim_num,1);
beta_pe = cell(sim_num,1);
beta_u = cell(sim_num,1);
tau_p = cell(sim_num,1);
tau_f = cell(sim_num,1);
drag_mult_p = cell(sim_num,1);
drag_term_p = cell(sim_num,1);
drag_mult_f = cell(sim_num,1);
drag_term_f = cell(sim_num,1);
dpdz = cell(sim_num,1);
d2pdz2 = cell(sim_num,1);
d3pdz3 = cell(sim_num,1);
dufdz = cell(sim_num,1);
d2ufdz2 = cell(sim_num,1);
d2updz2 = cell(sim_num,1);
dupdz = cell(sim_num,1);
dilatancy = cell(sim_num,1);
diffusion_term = cell(sim_num,1);
dpdt = cell(sim_num,1);
dphidt = cell(sim_num,1);
dphidz = cell(sim_num,1);
dupdt = cell(sim_num,1);
dufdt = cell(sim_num,1);

% loops over the simulations
for i = 1:sim_num
    buoyancy(i) = -(rho_p(i)-rho_f(i))*g(i)*phi_c(i)*cosd(theta(i))*h(i);

    z_pe{i,1} = linspace(1/(2*N(i)),1,N(i))';
    z_u{i,1} = linspace(0,1-1/(2*N(i)),N(i))';
    p_b{i,1} = (rho_p(i)-rho_f(i))*g(i)*phi_c(i)*cosd(theta(i))*(1-z_pe{i,1}).*h(i)./p_scale(i);
    
    vec = sim_list{1,i};
    % Processes the data differently depending on the simulation type
    if (sim_type(i) == "dil")
        p_e{i,1} = vec(:,1:N(i))'./p_scale(i);
        p_p{i,1} = p_b{i,1}-p_e{i,1};
        phi{i,1} = vec(:,N(i)+1:2*N(i))';
        u_f{i,1} = vec(:,2*N(i)+1:3*N(i))'./v_scale(i);
        u_p{i,1} = vec(:,3*N(i)+1:end)'./v_scale(i);
    elseif (sim_type(i) == "pcon")
        p_e{i,1} = s_frac(i)*p_b{i,1};
        p_p{i,1} = (1-s_frac(i))*p_b{i,1};
        % If looking at initial values for dilatancy sim, need to set phi
        % to zero here
        phi{i,1} = alpha(i)*p_p{i,1};
        u_f{i,1} = vec(:,1:N(i))'./v_scale(i);
        u_p{i,1} = vec(:,N(i)+1:end)'./v_scale(i);
    elseif (sim_type(i) == "pdriv")
        p_e{i,1} = vec(:,1:N(i))'./p_scale(i);
        p_p{i,1} = (p_b{i,1}-p_e{i,1});
        phi{i,1} = alpha(i)*p_p{i,1};
        u_f{i,1} = vec(:,N(i)+1:2*N(i))'./v_scale(i);
        u_p{i,1} = vec(:,2*N(i)+1:end)'./v_scale(i);
    elseif (sim_type(i) == "ucon")
        p_e{i,1} = vec(:,1:N(i))'./p_scale(i);
        p_p{i,1} = (p_b{i,1}-p_e{i,1})./p_scale(i);
        phi{i,1} = vec(:,N(i)+1:end)';
        initial_u = load("no_phi_no_p.txt");
        u_f{i,1} = initial_u(7,1:N(i))'./v_scale(i).*(ones(1,11));
        u_p{i,1} = initial_u(7,N(i)+1:end)'./v_scale(i).*(ones(1,11));
    end
    beta_pe{i,1} = 150*(phi_c(i) + phi{i,1}).^2.*eta_f(i)./((1-phi_c(i)-phi{i,1}).^3.*d(i)^2)*z_scale(i)^2./(p_scale(i)*t_scale(i));
    beta_u{i,1} = interp1(z_pe{i,1},beta_pe{i,1},z_u{i,1},'linear','extrap');
    % If the velocity evolves we need to define these quantities
    dpdz{i,1} = vertcat(zeros(1,n_times(i)),diff(p_e{i,1},1,1))./dz(i);
    d2pdz2{i,1} = vertcat(diff(dpdz{i,1},1,1),zeros(1,n_times(i)))./dz(i);
    d3pdz3{i,1} = vertcat(zeros(1,n_times(i)),diff(d2pdz2{i,1},1,1))./dz(i);
    dufdz{i,1} = vertcat(diff(u_f{i,1},1,1),zeros(1,n_times(i)))./dz(i);
    d2ufdz2{i,1} = vertcat(zeros(1,n_times(i)),diff(dufdz{i,1},1,1))./dz(i);
    dupdz{i,1} = vertcat(diff(u_p{i,1},1,1),zeros(1,n_times(i)))./dz(i);
    d2updz2{i,1} = vertcat(zeros(1,n_times(i)),diff(dupdz{i,1},1,1))./dz(i);
    Iv_temp = eta_f(i).*abs(dupdz{i,1})./(p_p{i,1}+1e-8)./z_scale(i).*v_scale(i)./p_scale(i);
    Iv_temp(N(i),:) = eta_f(i).*d2updz2{i,1}(N(i),:)./(-buoyancy(i)-dpdz{i,1}(N(i),:)).*v_scale(i)./p_scale(i);
    Iv{i,1} = Iv_temp;
    if (sim_type(i)~="ucon")
        tau_f{i,1} = eta_f(i)./(p_scale(i).*t_scale(i)).*d2ufdz2{i,1}./(1-phi_c(i));
    
        I{i,1} = 2.*d(i).*abs(dupdz{i,1}).*v_scale(i)./z_scale(i).*sqrt(rho_p(i))./sqrt(p_p{i,1}.*ones(1,n_times(i)))./sqrt(p_scale(i));
        mu{i,1} = mu_I_fn(I{i,1});
        tau_p{i,1} = 1/(phi_c(i).*density_ratio(i)).*vertcat(diff(mu{i,1}.*(p_p{i,1}.*ones(1,n_times(i))),1,1)./dz(i),zeros(1,n_times(i)));
    
        drag_mult_p{i,1} = (1-phi_c(i))^2.*beta_u{i,1}./(density_ratio(i)*phi_c(i)).*ones(1,n_times(i));
    
        drag_term_p{i,1} = drag_mult_p{i,1}.*(u_f{i,1}-u_p{i,1});
        
        drag_mult_f{i,1} = (1-phi_c(i)).*beta_u{i,1}.*ones(1,n_times(i));
    
        drag_term_f{i,1} = drag_mult_f{i,1}.*(u_f{i,1}-u_p{i,1});
        
        dupdt{i,1} = tau_p{i,1}+sind(theta(i))+drag_term_p{i,1};
        dufdt{i,1} = tau_f{i,1}+sind(theta(i))-drag_term_f{i,1};
    end
    
    
    % If the pressure evolves, have to define these quantities
    if ((sim_type(i)=="pdriv") || (sim_type(i)=="dil") || (sim_type(i)=="ucon"))
        dphidz{i,1} = vertcat(zeros(1,n_times(i)),diff(phi{i,1},1,1))./dz(i);
        dpdt{i,1} = phi_c(i)./(alpha(i).*p_scale(i)).*vertcat(diff(1./beta_pe{i,1}.*dpdz{i,1}),zeros(1,n_times(i)))./dz(i);
        if ((sim_type(i)=="dil") || (sim_type(i)=="ucon"))
            dilatancy{i,1} = 1./(alpha(i).*p_scale(i)).*dupdz{i,1}.*(phi{i,1}+sqrt(abs(Iv{i,1}))*phi_c(i)./(1+sqrt(abs(Iv{i,1}))));
            diffusion_term{i,1} = dpdt{i,1};
            dphidt{i,1} = -dpdt{i,1}.*p_scale(i).*alpha(i).*phi_c(i);
            dpdt{i,1} = dpdt{i,1} - dilatancy{i,1};
            d2pdzdt = vertcat(zeros(1,n_times(i)),diff(dpdt{i,1},1,1))./dz(i);
        end
    end
    
    subplot(sim_num,1,i);
    % Set the quantity to be plotted here
    plot_vec = dupdt{i,1};
    % Plots initial profile
%     plot(plot_vec(:,1),z_u{i,1}, 'DisplayName', 't=0.0');

    hold on
    % Plots profile at other times
    t1=100;
    nfig=5;
    for j=linspace(t1,t1+nfig-1,nfig)
        t_val = (1e-3*j)/t_scale(i);
        plot(plot_vec(1:end-1,j),z_u{i,1}(1:end-1),'DisplayName',['t=' num2str(t_val,4)]);
    end
    
    xlabel('dp_e/dt');
    ylabel("z/h");
    title("d = "+num2str(d(i), '%.2e'));   
    hold off  
end
legend();
% legend("0s","0.45s","0.9s","1.35s","1.8s","2.25s","2.7s","3.15s","3.6s","4.05s","4.5s", 'Location',"east");
% saveas(gcf,"dpdt_evolution.png")



make_animation = true;
ani_num = 1;

if (make_animation)
    nr_fr = 4000;
    t_vals = (1:nr_fr)*0.1/t_scale(ani_num);
    figure(1)
    animation_name = 'dil9_5_longt_small.gif';
    for i = 1 : nr_fr
%         plot(dilatancy{ani_num,1}(:,i), z_u{ani_num,1});
% %         xlim([0 2e-4]);
%         ylim([0 1]);
%         xlabel("Excess Pressure Diffusion");
%         ylabel("z");
%         legend("t="+num2str(t_vals(i)));
%         xlim([-0.1 0]);
%         ylim([0 1]);
%         xlabel("p_e");
%         ylabel("z");
        subplot(2,2,1)
        plot(p_e{ani_num,1}(:,i)./p_scale(ani_num), z_u{ani_num,1});
        legend("t="+num2str(t_vals(i)));
        xlim([-0.003 0]);
        ylim([0 1]);
        xlabel("p_e");
        ylabel("z");

        subplot(2,2,2)
        plot(phi{ani_num,1}(:,i), z_u{ani_num,1});
        xlim([-0.02,0])
        ylim([0 1])
        xlabel('\phi-\phi_c');
        ylabel("z");

        subplot(2,2,3)
        plot(u_f{ani_num,1}(:,i)/v_scale(ani_num), z_u{ani_num,1});
        xlim([0,0.5])
        ylim([0 1])
        xlabel("u_f");
        ylabel("z");

        subplot(2,2,4)
        plot(u_p{ani_num,1}(:,i)/v_scale(ani_num), z_u{ani_num,1});
        xlim([0,0.5])
        ylim([0 1])
        xlabel("u_p");
        ylabel("z");
        drawnow
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1
          imwrite(imind,cm,animation_name,'gif', 'LoopCount',Inf);
        else
          imwrite(imind,cm,animation_name,'gif','WriteMode','append','DelayTime', 0.02);
        end
    end
end

function beta_val=dimless_beta_fn(phihat)
    phi_c=0.6;
    beta_val = (phi_c + phihat).^2./((1-phi_c-phihat).^3);
end

function mu_val = mu_I_fn(I)
    mu1_I=0.342; % \frac{u+u_d}{\rho_f \phi_c}
    mu2_I=0.557; % 
    I_0 = 0.069;
    mu_val = tanh(1e10*I).*(mu1_I+(mu2_I-mu1_I)./(1+I_0./abs(I)));
end