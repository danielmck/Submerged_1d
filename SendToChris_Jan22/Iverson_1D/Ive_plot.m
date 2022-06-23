names = ["Ive_comp_13_deep.txt"]; % List of file names

% Loads the data to compare and puts them in a list of simulations
sim_num = size(names,2);
sim_list = cell(1,sim_num);
n_times = zeros(sim_num,1);

sim_type = strings(1,sim_num);
% Need to specify the type of simulation:
% dilatancy - "dil"
% diffusion only - "pdriv"
% constant pressure profile - "pcon"
% constant velocity profile - "ucon"

% Need to create lists storing the parameters in positions that corrrelate to
% the data in sim_list
N=zeros(1,sim_num);
h = zeros(1,sim_num);
d=zeros(1,sim_num);
dz = zeros(1,sim_num);
phi_c=zeros(1,sim_num);
eta_f_dl = zeros(1,sim_num);
g=9.81*ones(1,sim_num); % m/s^2
theta = zeros(1,sim_num); % deg
alpha_dl = zeros(1,sim_num); % 1/Pa
density_ratio = zeros(1,sim_num);

record = readtable('Results/result_record.csv');

% Extracts data from result_record.csv
for k=1:sim_num
    in_table = strcmp(record.Name, names(k));
    cd Results
    sim_list{k,1} = load(names(k));
    cd ../
    sim_type(k) = record.sim_type(in_table);
    N(k) = record.N(in_table);
    h(k) = record.h(in_table);
    d(k) = record.d(in_table);
    dz(k) = h/(N-0.5);
    phi_c(k) = record.phi_c(in_table);
    density_ratio(k) = record.rho_r(in_table);
    eta_f_dl(k) = record.eta_f(in_table);
    theta(k) = record.theta(in_table);
    alpha_dl(k) = record.alpha(in_table);
    n_times(k) = size(sim_list{k,1},1);
end

v_scale = sqrt(g.*h);
t_scale = sqrt(h./g);
z_scale = h;
rho = density_ratio.*phi_c+1-phi_c;

d_dl = d./z_scale;
dz_dl = dz./z_scale;

% Creates cells to store the matrices of values

buoyancy = zeros(sim_num);
t_vals = cell(sim_num,1);
z_pe = cell(sim_num,1);
z_u = cell(sim_num,1);
p_b = cell(sim_num,1);
p_e = cell(sim_num,1);
p_p = cell(sim_num,1);
phi = cell(sim_num,1);
u_f = cell(sim_num,1);
u_p = cell(sim_num,1);
w_f = cell(sim_num,1);
w_p = cell(sim_num,1);
I = cell(sim_num,1);
Iv = cell(sim_num,1);
phi_Iv = cell(sim_num,1);
mu_I = cell(sim_num,1);
mu_Iv = cell(sim_num,1);
beta_pe = cell(sim_num,1);
beta_u = cell(sim_num,1);
tau_p_I = cell(sim_num,1);
tau_p_Iv = cell(sim_num,1);
tau_f = cell(sim_num,1);
drag_mult_p = cell(sim_num,1);
drag_term_p = cell(sim_num,1);
drag_mult_f = cell(sim_num,1);
drag_term_f = cell(sim_num,1);
dpdz = cell(sim_num,1);
dppdz = cell(sim_num,1);
d2pdz2 = cell(sim_num,1);
d3pdz3 = cell(sim_num,1);
dufdz = cell(sim_num,1);
d2ufdz2 = cell(sim_num,1);
d2updz2 = cell(sim_num,1);
dupdz = cell(sim_num,1);
dilatancy = cell(sim_num,1);
diffusion_term = cell(sim_num,1);
tan_psi = cell(sim_num,1);
dpdt = cell(sim_num,1);
dphidt = cell(sim_num,1);
dphidz = cell(sim_num,1);
dupdt_I = cell(sim_num,1);
dupdt_Iv = cell(sim_num,1);
dufdt = cell(sim_num,1);
dudt = cell(sim_num,1);

% loops over the simulations and recreates all of the quantities used in
% the simulation, the names should be fairly self explanatory
for i = 1:sim_num
    buoyancy(i) = -(density_ratio(i)-1)*phi_c(i)*cosd(theta(i));
    z_pe{i,1} = linspace(1/(2*N(i)),1,N(i))';
    z_u{i,1} = linspace(0,1-1/(2*N(i)),N(i))';
    p_b{i,1} = (density_ratio(i)-1)*phi_c(i)*cosd(theta(i))*(1-z_pe{i,1});
    
    vec = sim_list{i,1};
    t_vals{i,1} = vec(:,1);
    vec = vec(:,2:end);
    % Processes the data differently depending on the simulation type
    if (sim_type(i) == "dil")
        p_e{i,1} = vec(:,1:N(i))';
        p_p{i,1} = p_b{i,1}-p_e{i,1};
        phi{i,1} = phi_c(i)+vec(:,N(i)+1:2*N(i))';
        u_f{i,1} = vec(:,2*N(i)+1:3*N(i))';
        u_p{i,1} = vec(:,3*N(i)+1:end)';
    elseif (sim_type(i) == "pcon")
        p_e{i,1} = s_frac(i)*p_b{i,1};
        p_p{i,1} = (1-s_frac(i))*p_b{i,1};
        % If looking at initial values for dilatancy sim, need to set phi
        % to phi_c here
        phi{i,1} = phi_c(i)+alpha_dl(i)*p_p{i,1};
        u_f{i,1} = vec(:,1:N(i))';
        u_p{i,1} = vec(:,N(i)+1:end)';
    elseif (sim_type(i) == "pdriv")
        p_e{i,1} = vec(:,1:N(i))';
        p_p{i,1} = (p_b{i,1}-p_e{i,1});
        phi{i,1} = phi_c(i)+alpha_dl(i)*p_p{i,1};
        u_f{i,1} = vec(:,N(i)+1:2*N(i))';
        u_p{i,1} = vec(:,2*N(i)+1:end)';
    elseif (sim_type(i) == "sin")
        p_e{i,1} = vec(:,1:N(i))';
        p_p{i,1} = (p_b{i,1}-p_e{i,1});
        phi{i,1} = phi_c(i)+vec(:,N(i)+1:2*N(i))';
        u_f{i,1} = vec(:,2*N(i)+1:end)';
        u_p{i,1} = vec(:,2*N(i)+1:end)';
    elseif (sim_type(i) == "ucon")
        p_e{i,1} = vec(:,1:N(i))';
        p_p{i,1} = (p_b{i,1}-p_e{i,1});
        phi{i,1} = phi_c(i)+vec(:,N(i)+1:end)';
        initial_u = load("no_phi_no_p.txt");
        u_f{i,1} = initial_u(7,1:N(i))'.*(ones(1,11));
        u_p{i,1} = initial_u(7,N(i)+1:end)'.*(ones(1,11));
    end
    
    if (sim_type(i)~="ucon")
        beta_pe{i,1} = 150*(phi{i,1}).^2.*eta_f_dl(i)./((1-phi{i,1}).^3.*d_dl(i)^2);
        beta_u{i,1} = interp1(z_pe{i,1},beta_pe{i,1},z_u{i,1},'linear','extrap');
        % If the velocity evolves we need to define these quantities
        dpdz{i,1} = vertcat(zeros(1,n_times(i)),diff(p_e{i,1}.*ones(N(i),n_times(i)),1,1))./dz_dl(i);
        d2pdz2{i,1} = vertcat(diff(dpdz{i,1},1,1),zeros(1,n_times(i)))./dz_dl(i);
        d3pdz3{i,1} = vertcat(zeros(1,n_times(i)),diff(d2pdz2{i,1},1,1))./dz_dl(i);
        dppdz{i,1} = vertcat(zeros(1,n_times(i)),diff(p_p{i,1}.*ones(N(i),n_times(i)),1,1))./dz_dl(i);
        dufdz{i,1} = vertcat(diff(u_f{i,1},1,1),zeros(1,n_times(i)))./dz_dl(i);
        d2ufdz2{i,1} = vertcat(zeros(1,n_times(i)),diff(dufdz{i,1},1,1))./dz_dl(i);
        dupdz{i,1} = vertcat(diff(u_p{i,1},1,1),zeros(1,n_times(i)))./dz_dl(i);
        d2updz2{i,1} = vertcat(zeros(1,n_times(i)),diff(dupdz{i,1},1,1))./dz_dl(i);
        Iv_temp = eta_f_dl(i).*abs(dupdz{i,1})./(p_p{i,1}+1e-8);
        Iv_temp(N(i),:) = eta_f_dl(i).*d2updz2{i,1}(N(i),:)./(-buoyancy(i)-dpdz{i,1}(N(i),:));
        Iv{i,1} = Iv_temp;
        w_p{i,1} = dpdz{i,1}./beta_u{i,1};
        w_f{i,1} = -phi_c(i)/(1-phi_c(i)).*w_p{i,1};
        tau_f{i,1} = eta_f_dl(i).*d2ufdz2{i,1}./(1-phi_c(i));
    
        I{i,1} = 2.*d_dl(i).*abs(dupdz{i,1}).*sqrt(density_ratio(i))./sqrt(p_p{i,1}.*ones(1,n_times(i))+1e-8);
        mu_I{i,1} = mu_I_fn(I{i,1});
        mu_Iv{i,1} = mu_Iv_fn(Iv{i,1});
        tau_p_I{i,1} = 1/(phi_c(i).*density_ratio(i)).*sign(dupdz{i,1}).*vertcat(zeros(1,n_times(i)),diff(mu_I{i,1}.*(p_p{i,1}.*ones(1,n_times(i))),1,1)./dz_dl(i));
        tau_p_Iv{i,1} = 1/(phi_c(i).*density_ratio(i)).*sign(dupdz{i,1}).*vertcat(zeros(1,n_times(i)),diff(mu_Iv{i,1}.*(p_p{i,1}.*ones(1,n_times(i))),1,1)./dz_dl(i));
        if (sim_type(i)=="sin")
            dudt{i,1} = ((phi_c(i).*density_ratio(i)).*tau_p_I{i,1}+(1-phi_c(i)).*tau_f{i,1})/rho(i)+sind(theta(i));
        else
            drag_mult_p{i,1} = (1-phi_c(i))^2.*beta_u{i,1}./(density_ratio(i)*phi_c(i)).*ones(1,n_times(i));

            drag_term_p{i,1} = drag_mult_p{i,1}.*(u_f{i,1}-u_p{i,1});

            drag_mult_f{i,1} = (1-phi_c(i)).*beta_u{i,1}.*ones(1,n_times(i));

            drag_term_f{i,1} = drag_mult_f{i,1}.*(u_f{i,1}-u_p{i,1});

            dupdt_I{i,1} = vertcat(zeros(1,n_times(i)),tau_p_I{i,1}(2:end,:)+sind(theta(i))+drag_term_p{i,1}(2:end,:));
            dupdt_Iv{i,1} = vertcat(zeros(1,n_times(i)),tau_p_Iv{i,1}(2:end,:)+sind(theta(i))+drag_term_p{i,1}(2:end,:));
            dufdt{i,1} = vertcat(zeros(1,n_times(i)),tau_f{i,1}(2:end,:)+sind(theta(i))-drag_term_f{i,1}(2:end,:));
        end
    end
    
    
    % If the pressure evolves, have to define these quantities
    if ((sim_type(i)=="pdriv") || (sim_type(i)=="dil") || (sim_type(i)=="ucon") || (sim_type(i)=="sin"))
        dphidz{i,1} = vertcat(zeros(1,n_times(i)),diff(phi{i,1},1,1))./dz_dl(i);
        dpdt{i,1} = 1./(alpha_dl(i)).*vertcat(diff(1./beta_pe{i,1}.*dpdz{i,1}),zeros(1,n_times(i)))./dz_dl(i);
        if ((sim_type(i)=="dil") || (sim_type(i)=="ucon") || (sim_type(i)=="sin"))
            phi_Iv{i,1} = phi_c(i)./(1+sqrt(abs(Iv{i,1})));
            tan_psi{i,1} = phi{i,1}-phi_Iv{i,1};
            dilatancy{i,1} = -1./(alpha_dl(i)).*dupdz{i,1}.*(tan_psi{i,1});
            diffusion_term{i,1} = dpdt{i,1};
            dphidt{i,1} = -dpdt{i,1}.*alpha_dl(i).*phi_c(i);
            dpdt{i,1} = dpdt{i,1} + dilatancy{i,1};
            d2pdzdt = vertcat(zeros(1,n_times(i)),diff(dpdt{i,1},1,1))./dz_dl(i);
        end
    end
end
%% 

for i=1:sim_num

    mu1=0.342;
    mu2 = 0.557;
    I0 = 0.069;
    theta0 = 13;
 
    crit_grad = -(density_ratio(i)*phi_c(i)+(1-phi_c(i)))*sind(theta(i))/0.32;
    crit_pe = p_b{i,1}+crit_grad*(1-z_pe{i,1});
    
    tan_psi_ss = alpha_dl(i).*eta_f_dl(i).^2./(p_p{i,1}.*Iv{i,1}.^2).*vertcat(zeros(1,n_times(i)),diff(tau_p_Iv{i,1})./dz_dl(i));
    % Steady state dilatancy angle after fast timescale
    
    dudz_ss_grad = rho(i).*sind(theta(i))./mu_Iv{i,1}.*(1-z_pe{i,1});
    dudz_ss = rho(i).*sind(theta(i).*(phi{i,1}(100,1)-phi_c(i)).^2)./(phi{i,1}(100,1).^2.*eta_f_dl(i).*mu_Iv{i,1}(100,1)).*(1-z_pe{i,1});
    fric_mult = (phi{i,1}(100,1).^2.*eta_f_dl(i).*mu_Iv{i,1}(100,1))./(phi{i,1}(100,1)-phi_c(i)).^2;
    dudz_init = (rho(i)-1).*cosd(theta0).*(phi{i,1}(100,1)-phi_c(i)).^2./(phi{i,1}(100,1).^2.*eta_f_dl(i));
    cos_const_term = (phi{i,1}(100,1)-phi_c(i)).^2./(phi{i,1}(100,1).^2.*eta_f_dl(i)).*((rho(i)-1).*cosd(theta0)-rho(i).*sind(theta(i))./mu_Iv{i,1}(100,1));
    dudz_medium = dudz_ss;
    for k =1:10
        dudz_medium = dudz_medium + 8./(2*k-1).^2./pi.^2.*cos_const_term.*exp(-fric_mult(:,1)./rho(i).*pi^2/4.*(2*k-1).^2.*t_vals{i,1}').*cos(pi/2.*(2*k-1).*z_pe{i,1});
    end
    % The medium timescale approximation
    
%     SetPaperSize(10,10);
%     subplot(1,2,1);

% Alters colurs used
%     C = brewermap(nfig+2, 'PuRd');
%     set(0, 'DefaultAxesColorOrder', C(3:end,:))
%      cmap = colororder();
     set(0,'DefaultAxesColorOrder','remove')

%     plot_names = ["ICs from 9.0 degree slope","ICs from 9.1 degree slope","ICs from 9.2 degree slope","ICs from 13 degree slope"];
    nline=1;
    plot_vec = mu_Iv{i,1};
    t_start = 500;
    t_step = 500;
    for j=linspace(1,nline,nline)
%         t_vals{j,1}(t1(j))
        t_max = t_start+(j-1)*t_step;
        t1=[sum(t_vals{i,1}<t_max)+1];
        hold on
        plot(plot_vec(1:end,t1),z_pe{i,1}(1:end),'DisplayName',"t="+num2str(t_max))

        ylabel("$z/h$");
        xlabel('$p_e$');
        box on
    end
%     plot(tan_psi_ss(1:end-1,t1),z_pe{i,1}(1:end-1),'DisplayName',"Steady State")
%     plot(crit_pe,z_pe{i,1}(1:end),"k--",'DisplayName',"Critical $p_e$");
    
    legend('Location', "best",'UserData', 5);

    title('The evolution of $p_e$ showing the expansion of the boundary layer')
    pdfname = 'Ive_pe_evo_post_medium';    
    
%     PrintFig(pdfname)


% ///// Below is code allowing comparison to depth averaged and reduced depth
% averaged models ////
%     da_vals = load("../Iverson_DA/DA_Results/Ive_da_5_deep_alpha4.txt");
%     da_times = da_vals(:,1);
%     da_phi = da_vals(:,3);
%     da_u = da_vals(:,4);
%     da_pe = da_vals(:,5)-cosd(5);
%     da_Iv = 3.*eta_f_dl.*da_u./((rho-1).*cosd(theta)-da_pe);
%     da_phi_Iv = phi_c./(1+sqrt(da_Iv));
%     
%     red_vals = load("../Iverson_DA/DA_Results/Ive_da_red_5_deep_13_start.txt");
%     red_times = red_vals(:,1);
%     red_h = red_vals(:,2);
%     red_phi = red_vals(:,3);
%     red_u = red_vals(:,4);
%     red_pp = (red_phi.^2./(red_phi-phi_c).^2).*3.*red_u.*eta_f_dl./red_h;
% 
%     u_ave1 = depth_average(u_p{1,1},N(1),n_times);
%      SetPaperSize(12,10);


%      tmin = 400;
%      tmax = 15000;
%      plot(t_vals{1,1}(t_vals{1,1}<tmax),u_ave1(t_vals{1,1}<tmax),"blue",'DisplayName', "Full Model")
%      plot(da_times(da_times<tmax),da_u(da_times<tmax),"red",'DisplayName', "Depth Averaged Model")
%      plot(red_times(red_times<tmax),red_u(red_times<tmax),"green",'DisplayName', "Reduced Depth Averaged Model")
%      xlabel("$p_e$");

     
%      ylabel('$\Bar{u}$');
%      xlabel('$t$');
%      legend('Location', "best",'UserData', 21);
%      title("Late evolution of different models from $13$ degree slope ICs")
     
%     pdfname = 'Ive_da_vs_full_vs_red_u_comp_13';    
    % PrintFig(pdfname)
end

% ///// Really old code to make videos ////
% ani_num = 1;
% make_avi=false;
% if (make_avi) 
%     nframes=6000;
%     t_vals = ((0:nframes))*10;
%     Frames=moviein(nframes);
%     figure
%     v = VideoWriter('dil9_5_dilatant_depo.avi');
%     v.FrameRate=50;
%     open(v)
%     ax = gca();
%     for i=2000:6000
%         subplot(2,2,1)
%         plot(p_e{ani_num,1}(:,i), z_pe{ani_num,1});
%         hold on
%         plot((-buoyancy(ani_num)-(density_ratio(ani_num)*phi_c(ani_num)+(1-phi_c(ani_num)))*sind(theta(ani_num))/0.342)*(1-z_pe{ani_num,1}), z_pe{1,1},"--");
%         hold off
%         legend("t="+num2str(t_vals(i-1999),'%.1f'),"Critical $p_e$","Location",'northeast','Interpreter','latex');
% %         legend();
%         xlim([0 0.06]);
%         ylim([0 1]);
%         xlabel("$p_e$",'Interpreter','latex');
%         ylabel("$z$",'Interpreter','latex');
% 
%         subplot(2,2,2)
%         plot(phi_c(1)+phi{ani_num,1}(:,i), z_u{ani_num,1});
% %         hold on
% %         plot(phi{2,1}(:,i+29), z_u{2,1});
% %         hold off
%         xlim([0.55,0.6])
%         ylim([0 1])
%         xlabel('$\phi$','Interpreter','latex');
%         ylabel("$z$",'Interpreter','latex');
% 
%         subplot(2,2,3)
%         plot(u_f{ani_num,1}(:,i), z_u{ani_num,1});
% %         hold on
% %         plot(u_f{2,1}(:,i+29), z_u{2,1});
% %         hold off
%         xlim([0,0.5])
%         ylim([0 1])
%         xlabel("$u_f$",'Interpreter','latex');
%         ylabel("$z$",'Interpreter','latex');
% 
%         subplot(2,2,4)
%         plot(u_p{ani_num,1}(:,i), z_u{ani_num,1});
% %         hold on
% %         plot(u_p{2,1}(:,i+29), z_u{2,1});
% %         hold off
%         xlim([0,0.5])
%         ylim([0 1])
%         xlabel("$u_p$",'Interpreter','latex');
%         ylabel("$z$",'Interpreter','latex','Interpreter','latex');
%         drawnow();
%         writeVideo(v,getframe(gcf))
%     end
%     close all
%     close(v)
% end

function beta_val=dimless_beta_fn(phihat)
    phi_c=0.6;
    beta_val = (phi_c + phihat).^2./((1-phi_c-phihat).^3);
end

function mu_val = mu_I_fn(I)
    mu1_I=0.342; % \frac{u+u_d}{\rho_f \phi_c}
    mu2_I=0.557; % 
    I_0 = 0.069;
    mu_val = tanh(1e8*I).*(mu1_I+(mu2_I-mu1_I)./(1+I_0./abs(I)));
end

function mu_val = mu_Iv_fn(Iv)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    phi_c=0.6;
    mu_val = tanh(1e8*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
end