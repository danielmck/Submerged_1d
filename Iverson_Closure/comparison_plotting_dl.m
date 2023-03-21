% opengl firmware
names = ["Ive_5deg_13init.txt"]; %"Ive_comp_4_deep_1_flux.txt"];%,"Ive_comp_4_deep_9_2_start.txt""Ive_comp_4_deep_9_1_start.txt","
% names = ["Ive_comp_5_deep_custom_time_IC.txt"]; %,"Ive_comp_5_deep_custom_time_IC.txt"];%
% Loads the data to compare and puts them in a list of simulations

sim_num = size(names,2);
sim_list = cell(1,sim_num);
n_times = zeros(sim_num,1);
% for j=1:sim_num
%     n_times(1,j) = size(sim_list{1,j},1);
% end
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
s_frac = [0.6,0.6,0.6];
density_ratio = zeros(1,sim_num);
t_step = zeros(1,sim_num);

record = readtable('Results/result_record.csv');

for k=1:sim_num
    in_table = strcmp(record.Name, names(k));
    cd Results
    sim_list{k,1} = load(names(k));
    cd ../
    sim_type(k) = record.sim_type(in_table);
    N(k) = record.N(in_table);
    h(k) = record.h(in_table);
    d(k) = record.d(in_table);
    dz(k) = h(k)/(N(k)-0.5);
    phi_c(k) = record.phi_c(in_table);
    density_ratio(k) = record.rho_r(in_table);
    eta_f_dl(k) = record.eta_f(in_table);
    theta(k) = record.theta(in_table);
    alpha_dl(k) = record.alpha(in_table);
    t_step(k) = record.t_step(in_table);
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
Iv_phi = cell(sim_num,1);
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
d2updtdz = cell(sim_num,1);
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
flow_depth = cell(sim_num,1);

% loops over the simulations
for i = 1:sim_num
    buoyancy(i) = -(density_ratio(i)-1)*phi_c(i)*cosd(theta(i));
    z_pe{i,1} = linspace(1/(2*N(i)),1,N(i))';
    z_u{i,1} = linspace(0,1-1/(2*N(i)),N(i))';
    p_b{i,1} = (density_ratio(i)-1)*phi_c(i)*cosd(theta(i))*(1-z_pe{i,1});
    
    vec = sim_list{i,1};
    if (mod(size(vec,2),N(i)))
        t_vals{i,1} = vec(:,1);
        vec = vec(:,2:end);
    else
        t_vals{i,1} = linspace(0,(n_times(i)-1)*t_step(i),n_times(i));
    end
%     tlength=[sum(t_vals{i,1}<1000)];
%     t_vals{i,1} = t_vals{i,1}(1:tlength+1);
%     vec = vec(1:tlength+1,:);
%     n_times(i) = tlength+1;
    
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
        % to zero here
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
    
    flow_depth{i,1} = sum(u_p{i,1} > 1e-3,1)/N(i);
    
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
            d2updtdz{i,1} = vertcat(diff(dupdt_Iv{i,1},1,1)./dz_dl(i), zeros(1,n_times(i)));
        end
    end
    
    
    % If the pressure evolves, have to define these quantities
    if ((sim_type(i)=="pdriv") || (sim_type(i)=="dil") || (sim_type(i)=="ucon") || (sim_type(i)=="sin"))
        dphidz{i,1} = vertcat(zeros(1,n_times(i)),diff(phi{i,1},1,1))./dz_dl(i);
        dpdt{i,1} = 1./(alpha_dl(i)).*vertcat(diff(1./beta_pe{i,1}.*dpdz{i,1}),zeros(1,n_times(i)))./dz_dl(i);
        if ((sim_type(i)=="dil") || (sim_type(i)=="ucon") || (sim_type(i)=="sin"))
            phi_Iv{i,1} = phi_c(i)./(1+sqrt(abs(Iv{i,1})));
            Iv_phi{i,1} = (phi{i,1}-phi_c(i)).^2./phi{i,1}.^2;
            tan_psi{i,1} = phi{i,1}-phi_Iv{i,1};
            dilatancy{i,1} = -1./(alpha_dl(i)).*dupdz{i,1}.*(phi{i,1}-phi_Iv{i,1});
            diffusion_term{i,1} = dpdt{i,1};
            dphidt{i,1} = -dpdt{i,1}.*alpha_dl(i).*phi_c(i);
            dpdt{i,1} = dpdt{i,1} + dilatancy{i,1};
            d2pdzdt = vertcat(zeros(1,n_times(i)),diff(dpdt{i,1},1,1))./dz_dl(i);
        end
    end
end
%% 

for i=1:sim_num

    crit_grad = -(density_ratio(i)*phi_c(i)+(1-phi_c(i)))*sind(theta(i))/0.32;
    crit_pe = p_b{i,1}+crit_grad*(1-z_pe{i,1});
    

    dIvdt_u = vertcat(diff(tau_p_Iv{i,1},1,1)./dz_dl(i),zeros(1,n_times(i))).*eta_f_dl(i)./p_p{i,1};
    dIvdt_pp = 1./p_p{i,1}.*Iv{i,1}.*dpdt{i,1};
    dIvdt = dIvdt_u + dIvdt_pp;
    dphidt_Iv = phi{i,1}.^2./(phi_c(i).*sqrt(Iv{i,1})).*dIvdt;

    theta0 = 13;
    tan_psi_ss = alpha_dl(i).*eta_f_dl(i).^2./(p_p{i,1}.*Iv{i,1}.^2).*vertcat(zeros(1,n_times(i)),diff(tau_p_Iv{i,1})./dz_dl(i));
    t_max = 70;
    
    dudz_ss_grad = rho(i).*sind(theta(i))./mu_Iv{i,1}.*(1-z_pe{i,1});
    dudz_ss = rho(i).*sind(theta(i).*(phi{i,1}(100,1)-phi_c(i)).^2)./(phi{i,1}(100,1).^2.*eta_f_dl(i).*mu_Iv{i,1}(100,1)).*(1-z_pe{i,1});
    fric_mult = (phi{i,1}(100,1).^2.*eta_f_dl(i).*mu_Iv{i,1}(100,1))./(phi{i,1}(100,1)-phi_c(i)).^2;
    dudz_init = (rho(i)-1).*cosd(theta0).*(phi{i,1}(100,1)-phi_c(i)).^2./(phi{i,1}(100,1).^2.*eta_f_dl(i));
    cos_const_term = (phi{i,1}(100,1)-phi_c(i)).^2./(phi{i,1}(100,1).^2.*eta_f_dl(i)).*((rho(i)-1).*cosd(theta0)-rho(i).*sind(theta(i))./mu_Iv{i,1}(100,1));
%     dudz_medium = dudz_ss;
%     for k =1:10
%         dudz_medium = dudz_medium + 8./(2*k-1).^2./pi.^2.*cos_const_term.*exp(-fric_mult(:,1)./rho(i).*pi^2/4.*(2*k-1).^2.*t_vals{i,1}').*cos(pi/2.*(2*k-1).*z_pe{i,1});
%     end
    %,sum(t_vals{2,1}<t_max),sum(t_vals{3,1}<t_max),5001];
%     pp_medium = dudz_medium.*eta_f_dl.*phi{i,1}(1,100)^2./(phi{i,1}(1,100)-phi_c(i))^2;
%     pe_medium = p_b{i,1} - pp_medium;
    

%     SetPaperSize(10,10);
    if (i==1)
        plot_name = "Short Run";
        colour = 'Purples';
    elseif (i==2)
        colour = 'Reds';
    else
        plot_name = "Long Run";
        colour = 'Greens';
    end
    
%     plot(linspace(0,2e-4,201),mu_Iv_fn(linspace(0,0.02,201)))
    nfig=1;
    C = brewermap(nfig+2, colour);
    plot_vec = u_p{i,1};
    da_Iv = depth_average(Iv{i,1},N(i),n_times(i));
    da_phi_Iv = depth_average(phi_Iv{i,1},N(i),n_times(i));
    for j=linspace(1,nfig,nfig)
%         t_vals{j,1}(t1(j))
        t_max = 1800+150*(j-1);
%         t_max = [50000,100000];
%         t_val = (k-1)*t_step(i);
        t1=[sum(t_vals{i,1}<t_max)+1];
        hold on
        C = brewermap(nfig+2, colour);
%         plot(plot_vec(1:end,t1),z_pe{i,1}(1:end),'DisplayName','Full Profile')%,'color',C(j+2,:)"t="+num2str(t_max))
%         plot(plot_vec(end,t1)*(1-(1-z_u{i,1}).^2),z_pe{i,1}(1:end),'color','green','DisplayName',"t="+num2str(t_max))
%         plot(da_phi_Iv(t1)*ones(N(i),1),z_pe{i,1}(1:end),'DisplayName','Depth Averaged Value')
%         plot((phi_c/(1+sqrt(da_Iv(t1))))*ones(N(i),1),z_pe{i,1}(1:end),'DisplayName','Value with DA $I_v$')
        
%         depth_average(plot_vec(1:end,t1),N(i),1)
        ylabel("$z/h$");
        xlabel('$\phi(I_v)$');
%         xlabel('$\frac{\partial u_p}{\partial z}$');
        box on
%         C = brewermap(nfig+2, 'Reds');
%         xlim([0,25])
%         plot(diffusion_term{i,1}(1:end,j),z_pe{i,1}(1:end),'DisplayName',"Diffusion Term")
%         plot(pe_medium(:,t1),z_pe{i,1},'color',C(j+2,:),'DisplayName',"Approx"+" t="+num2str(t_max))
    end
    
    eps_bl = 0.4;
%     plot(tand(13).*(density_ratio*phi_c(i)+1-phi_c(i)).*ones(200,1),z_pe{i,1}(1:end),'DisplayName',"Steady State")
%     plot(sind(theta)*(density_ratio-1)*phi_c(i)/tand(13).*(1-z_pe{i,1}),z_pe{i,1}(1:end),'DisplayName',"Steady State")
%     const_bl = ((density_ratio(i)-1)*phi_c(i)*cosd(theta(i))-sind(theta)*(density_ratio-1)*phi_c(i)/tand(13))*(1-eps_bl);
%     cos_coeff_bl = ((density_ratio(i)-1)*phi_c(i)*cosd(theta(i))-sind(theta)*(density_ratio-1)*phi_c(i)/tand(13))*2*eps_bl/pi;
%     plot(const_bl+(cos_coeff_bl).*(1-(z_pe{i,1}/0.4).^2),z_pe{i,1}(1:end),'DisplayName',"Steady State")
%     plot(const_bl+(cos_coeff_bl).*(cos(z_pe{i,1}/eps_bl*pi/2)),z_pe{i,1}(1:end),'DisplayName',"Steady State")
%     plot(crit_pe,z_pe{i,1}(1:end),"k--",'DisplayName',"Critical $p_e$");
    legend('Location', "best",'UserData', 14);
%     ax = gca;
%     annotation('arrow',[0.306701030927835 0.268041237113402],...
%     [0.608164420485175 0.746630727762803]);
    title('Values of $\phi(I_v)$ with $\mu(I)$ Rheology') 
%     cd Figures/ShallowSlowing
    pdfname = 'Ive_mu_I_DA_phi_Iv';
%     PrintFig(pdfname)
%     cd ../../



%     plot((density_ratio(i)-1)*phi_c(i)*sind(theta(i))/tand(13)*(1-z_pe{i,1}),z_pe{i,1}(1:end),'DisplayName',"Analytic Min Value")
%     plot(crit_pe,z_pe{i,1},"k--",'DisplayName',"Critical $p_e$")
%     legend('Location', "best",'UserData', 12);
%     ylabel("z/h");
%     xlabel('$\frac{\partial p_e}{\partial t}$');
% % % 
% % % %     
%     da_vals = load("../Iverson_DA/Results/Ive_da_5_deep_alpha4.txt");
%     da_times = da_vals(:,1);
%     da_phi = da_vals(:,3);
%     da_u = da_vals(:,4);
%     da_pe = da_vals(:,5)-cosd(5);
%     da_Iv = 3.*eta_f_dl.*da_u./((rho-1).*cosd(theta)-da_pe);
%     da_phi_Iv = phi_c./(1+sqrt(da_Iv));
% %     plot(da_u(sum(da_times<300)),0,"rx",'DisplayName',"DA Model Basal Value")
%     
%     red_vals = load("../Iverson_DA/Results/Ive_da_red_5_deep_13_start.txt");
%     red_times = red_vals(:,1);
%     red_h = red_vals(:,2);
%     red_phi = red_vals(:,3);
%     red_u = red_vals(:,4);
%     % red_pe = red_vals(:,5)-cosd(theta);
%     red_pp = (red_phi.^2./(red_phi-phi_c).^2).*3.*red_u.*eta_f_dl./red_h;

%     u_ave1 = depth_average(u_p{1,1},N(1),n_times);
%      SetPaperSize(12,10);

%      [peak_pos, peak_ind] = min(dilatancy{i,1},[],1);
%      yyaxis left;
     tmin = 10000;
     tmax = 50000;
%      [~,da_peak] = max(da_pe);
%      [~,full_peak] = max(p_e{i,1}(1,:));
     
%      plot(t_vals{1,1}(max(sum(t_vals{1,1}<tmin),1):sum(t_vals{1,1}<tmax)),u_ave1(max(sum(t_vals{1,1}<tmin),1):sum(t_vals{1,1}<tmax)),"blue",'DisplayName', "Full Model")
     
%      tmin = 400;
%      tmax = 15000;
%      plot(da_u(min(sum(da_times<tmin),1):sum(da_times<tmax)),da_times(min(sum(da_times<tmin),1):sum(da_times<tmax)),"red",'DisplayName', "Depth Averaged Model")
%      plot(t_vals{1,1}(t_vals{1,1}<tmax),u_ave1(t_vals{1,1}<tmax),"blue",'DisplayName', "Full Model")
%      plot(da_times(da_times<tmax),da_u(da_times<tmax),"red",'DisplayName', "Depth Averaged Model")
%      plot(red_times(red_times<tmax),red_u(red_times<tmax),"green",'DisplayName', "Reduced Depth Averaged Model")
%      xlabel("$p_e$");
%      yyaxis right;
%      plot(linspace(1,5000,5000+1),p_e{i,1}(1,:),'DisplayName', "Peak Magnitude")
%      ylabel("Basal $p_e$");
     
%      set(gca,'ytick',[])
%      ax.YAxis(1).Exponent = 0;
%      ax.XAxis.Exponent = 0;
%      xlim([10000 20000]);
%      ylabel('$\Bar{u}$');
%      xlabel('$t$');
%      legend('Location', "best",'UserData', 21);
%      title("Late evolution of different models from $13$ degree slope ICs")
     
%      xtickfomat('%5.1e');
% %     hold off
%     pdfname = 'Ive_da_vs_full_vs_red_u_comp_13';    
    
%     movefile([pdfname '.pdf'], 'Figures/FirstYearReport/')
end
% PrintFig(pdfname)
% movefile([pdfname '.pdf'], 'Figures/FirstYearReport/')
% 
% SetPaperSize(10,10);
% da_ave = load("Results/dil_9deg_depth_averaged_slowing_unmod_shear.txt");
% u_ave1 = depth_average(u_p{1,1},N(1),n_times);
% plot(linspace(1,5001,5001),u_ave1,'DisplayName',"Full Model")
% phi_only_ave = depth_average(phi_only',N(1),3500);
% plot(linspace(1500,5000,3500),phi_c(i)+phi_only_ave,'DisplayName',"Full Model")
% save("slowing_depth_ave_phi.txt", 'u_ave1','-ascii')
% hold on
% plot(linspace(1,501,501),da_ave,'DisplayName',"Depth Averaged Model")
% xlabel("Time")
% ylabel("Depth averaged velocity $\bar{u}$")
% title("Comparison of Slowing Velocity between Different Models")
% legend('Location', "northeast",'UserData', 15)
% PrintFig('dil_9deg_depth_ave_u_comp')
% movefile('dil_9deg_depth_ave_u_comp.pdf', 'Figures/Dilatant/')

% SetPaperSize(10,10);
% ave1 = depth_average(u_p{1,1},N(1),3000);
% ave2 = depth_average(u_p{2,1},N(2),3000);
% semilogx(10*linspace(1,3000,3000),ave1,'DisplayName','Standard Case');
% hold on
% semilogx(10*linspace(1,3000,3000),ave2,'DisplayName','Large $\alpha$ Case');
% xlabel('$t$');
% ylabel('$\bar{p_e}$');
% title('Relative $\bar{p_e}$ Evolution of Large $\alpha$ Case');
% legend('Location', "northeast",'UserData', 8);
% ax = gca;
% ax.XAxis.Exponent = 0;
% ax.YAxis.Exponent = 0;
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',8)
% a = get(gca,'YTickLabel');
% set(gca,'YTickLabel',a,'fontsize',8)
% ytickformat('%5.1e');
% % xtickformat('%5.1e');
% pdfname = 'dil9deg_assymp_large_alpha_up';
% PrintFig(pdfname)

% The command from Chris to clean up the plots
% ylabel('stuff','Position',[x y]) 

ani_num = 1;
make_avi=true;
C = viridis(4);
if (make_avi) 
    nframes=100;
%     t_vals = ((0:nframes))*10;
%     Frames=moviein(nframes);
%     figure
%     v = VideoWriter('Ive_small_flux_slowing_4deg.avi');
%     v.FrameRate=50;
%     open(v)
    ax = gca();
    workingDir = 'fluxvids';
%     mkdir(workingDir)
%     mkdir(workingDir,'images')
    cd fluxvids/images
    width = 12;
    height = 7;
    width_p = ceil(width*37.7952755906);
    height_p = ceil(height*37.7952755906);
    f = figure;
%     set(f,'Position',[15 15 1200 700])
    for i=1:1
        clf;
        % Set a size if desired
        
%         SetPaperSize(12,7);
        t_max = 0.0+(i-1)*5;
        
        subplot(1,2,2)
        for j=1:sim_num
            t_index = sum(t_vals{j,1}<t_max)+1;
            plot(p_e{j,1}(:,t_index), z_pe{ani_num,1},'color',C(1,:));
            if (j==1) 
                hold on 
            end
        end
%         plot(crit_pe, z_pe{1,1},"--");
        hold off
%         legend("Initial Flux = 0.6","Initial Flux = 0.7","Initial Flux = 0.8","Critical $p_e$","Location",'northeast','Interpreter','latex');
%         legend();
%         text(0.225,0.9, {"t="+num2str(t_max,'%11.1f'),"Rate = $\times 10$"},'Interpreter','latex')
%         annotation('textbox', [0.9, 0.7, 0.15, 0.1], 'String', "t="+num2str(t_max))
        xlim([0 0.5]);
        ylim([0 1]);
        xlabel("$p_e$",'Interpreter','latex','FontSize',12);
%         ylabel("$z$",'Interpreter','latex');

%         subplot(2,2,2)
%         plot(phi_c(1)+phi{ani_num,1}(:,i), z_u{ani_num,1});
%         hold on
%         plot(phi{2,1}(:,i+29), z_u{2,1});
%         hold off
%         xlim([0.55,0.6])
%         ylim([0 1])
%         xlabel('$\phi$','Interpreter','latex');
%         ylabel("$z$",'Interpreter','latex');

        subplot(1,2,1)
        col = [0.9 0.1 0.1];
        for j=1:ani_num
            t_index = sum(t_vals{j,1}<t_max)+1;
            plot(u_p{j,1}(:,t_index), z_pe{ani_num,1},'color',C(1,:));
            if (j==1) 
                hold on
                col = [0.1 0.1 0.9];
            end
        end
        xlim([0 15])
        ylim([0 1])
        xlabel("$u$",'Interpreter','latex','FontSize',12);
        ylabel("$z$",'Interpreter','latex','FontSize',12);
        sgtitle("t="+num2str(t_max))      
%         subplot(2,1,2)
%         col = [0.9 0.1 0.1];
%         for j=1:sim_num
%             t_index = sum(t_vals{j,1}<t_max)+1;
%             plot(t_vals{j,1}(1:t_index), flow_depth{j,1}(1:t_index),'LineWidth',2, 'Color', col);
%             if (j==1) 
%                 hold on 
%             end
%             plot(t_vals{j,1}(t_index), flow_depth{j,1}(t_index),'.','MarkerSize',15, 'Color', col)
%             if (j==1) 
% %                 hold on
%                 col = [0.1 0.1 0.9];
%             end
%             xlim([0 300])
%             ylim([0 1])
%         end
%         xlabel("$t$",'Interpreter','latex','FontSize',16);
%         ylabel("Flow Height",'Interpreter','latex','FontSize',16);
%         hold on
%         plot(u_f{2,1}(:,i+29), z_u{2,1});
        hold off
%         xlim([0,1.2])
%         ylim([0 1])
        

%         subplot(2,2,4)
%         plot(u_p{ani_num,1}(:,i), z_u{ani_num,1});
% %         hold on
% %         plot(u_p{2,1}(:,i+29), z_u{2,1});
% %         hold off
%         xlim([0,0.5])
%         ylim([0 1])
%         xlabel("$u_p$",'Interpreter','latex');
%         ylabel("$z$",'Interpreter','latex','Interpreter','latex');
%         sgtitle("Slowing of Small Flux Flows on 4 degree Slope",'Interpreter','latex')
%         drawnow();
%         writeVideo(v,getframe(gcf))
%         set(f,'Position',[15 15 1200 700])
%         exportgraphics(f,"img"+num2str(i,'%04d')+".jpg","Resolution",300)
        figsave(f,"img"+num2str(i,'%04d')+".jpg",[1200 700])
        
    end
    close(f)
%     close all
%     close(v)
    cd ../../
end

%%
imageNames = dir(fullfile(workingDir,'images','*.jpg'));
imageNames = {imageNames.name}';

outputVideo = VideoWriter(fullfile(workingDir,'shuttle_out.avi'));
outputVideo.FrameRate = 30;
open(outputVideo)
for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir,'images',imageNames{ii}));
   writeVideo(outputVideo,img)
end
close(outputVideo)

function beta_val=dimless_beta_fn(phi)
    phi_c=0.6;
    beta_val = (phi).^2./((1-phi).^3);
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
    mu_val = tanh(1e8*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(abs(Iv)));
end

function figsave(fig,file,rez,txt,bak)
    %Save figure as image, with custom resolutions and text scaling.
    % figsave(fig,file,rez,txt,bak)
    %
    %Example:
    % clf,text(0.1,0.5,{'This text should be';'50 pixels high';'and the image';'900W x 600H pix'},'FontSize',50)
    % figsave(gcf,'Figure.jpg',[900 600])
    if nargin<1 || isempty(fig),  fig  = gcf; end          %figure handle
    if nargin<2 || isempty(file), file = 'Figure.jpg'; end %output file name
    if nargin<3 || isempty(rez),  rez  = [900 600]; end    %resolution [horizontal vertical]
    if nargin<4 || isempty(txt),  txt  = 1; end            %text scale factor
    if nargin<5 || isempty(bak),  bak  = 1; end            %preserve background colour
    set(fig,'PaperPosition',[0 0 rez/(100*txt)],'PaperUnits','inches'); %set paper size (does not affect display)
    if bak
        set(fig,'InvertHardcopy','off'); %preserve background colour
    end
    imwrite(print(gcf,'-RGBImage',['-r' num2str(100*txt,'%f')]),file) %print RGB image to file (slow)
end