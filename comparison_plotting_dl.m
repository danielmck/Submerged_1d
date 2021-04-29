opengl software
names = ["dil9_5deg_init_elastic.txt","dil9deg_slowing_flow.txt"];

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
s_frac = [0.0,0,0.7];
density_ratio = zeros(1,sim_num);
t_step = zeros(1,sim_num);

record = readtable('Results/result_record.csv');

for k=1:sim_num
    in_table = strcmp(record.Name, names(k));
    cd Results
    sim_list{k,1} = load(names(k));
    cd ../
    sim_type(k) = record.Type(in_table);
    N(k) = record.N(in_table);
    h(k) = record.h(in_table);
    d(k) = record.d(in_table);
    dz(i) = h/(N-0.5);
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


d_dl = d./z_scale;
dz_dl = dz./z_scale;

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
dppdz = cell(sim_num,1);
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
    buoyancy(i) = -(density_ratio(i)-1)*phi_c(i)*cosd(theta(i));
    z_pe{i,1} = linspace(1/(2*N(i)),1,N(i))';
    z_u{i,1} = linspace(0,1-1/(2*N(i)),N(i))';
    p_b{i,1} = (density_ratio(i)-1)*phi_c(i)*cosd(theta(i))*(1-z_pe{i,1});
    
    vec = sim_list{i,1};
    % Processes the data differently depending on the simulation type
    if (sim_type(i) == "dil")
        p_e{i,1} = vec(:,1:N(i))';
        p_p{i,1} = p_b{i,1}-p_e{i,1};
        phi{i,1} = vec(:,N(i)+1:2*N(i))';
        u_f{i,1} = vec(:,2*N(i)+1:3*N(i))';
        u_p{i,1} = vec(:,3*N(i)+1:end)';
    elseif (sim_type(i) == "pcon")
        p_e{i,1} = s_frac(i)*p_b{i,1};
        p_p{i,1} = (1-s_frac(i))*p_b{i,1};
        % If looking at initial values for dilatancy sim, need to set phi
        % to zero here
        phi{i,1} = alpha_dl(i)*p_p{i,1};
        u_f{i,1} = vec(:,1:N(i))';
        u_p{i,1} = vec(:,N(i)+1:end)';
    elseif (sim_type(i) == "pdriv")
        p_e{i,1} = vec(:,1:N(i))';
        p_p{i,1} = (p_b{i,1}-p_e{i,1});
        phi{i,1} = alpha_dl(i)*p_p{i,1};
        u_f{i,1} = vec(:,N(i)+1:2*N(i))';
        u_p{i,1} = vec(:,2*N(i)+1:end)';
    elseif (sim_type(i) == "ucon")
        p_e{i,1} = vec(:,1:N(i))';
        p_p{i,1} = (p_b{i,1}-p_e{i,1});
        phi{i,1} = vec(:,N(i)+1:end)';
        initial_u = load("no_phi_no_p.txt");
        u_f{i,1} = initial_u(7,1:N(i))'.*(ones(1,11));
        u_p{i,1} = initial_u(7,N(i)+1:end)'.*(ones(1,11));
    end
    beta_pe{i,1} = 150*(phi_c(i) + phi{i,1}).^2.*eta_f_dl(i)./((1-phi_c(i)-phi{i,1}).^3.*d_dl(i)^2);
    beta_u{i,1} = interp1(z_pe{i,1},beta_pe{i,1},z_u{i,1},'linear','extrap');
    % If the velocity evolves we need to define these quantities
    dpdz{i,1} = vertcat(zeros(1,n_times(i)),diff(p_e{i,1},1,1))./dz_dl(i);
    d2pdz2{i,1} = vertcat(diff(dpdz{i,1},1,1),zeros(1,n_times(i)))./dz_dl(i);
    d3pdz3{i,1} = vertcat(zeros(1,n_times(i)),diff(d2pdz2{i,1},1,1))./dz_dl(i);
    dppdz{i,1} = vertcat(zeros(1,n_times(i)),diff(p_p{i,1},1,1))./dz_dl(i);
    dufdz{i,1} = vertcat(diff(u_f{i,1},1,1),zeros(1,n_times(i)))./dz_dl(i);
    d2ufdz2{i,1} = vertcat(zeros(1,n_times(i)),diff(dufdz{i,1},1,1))./dz_dl(i);
    dupdz{i,1} = vertcat(diff(u_p{i,1},1,1),zeros(1,n_times(i)))./dz_dl(i);
    d2updz2{i,1} = vertcat(zeros(1,n_times(i)),diff(dupdz{i,1},1,1))./dz_dl(i);
    Iv_temp = eta_f_dl(i).*abs(dupdz{i,1})./(p_p{i,1}+1e-8);
    Iv_temp(N(i),:) = eta_f_dl(i).*d2updz2{i,1}(N(i),:)./(-buoyancy(i)-dpdz{i,1}(N(i),:));
    Iv{i,1} = Iv_temp;
    if (sim_type(i)~="ucon")
        tau_f{i,1} = eta_f_dl(i).*d2ufdz2{i,1}./(1-phi_c(i));
    
        I{i,1} = 2.*d_dl(i).*abs(dupdz{i,1}).*sqrt(density_ratio(i))./sqrt(p_p{i,1}.*ones(1,n_times(i))+1e-8);
        mu{i,1} = mu_I_fn(I{i,1});
        tau_p{i,1} = 1/(phi_c(i).*density_ratio(i)).*sign(dupdz{i,1}).*vertcat(zeros(1,n_times(i)),diff(mu{i,1}.*(p_p{i,1}.*ones(1,n_times(i))),1,1)./dz_dl(i));
    
        drag_mult_p{i,1} = (1-phi_c(i))^2.*beta_u{i,1}./(density_ratio(i)*phi_c(i)).*ones(1,n_times(i));
    
        drag_term_p{i,1} = drag_mult_p{i,1}.*(u_f{i,1}-u_p{i,1});
        
        drag_mult_f{i,1} = (1-phi_c(i)).*beta_u{i,1}.*ones(1,n_times(i));
    
        drag_term_f{i,1} = drag_mult_f{i,1}.*(u_f{i,1}-u_p{i,1});
        
        dupdt{i,1} = vertcat(zeros(1,n_times(i)),tau_p{i,1}(2:end,:)+sind(theta(i))+drag_term_p{i,1}(2:end,:));
        dufdt{i,1} = vertcat(zeros(1,n_times(i)),tau_f{i,1}(2:end,:)+sind(theta(i))-drag_term_f{i,1}(2:end,:));
    end
    
    
    % If the pressure evolves, have to define these quantities
    if ((sim_type(i)=="pdriv") || (sim_type(i)=="dil") || (sim_type(i)=="ucon"))
        dphidz{i,1} = vertcat(zeros(1,n_times(i)),diff(phi{i,1},1,1))./dz_dl(i);
        dpdt{i,1} = 1./(alpha_dl(i)).*vertcat(diff(1./beta_pe{i,1}.*dpdz{i,1}),zeros(1,n_times(i)))./dz_dl(i);
        if ((sim_type(i)=="dil") || (sim_type(i)=="ucon"))
            dilatancy{i,1} = -1./(alpha_dl(i)).*dupdz{i,1}.*(phi{i,1}+sqrt(abs(Iv{i,1}))*phi_c(i)./(1+sqrt(abs(Iv{i,1}))));
            diffusion_term{i,1} = dpdt{i,1};
            dphidt{i,1} = -dpdt{i,1}.*alpha_dl(i).*phi_c(i);
            dpdt{i,1} = dpdt{i,1} + dilatancy{i,1};
            d2pdzdt = vertcat(zeros(1,n_times(i)),diff(dpdt{i,1},1,1))./dz_dl(i);
        end
    end
end
%% 

for i=1:sim_num-2

    % Set the quantity to be plotted here

%     test = 1/(phi_c(i).*density_ratio(i)).*sign(dupdz{i,1}).*mu{i,1}.*vertcat(zeros(1,n_times(i)),diff((p_p{i,1}.*ones(1,n_times(i))),1,1)./dz_dl(i));
%     plot_vec = vertcat(zeros(1,n_times(i)),test(2:end,:)+sind(theta(i))+drag_term_p{i,1}(2:end,:));

%     [M,Index] = min(dilatancy{i,1});
%     crit_pp = get_critical_value(u_p{i,1}',p_e{i,1}',true,1e-5);
%     crit_dppdz = get_critical_value(u_p{i,1}',dppdz{i,1}',false,1e-5);
%     crit_tau_f = get_critical_value(u_p{i,1}',tau_f{i,1}'*(1-phi_c(i)),false,1e-5);
%     crit_p_deriv = get_critical_value(u_p{i,1}',(vertcat(zeros(1,n_times(i)),diff(p_p{i,1},1,1))./dz_dl(i).*mu{i,1})',false,1e-5);
%     crit_mu_deriv = get_critical_value(u_p{i,1}',(vertcat(zeros(1,n_times(i)),diff(mu{i,1},1,1))./dz_dl(i).*p_p{i,1})',false,1e-5);
%     

%     s_c = 1-(1-phi_c(i)+density_ratio(i)*phi_c(i))/((density_ratio(i)-1)*phi_c(i))*tand(theta(i))/0.342;
%      plot(crit_tau_f,z_u{i,1},'DisplayName','Actual Value')
%     hold on
%     plot(crit_p_deriv,z_u{i,1},'DisplayName','Actual Value')
%     plot(crit_mu_deriv,z_u{i,1},'DisplayName','Actual Value')
%     plot(crit_dppdz,z_u{i,1},'DisplayName','Actual Value')
    a=0.0107;
    b=-0.0101;
    p=0.9169;%0.9078;
    c=0.01128;
%     plot(-crit_grad*(1-z_u{i,1}),z_u{i,1})
    approx_up = (-0.25*p*a^2*(1-z_u{i,1}).^4-2/3*a*b*p*(1-z_u{i,1}).^3-b^2*p/2*(1-z_u{i,1}).^2)/(eta_f_dl(i)*phi_c(i)^2)+c;
%     comp_vec = p_p{i,1}.*(phi{i,1}).^2/(phi_c(i)^2*eta_f_dl(i));
    approx_dupdz = (a*(1-z_u{i,1})+b).^2.*p.*(1-z_u{i,1})./(phi_c(i)^2*eta_f_dl(i));
%     plot(approx_up,z_pe{i,1}(1:end))
%     test_vec = a*(1-z_u{i,1})+b;
    crit_grad = -(density_ratio(i)*phi_c(i)+(1-phi_c(i)))*sind(theta(i))/0.342;
    M = 3.1159*2.*d_dl(i).*sqrt(density_ratio(i));
    approx_f = (eta_f_dl(i)+M*sqrt(-crit_grad*(1-z_u{i,1}))).*approx_dupdz/(0.342);

%     erosion_point = get_erosion_time(u_p{i,1}',5e-6);
%     pe_grad = zeros(6001,1);
%     for k=1:6001
%         pe_grad(k,1) = dpdz{i,1}(round(max(erosion_point(k)*N(i)-3,1)),k);
%     end
%     ep_diff = erosion_point(430:2030)-erosion_point(400:2000);
%     ep_grad = zeros(1601,1);
%     ep_grad(1)=(ep_diff(1)+ep_diff(2)+ep_diff(3))/3;
%     ep_grad(2)=(ep_diff(1)+ep_diff(2)+ep_diff(3)+ep_diff(4))/4;
%     for l=3:1599
%         ep_grad(l)=(ep_diff(l-2)+ep_diff(l-1)+ep_diff(l)+ep_diff(l+1)+ep_diff(l+2))/5;
%     end
%     ep_grad(1600)=(ep_diff(1598)+ep_diff(1599)+ep_diff(1600)+ep_diff(1601))/4;
%     ep_grad(1601)=(ep_diff(1599)+ep_diff(1600)+ep_diff(1601))/3;
%     plot_vec = -0.342*((crit_grad+8e-5).*(1-z_pe{i,1}.*ones(1,6001))+p_p{i,1})./(eta_f_dl(i)+M*sqrt(p_p{i,1}));%vertcat(zeros(1,n_times(i)),diff(mu{i,1},1,1))./dz_dl(i).*p_p{i,1};
    
    approx_pe = p_b{i,1}+(crit_grad+4e-5*(2.5*0.6+0.4)).*(1-z_pe{i,1})+approx_f;
    % Plots initial profile  
%     subplot(sim_num,1,i);
    % Plots profile at other times
%     SetPaperSize(10,10);
    
%     plot(linspace(40,200,1601),dupdt{i,1}(155,400:2000))

    p_p_deriv = (vertcat(zeros(1,n_times(i)),diff(p_p{i,1},1,1))./dz_dl(i).*mu{i,1});
    mu_deriv = (vertcat(zeros(1,n_times(i)),diff(mu{i,1},1,1))./dz_dl(i).*p_p{i,1});
    
    plot_vec = u_p{i,1};
    t1=[110,5,2000];
    nfig=5;
    for j=linspace(t1(i),t1(i)+(nfig-1)*10,nfig)
        t_val = j;
        plot(plot_vec(1:end,j),z_pe{i,1},'DisplayName',"t="+num2str(t_val));
        hold on
%         plot(approx_pe,z_pe{i,1},'DisplayName',"Approximated Profile");
%         plot(plot_vec(1:end,j)+(crit_grad-8e-5)*(1-z_u{i,1}),z_u{i,1},'DisplayName',"t="+num2str(t_val));       
        
%         plot(-13.26*(1-z_u{i,1}).^(3/2)+13.26,z_u{i,1})
%         plot(diffusion_term{i,1}(1:end,j),z_pe{i,1}(1:end),'DisplayName',"t="+num2str(t_val));
%         plot(tau_f{i,1}(:,j),z_pe{i,1}(1:end),'DisplayName','$\frac{\partial \tau_f}{\partial z}$')
%         plot(p_p_deriv(:,j),z_pe{i,1}(1:end),'DisplayName','$\frac{\partial p_p}{\partial z} \mu$')
%         plot(mu_deriv(:,j),z_pe{i,1}(1:end),'DisplayName','$\frac{\partial \mu}{\partial z} p_p$')
    end
%     plot((-buoyancy(i)+crit_grad)*(1-z_u{i,1}),z_u{i,1},"r--")
    legend('Location', "northwest",'UserData', 13);
    ylabel("z/h");
    xlabel('$p_e$');

    
    title('Excess Pressure and the Approximated Profile');

%      ax = gca;
% %      ax.YAxis(2).Exponent = 0;
%      ax.XAxis.Exponent = 0;
% %      ytickformat('%.0e');
%      xtickformat('%.1e');
% %     hold off
%     pdfname = 'dil9_5deg_pe_appox';
%     PrintFig(pdfname)
%     movefile([pdfname '.pdf'], 'Figures/Dilatant/')
end

ani_num = 1;
make_avi=false;
if (make_avi) 
    nframes=1500;
    t_vals = (0:nframes)*0.1;
    Frames=moviein(nframes);
    figure
    v = VideoWriter('dil9_5_low_phi_erosion.avi');
    v.FrameRate=50;
    open(v)
    ax = gca();
    for i=500:nframes+500
        subplot(2,2,1)
        plot(p_e{ani_num,1}(:,i), z_pe{ani_num,1});
        hold on
        plot((-buoyancy(ani_num)-(density_ratio(ani_num)*phi_c(ani_num)+(1-phi_c(ani_num)))*sind(theta(ani_num))/0.342)*(1-z_pe{ani_num,1}), z_u{2,1},"--");
        hold off
        legend("t="+num2str(t_vals(i-499)+50,'%.1f'),"Critical $p_e$","Location",'northwest','Interpreter','latex');
%         legend();
        xlim([-0.05 0]);
        ylim([0 1]);
        xlabel("$p_e$",'Interpreter','latex');
        ylabel("$z$",'Interpreter','latex');

        subplot(2,2,2)
        plot(phi{ani_num,1}(:,i), z_u{ani_num,1});
%         hold on
%         plot(phi{2,1}(:,i+29), z_u{2,1});
%         hold off
        xlim([-0.03,0])
        ylim([0 1])
        xlabel('$\hat{\phi}$','Interpreter','latex');
        ylabel("$z$",'Interpreter','latex');

        subplot(2,2,3)
        plot(u_f{ani_num,1}(:,i), z_u{ani_num,1});
%         hold on
%         plot(u_f{2,1}(:,i+29), z_u{2,1});
%         hold off
        xlim([0,0.05])
        ylim([0 1])
        xlabel("$u_f$",'Interpreter','latex');
        ylabel("$z$",'Interpreter','latex');

        subplot(2,2,4)
        plot(u_p{ani_num,1}(:,i), z_u{ani_num,1});
%         hold on
%         plot(u_p{2,1}(:,i+29), z_u{2,1});
%         hold off
        xlim([0,0.05])
        ylim([0 1])
        xlabel("$u_p$",'Interpreter','latex');
        ylabel("$z$",'Interpreter','latex','Interpreter','latex');
        drawnow();
        writeVideo(v,getframe(gcf))
    end
    close all
    close(v)
end

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