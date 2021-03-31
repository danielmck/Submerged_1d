opengl software
cd Results
sim1 = load("dil9_5deg_low_phi_small_early.txt");
% sim2 = load("dil9_5deg_dl_small.txt");
% sim3 = load("sim3_out.txt");
cd ../
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
sim_type = ["dil","dil","pdriv"];

N=[200,400,200];
h = [4e-3,4e-3,4e-3];
d=[1.43e-5, 1.43e-5,1.43e-4];
dz = h./(N-0.5); % z spacing
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

d_dl = d./z_scale;
dz_dl = dz./z_scale;
eta_f_dl = eta_f./(p_scale.*t_scale);
alpha_dl = alpha.*p_scale;

% Creates cells to store the matrices of values

buoyancy = zeros(sim_num);
buoyancy_dl = zeros(sim_num);
z_pe = cell(sim_num,1);
z_u = cell(sim_num,1);
p_b = cell(sim_num,1);
p_e = cell(sim_num,1);
p_p = cell(sim_num,1);
phi = cell(sim_num,1);
u_f = cell(sim_num,1);
u_p = cell(sim_num,1);
z_pe_dl = cell(sim_num,1);
z_u_dl = cell(sim_num,1);
p_b_dl = cell(sim_num,1);
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
    z_pe_dl{i,1} = z_pe{i,1}./z_scale(i);
    z_u_dl{i,1} = z_u{i,1}./z_scale(i);
    p_b_dl{i,1} = p_b{i,1}/p_scale(i);
    buoyancy_dl(i) = buoyancy(i)/(p_scale(i));
    
    vec = sim_list{1,i};
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
        phi{i,1} = alpha(i)*p_p{i,1};
        u_f{i,1} = vec(:,1:N(i))';
        u_p{i,1} = vec(:,N(i)+1:end)';
    elseif (sim_type(i) == "pdriv")
        p_e{i,1} = vec(:,1:N(i))';
        p_p{i,1} = (p_b{i,1}-p_e{i,1});
        phi{i,1} = alpha(i)*p_p{i,1};
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
    beta_u{i,1} = interp1(z_pe_dl{i,1},beta_pe{i,1},z_u_dl{i,1},'linear','extrap');
    % If the velocity evolves we need to define these quantities
    dpdz{i,1} = vertcat(zeros(1,n_times(i)),diff(p_e{i,1},1,1))./dz_dl(i);
    d2pdz2{i,1} = vertcat(diff(dpdz{i,1},1,1),zeros(1,n_times(i)))./dz_dl(i);
    d3pdz3{i,1} = vertcat(zeros(1,n_times(i)),diff(d2pdz2{i,1},1,1))./dz_dl(i);
    dufdz{i,1} = vertcat(diff(u_f{i,1},1,1),zeros(1,n_times(i)))./dz_dl(i);
    d2ufdz2{i,1} = vertcat(zeros(1,n_times(i)),diff(dufdz{i,1},1,1))./dz_dl(i);
    dupdz{i,1} = vertcat(diff(u_p{i,1},1,1),zeros(1,n_times(i)))./dz_dl(i);
    d2updz2{i,1} = vertcat(zeros(1,n_times(i)),diff(dupdz{i,1},1,1))./dz_dl(i);
    Iv_temp = eta_f_dl(i).*abs(dupdz{i,1})./(p_p{i,1}+1e-8);
    Iv_temp(N(i),:) = eta_f_dl(i).*d2updz2{i,1}(N(i),:)./(-buoyancy_dl(i)-dpdz{i,1}(N(i),:));
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

for i=1:sim_num
    % Set the quantity to be plotted here

%     test = 1/(phi_c(i).*density_ratio(i)).*sign(dupdz{i,1}).*mu{i,1}.*vertcat(zeros(1,n_times(i)),diff((p_p{i,1}.*ones(1,n_times(i))),1,1)./dz_dl(i));
%     plot_vec = vertcat(zeros(1,n_times(i)),test(2:end,:)+sind(theta(i))+drag_term_p{i,1}(2:end,:));

%     [M,Index] = min(dilatancy{i,1});

    s_c = 1 - (1-phi_c(i)+density_ratio(i)*phi_c(i))/((density_ratio(i)-1)*phi_c(i))*tand(theta(i))/0.342;
    plot_vec = u_p{i,1};%vertcat(zeros(1,n_times(i)),diff(mu{i,1},1,1))./dz_dl(i);
    % Plots initial profile
%     plot(plot_vec(:diffusion_term,1),z_u{i,1}, 'DisplayName', 't=0.0'); 
%     SetPaperSize(15,10)
    subplot(sim_num,1,i);
    hold on
    % Plots profile at other times
    t1=[1500];
    nfig=1;
    a=0.0107;
    b=-0.0101;
    p=0.9078;
    c=0.01128;
    approx_up = 0.25*p*(-a^2*(1-z_u{i,1}).^4-2/3*a*b*(1-z_u{i,1}).^3-b^2*(1-z_u{i,1}).^2)/(eta_f_dl(i)*phi_c(i)^2)+c;
    comp_vec = p_p{i,1}.*(phi{i,1}).^2/(phi_c(i)^2*eta_f_dl(i));
    approx_dupdz = (a*(1-z_u{i,1})+b).^2.*p.*(1-z_u{i,1})./(phi_c(i)^2*eta_f_dl(i));
    test_vec = a*(1-z_u{i,1})+b;
    for j=linspace(t1(i),t1(i)+(nfig-1)*300,nfig)
        t_val = j;
%         plot(dilatancy{i,1}(1:end,j),z_u{i,1}(1:end),'DisplayName','Dilatant Term');
%         plot(diffusion_term{i,1}(1:end,j),z_u{i,1}(1:end),'DisplayName','Diffusion Term');
%         plot(linspace(-5e-5,5e-5,100),ones(1,100)*shear_end(86)/N(i))
%         plot(dpdt{i,1}(1:end,j),z_u{i,1}(1:end),'DisplayName','$\frac{\partial p_e}{\partial t}$');
        plot(plot_vec(110:end,j),z_u{i,1}(110:end),'DisplayName',"t="+num2str(t_val));
%         plot(comp_vec(110:end,j),z_u{i,1}(110:end),'DisplayName',"t="+num2str(t_val));
%         plot(approx_dupdz(110:end),z_u{i,1}(110:end),'DisplayName',"t="+num2str(t_val));
        plot(approx_up(110:end),z_u{i,1}(110:end))
%         plot((dupdt{i,1}(:,j)>2e-6)-(p_e{i,1}(:,j)>s_c*(1-z_u{i,1})) ,z_u{i,1})
    end
%     plot(-sind(theta(i))*(phi_c(i)*density_ratio(i)+1-phi_c(i))*ones(N(i),1),z_u{i,1},'--') % 
%     yyaxis left
%     plot(linspace(15,1000,986).*50,Index(15:1000)./N(i),'DisplayName', 'Peak Position')
    ylabel("z/h");
    xlabel('$\frac{\partial p_e}{\partial t}$');
    legend();
%     xlim([-0.22 -0.2])
%     hold on
%     yyaxis right
%     plot(linspace(15,1000,986).*50,M(15:1000),'DisplayName', 'Peak Magnitude')
%     hold off
%     ylabel('Dilatancy')
%     xlabel('t');
    
%     xlim([0,0.5])
    title('Contributions of the Dilatant and Diffusive Terms at t=1050');
%     legend('UserData',12,'Location','NorthEast');
%      ax = gca;
% %      ax.YAxis(2).Exponent = 0;
%      ax.XAxis.Exponent = 0;
% %      ytickformat('%.0e');
%      xtickformat('%.1e');
% %     hold off
%     pdfname = 'dil9_5_dil_diff_dpdt_low_phi_small';
%     PrintFig(pdfname)
%     movefile([pdfname '.pdf'], 'Figures/Dilatant/')
end
% legend("0s","0.45s","0.9s","1.35s","1.8s","2.25s","2.7s","3.15s","3.6s","4.05s","4.5s", 'Location',"east");
% saveas(gcf,"dpdt_evolution.png")



make_animation = false;
ani_num = 1;

% if (make_animation)
%     axis tight manual
%     nr_fr = 5000;
%     t_vals = (1:nr_fr)*0.5;
%     figure(1)
%     animation_name = 'dil9_5_erosion_small.gif';
%     for i = 1 : nr_fr
% %         plot(dilatancy{ani_num,1}(:,i), z_u{ani_num,1});
% % %         xlim([0 2e-4]);
% %         ylim([0 1]);
% %         xlabel("Excess Pressure Diffusion");
% %         ylabel("z");
% %         legend("t="+num2str(t_vals(i)));
% %         xlim([-0.1 0]);
% %         ylim([0 1]);
% %         xlabel("p_e");
% %         ylabel("z");
% 
%         subplot(2,2,1)
%         plot(p_e{ani_num,1}(:,i), z_u{ani_num,1}, 'DisplayName','High u_p IC');
% %         hold on
% %         plot(p_e{2,1}(:,i+29), z_u{2,1}, 'DisplayName', 'High p_p IC');
% %         hold off
%         legend("t="+num2str(t_vals(i),'%.1f'));
% %         legend();
%         xlim([-0.05 0]);
%         ylim([0 1]);
%         xlabel("p_e");
%         ylabel("z");
% 
%         subplot(2,2,2)
%         plot(phi{ani_num,1}(:,i), z_u{ani_num,1});
% %         hold on
% %         plot(phi{2,1}(:,i+29), z_u{2,1});
% %         hold off
%         xlim([-0.03,0])
%         ylim([0 1])
%         xlabel('\phi-\phi_c');
%         ylabel("z");
% 
%         subplot(2,2,3)
%         plot(u_f{ani_num,1}(:,i), z_u{ani_num,1});
% %         hold on
% %         plot(u_f{2,1}(:,i+29), z_u{2,1});
% %         hold off
%         xlim([0,0.1])
%         ylim([0 1])
%         xlabel("u_f");
%         ylabel("z");
% 
%         subplot(2,2,4)
%         plot(u_p{ani_num,1}(:,i), z_u{ani_num,1});
% %         hold on
% %         plot(u_p{2,1}(:,i+29), z_u{2,1});
% %         hold off
%         xlim([0,0.1])
%         ylim([0 1])
%         xlabel("u_p");
%         ylabel("z");
%         drawnow
%         frame = getframe(1);
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
%         if i == 1
%           imwrite(imind,cm,animation_name,'gif', 'LoopCount',Inf);
%         else
%           imwrite(imind,cm,animation_name,'gif','WriteMode','append','DelayTime', 0.05);
%         end
%     end
% end
make_avi=false;
if (make_avi) 
    nframes=6000;
    t_vals = (1:nframes);
    Frames=moviein(nframes);
    figure
    v = VideoWriter('dil9_5_low_phi_small.avi');
    v.FrameRate=50;
    open(v)
    ax = gca();
    for i=1:nframes
        subplot(2,2,1)
        plot(p_e{ani_num,1}(:,i), z_u{ani_num,1});
%         hold on
%         plot(p_e{2,1}(:,i+29), z_u{2,1}, 'DisplayName', 'High p_p IC');
%         hold off
        legend("t="+num2str(t_vals(i),'%.1f'));
%         legend();
        xlim([-0.05 0]);
        ylim([0 1]);
        xlabel("p_e");
        ylabel("z");

        subplot(2,2,2)
        plot(phi{ani_num,1}(:,i), z_u{ani_num,1});
%         hold on
%         plot(phi{2,1}(:,i+29), z_u{2,1});
%         hold off
        xlim([-0.03,0])
        ylim([0 1])
        xlabel('\phi-\phi_c');
        ylabel("z");

        subplot(2,2,3)
        plot(u_f{ani_num,1}(:,i), z_u{ani_num,1});
%         hold on
%         plot(u_f{2,1}(:,i+29), z_u{2,1});
%         hold off
        xlim([0,0.1])
        ylim([0 1])
        xlabel("u_f");
        ylabel("z");

        subplot(2,2,4)
        plot(u_p{ani_num,1}(:,i), z_u{ani_num,1});
%         hold on
%         plot(u_p{2,1}(:,i+29), z_u{2,1});
%         hold off
        xlim([0,0.1])
        ylim([0 1])
        xlabel("u_p");
        ylabel("z");

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