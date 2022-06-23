names = ["Rauter_5_deep.txt"];
sim_num=size(names,2);
sim_list = cell(1,sim_num);
sim_title = ["$5$ degree Slope"];

N=zeros(1,sim_num);
h = zeros(1,sim_num);
d=zeros(1,sim_num);
dz = zeros(1,sim_num);
phi_c=zeros(1,sim_num);
eta_f_dl = zeros(1,sim_num);
g=9.81*ones(1,sim_num); % m/s^2
theta = zeros(1,sim_num); % deg
a_dl = zeros(1,sim_num); % 1/Pa
density_ratio = zeros(1,sim_num);
shear_lim_dl = zeros(1,sim_num);
phi_rcp = zeros(1,sim_num);
creep_type = zeros(1,sim_num);
phi_rlp = zeros(1,sim_num);
t_step = zeros(1,sim_num);
n_times = zeros(sim_num,1);
crit_grad = zeros(sim_num,1);

record = readtable('EqnOfState_Results/result_record.csv');

for k=1:sim_num
    in_table = strcmp(record.Name, names(k));
    cd EqnOfState_Results
    sim_list{k,1} = load(names(k));
    cd ../
    N(k) = record.N(in_table);
    h(k) = record.h(in_table);
    d(k) = record.d(in_table);
    dz(k) = h/(N-0.5);
    phi_c(k) = record.phi_c(in_table);
    phi_rcp(k) = record.phi_rcp(in_table);
    phi_rlp(k) = record.phi_rlp(in_table);
    creep_type(k) = record.creep_type(in_table);
    shear_lim_dl(k) = record.S0(in_table);
    density_ratio(k) = record.rho_r(in_table);
    eta_f_dl(k) = record.eta_f(in_table);
    theta(k) = record.theta(in_table);
    a_dl(k) = record.a(in_table);
    t_step(k) = record.t_step(in_table);
    n_times(k) = size(sim_list{k,1},1);
end

v_scale = sqrt(g.*h);
t_scale = sqrt(h./g);
z_scale = h;
dz = h./(N-0.5); % z spacing

d_dl = d./z_scale;
dz_dl = dz./z_scale;

z_pe = cell(sim_num,1);
z_u = cell(sim_num,1);
z_pe_dl = cell(sim_num,1);
z_u_dl = cell(sim_num,1);
p_b_dl = cell(sim_num,1);
p_e = cell(sim_num,1);
p_p = cell(sim_num,1);
p_c = cell(sim_num,1);
p_i = cell(sim_num,1);
phi_hat = cell(sim_num,1);
phi = cell(sim_num,1);
u_f = cell(sim_num,1);
u_p = cell(sim_num,1);
Iv = cell(sim_num,1);
tan_psi = cell(sim_num,1);
Iv_phi = cell(sim_num,1);
mu_Iv = cell(sim_num,1);
beta_pe = cell(sim_num,1);
beta_u = cell(sim_num,1);
tau_p = cell(sim_num,1);
tau_f = cell(sim_num,1);
drag_force = cell(sim_num,1);
dphidz = cell(sim_num,1);
phi_Iv = cell(sim_num,1);
dpdz = cell(sim_num,1);
d2pdz2 = cell(sim_num,1);
dufdz = cell(sim_num,1);
d2ufdz2 = cell(sim_num,1);
d2updz2 = cell(sim_num,1);
dupdz = cell(sim_num,1);
dpidt = cell(sim_num,1);
dpcdt = cell(sim_num,1);
dphidt = cell(sim_num,1);
dupdt = cell(sim_num,1);
dufdt = cell(sim_num,1);
dudt = cell(sim_num,1);
smooth_dpidt = cell(sim_num,1);
smooth_dpcdt = cell(sim_num,1);
smooth_dphidt = cell(sim_num,1);
smooth_dupdt = cell(sim_num,1);
smooth_dpedt = cell(sim_num,1);
smooth_dgddt = cell(sim_num,1);
dpidt_gd = cell(sim_num,1);
dpidt_phi = cell(sim_num,1);
phi_m = cell(sim_num,1);
excess_p_c = cell(sim_num,1);
p_crit = cell(sim_num,1);
phi_crit = cell(sim_num,1);
dpidphi = cell(sim_num,1);
dpidgd = cell(sim_num,1);
pe_ave = cell(sim_num,1);
pi_ave = cell(sim_num,1);
pc_ave = cell(sim_num,1);
up_ave = cell(sim_num,1);
phi_ave = cell(sim_num,1);
denom_ave = cell(sim_num,1);

for i = 1:sim_num
    z_pe{i,1} = linspace(dz(i)/2,h(i),N(i));
    z_u{i,1} = linspace(0,h(i)-dz(i)/2,N(i));
        
    z_pe_dl{i,1} = z_pe{i,1}./z_scale(i);
    z_u_dl{i,1} = z_u{i,1}./z_scale(i);
    p_b_dl{i,1} = (density_ratio(i)-1).*phi_c(i).*cosd(theta(i)).*(1-z_pe_dl{i,1});
    phi_hat{i,1} = sim_list{i,1}(:,1:200);
    u_f{i,1} = sim_list{i,1}(:,201:400);
    u_p{i,1} = sim_list{i,1}(:,401:600);
    phi{i,1}=phi_c(i)+phi_hat{i,1};

    dufdz{i,1}=diff(u_f{i,1},1,2)./dz_dl(i);
    dupdz{i,1}=diff(u_p{i,1},1,2)./dz_dl(i);

    p_c{i,1} = a_dl(i)*(phi{i,1}-phi_rlp(i))./(phi_rcp(i)-phi{i,1});
    p_c{i,1}(p_c{i,1}<0) = 0;
    
    if (creep_type(i) == 0)
        phi_m{i,1} = phi_c(i);
    else
        phi_m{i,1} = phi_c(i)+(phi_rcp(i)-phi_c(i)).*(abs(dupdz{i,1})<shear_lim_dl(i)).*(shear_lim_dl(i)-abs(dupdz{i,1})).^creep_type(i)/shear_lim_dl(i)^creep_type(i);
    end
    p_i{i,1} = eta_f_dl(i).*abs(dupdz{i,1})./((phi_m{i,1}./(phi{i,1}(:,1:end-1))-1).^2);

    p_p{i,1} = horzcat(p_c{i,1}(:,1:end-1) + p_i{i,1},zeros(n_times(i),1));
    p_e{i,1} = p_b_dl{i,1}.*ones(n_times(i),1)-p_p{i,1}; 

    Iv{i,1} = eta_f_dl(i).*abs(dupdz{i,1})./(p_p{i,1}(:,1:end-1));
    phi_Iv{i,1} = phi_c(i)./(1+sqrt(abs(Iv{i,1})));
    tan_psi{i,1} = phi{i,1}-horzcat(phi_Iv{i,1}, phi_rcp(i)*ones(n_times(i),1));
    crit_grad(i) = -(density_ratio(i).*phi_c(i)+(1-phi_c(i))).*sind(theta(i))/0.32;
    excess_p_c{i,1} = p_b_dl{i,1}+crit_grad(i)*(1-z_pe_dl{i,1});
    p_crit{i,1} = -crit_grad(i)*(1-z_pe_dl{i,1});
    phi_crit{i,1} = (phi_rcp(i)*p_crit{i,1}+a_dl(i)*phi_rlp(i))./(a_dl(i)+p_crit{i,1});

    dufdz{i,1}=diff(u_f{i,1},1,2)./dz_dl(i);
    dupdz{i,1}=diff(u_p{i,1},1,2)./dz_dl(i);

    dpdz{i,1} = horzcat(zeros(n_times(i),1), diff(p_e{i,1},1,2)./dz_dl(i));
    d2pdz2{i,1} = horzcat(diff(dpdz{i,1},1,2)./dz_dl(i),zeros(n_times(i),1));
    dphidz{i,1} = diff(phi{i,1},1,2)./dz_dl(i);

    Iv_phi{i,1} = (phi_c(i)./phi{i,1}-1).^2; 

    beta_pe{i,1} = beta_fn(phi_hat{i,1})*150.*eta_f_dl(i)./(d_dl(i)^2);
    beta_u{i,1} = interp1(z_pe{i,1},beta_pe{i,1}',z_u{i,1},'linear','extrap')';

    dphidt{i,1} = -phi_c(i).*diff((1./beta_u{i,1}.*dpdz{i,1}),1,2)./dz_dl(i);
    dphidt{i,1} = interp1(z_pe{i,1}(1:end-1),dphidt{i,1}',z_pe{i,1},'linear','extrap')';

    mu_Iv{i,1} = mu_Iv_fn(Iv{i,1});
    tau_f{i,1} = horzcat(eta_f_dl(i).* dufdz{i,1}, zeros(n_times(i),1));
    tau_p{i,1} = horzcat(mu_Iv{i,1}.*sign(dupdz{i,1}).*p_p{i,1}(:,1:end-1), zeros(n_times(i),1));
    %         diff_tau_p = (mu_I_fn(I).*[0 diff(sign(dupdz).*p_p)]+(abs(I_u)>10/reg_param).*sign(dupdz).*p_p.*[0 diff(mu_I_fn(I))])./dz_dl;
    drag_force{i,1} = (1-phi_c(i))^2.*beta_u{i,1}.*(u_f{i,1}-u_p{i,1});
    dufdt{i,1} = horzcat(zeros(n_times(i),1), (1/(1-phi_c(i)).*diff(tau_f{i,1},1,2)./dz_dl(i) + sind(theta(i))-drag_force{i,1}(:,2:end)./(1-phi_c(i))));
    dupdt{i,1} = horzcat(zeros(n_times(i),1), (1/(density_ratio(i)*phi_c(i)).*diff(tau_p{i,1},1,2)./dz_dl(i) + sind(theta(i))+drag_force{i,1}(:,2:end)./(density_ratio(i)*phi_c(i))));

    dpidt{i,1} = a_dl(i)*(phi_rcp(i)-phi_rlp(i))./(phi_rcp(i)-phi{i,1}).^2.*dphidt{i,1};
    dpcdt{i,1} = horzcat(2*phi_c(i)*eta_f_dl(i).*abs(dupdz{i,1})./((phi_m{i,1}./(phi{i,1}(:,1:end-1))-1).^3).*dphidt{i,1}(:,1:end-1)+eta_f_dl(i)./((phi_m{i,1}./(phi{i,1}(:,1:end-1))-1).^2).*diff(dupdt{i,1},1,2), zeros(n_times(i),1));

    smooth_dpidt{i,1} = smooth_diff(p_i{i,1},5,t_step(i));
    smooth_dpcdt{i,1} = smooth_diff(p_c{i,1},5,t_step(i));
    smooth_dpedt{i,1} = smooth_diff(p_e{i,1},5,t_step(i));
    smooth_dgddt{i,1} = smooth_diff(dupdz{i,1},5,t_step(i));
    smooth_dphidt{i,1} = smooth_diff(phi_hat{i,1},5,t_step(i));
    smooth_dupdt{i,1} = smooth_diff(u_p{i,1},5,t_step(i));
    
    pe_ave{i,1} = depth_average(p_e{i,1}',N(i),n_times(i));

    dpidphi{i,1} = eta_f_dl(i).*abs(dupdz{i,1}).*phi_m{i,1}.*(phi{i,1}(:,1:end-1))./((phi_m{i,1}-(phi{i,1}(:,1:end-1))).^3);
    dpidgd{i,1} = eta_f_dl(i)./((phi_m{i,1}./(phi{i,1}(:,1:end-1))-1).^2)+4.*eta_f_dl(i).*abs(dupdz{i,1})./((phi_m{i,1}./(phi{i,1}(:,1:end-1))-1).^3).*(abs(dupdz{i,1})<shear_lim_dl(i))./(phi{i,1}(:,1:end-1)).*(shear_lim_dl(i)-abs(dupdz{i,1}))/shear_lim_dl(i)^2;

    dpidt_phi{i,1} = dpidphi{i,1}.*smooth_dphidt{i,1}(:,1:end-1);
    dpidt_gd{i,1} = dpidgd{i,1}.*smooth_dgddt{i,1}(:,1:end);
    
    pi_ave{i,1} = depth_average(p_i{i,1}',N(i)-1,n_times(i));
    pc_ave{i,1} = depth_average(p_c{i,1}',N(i),n_times(i));
    up_ave{i,1} = depth_average(u_p{i,1}',N(i),n_times(i));
    phi_ave{i,1} = phi_c+depth_average(phi_hat{i,1}',N(i),n_times(i));
    denom_ave{i,1} = depth_average((phi_c(i)./(phi{i,1}')-1),N(i),n_times(i));
end
%% 

nfig = 1;
% axes('ColorOrder',brewermap(nfig+3,'*Purples'))
C = brewermap(nfig+2, 'PuRd');
set(0, 'DefaultAxesColorOrder', C(3:end,:))
% subplot(2,1,1);

hold on

% SetPaperSize(10,10);
plot_times = [51,300];

for i= 1:sim_num
    for j=linspace(1,nfig,nfig)
        k = plot_times(i)+(j-1)*100;
        t_val = (k-1)*t_step;
        plot(phi{i,1}(k,1:end),z_pe_dl{i,1}(1:end),'DisplayName',"t="+num2str(t_val)); 
        
    % Plots the depth averaged value ver time
    %     plot(linspace(1,1500,1500),up_ave(1:1500)')
    end
%     plot(p_b_dl{i,1}-p_crit{i,1},z_pe_dl{i,1}(1:end),"k--",'DisplayName',"Critical $p_e$");
end

ylabel("$z/h$");
xlabel('$p_e$');
% box on
% xlim([0, 5])
title("Decrease of $u_p$ at Late Time in Rauter Model");

% PrintFig('Rau_5deg_up_late')

% Old code to make animations
% make_avi=false;
% if (make_avi) 
%     nframes=3001;
%     t_vals = ((0:nframes))*10;
%     Frames=moviein(nframes);
%     figure
%     v = VideoWriter('Rauter_comp_deep_8_5deg.avi');
%     v.FrameRate=50;
%     open(v)
%     ax = gca();
%     for i=1:1000
%         sgtitle("t="+num2str(t_vals(i))); 
%         subplot(2,2,1)
%         plot(old_p_e(i,:), z_pe_dl{1,1});
%         hold on
%         plot(p_e{1,1}(i,:), z_pe_dl{1,1});
% %         plot((-buoyancy(ani_num)-(density_ratio(ani_num)*phi_c(ani_num)+(1-phi_c(ani_num)))*sind(theta(ani_num))/0.342)*(1-z_pe{ani_num,1}), z_pe{1,1},"--");
%         hold off
%         legend("Iverson","Rauter","Location",'northeast','Interpreter','latex');
% %         legend();
%         xlim([0 0.25]);
%         ylim([0 1]);
%         xlabel("$p_e$",'Interpreter','latex');
%         ylabel("$z$",'Interpreter','latex');
% 
%         subplot(2,2,2)
%         plot(phi_c(1)+old_phi_hat(i,:), z_pe_dl{1,1});
%         hold on
%         plot(phi{1,1}(i,:), z_pe_dl{1,1});
%         hold off
%         xlim([0.53,0.63])
%         ylim([0 1])
%         xlabel('$\phi$','Interpreter','latex');
%         ylabel("$z$",'Interpreter','latex');
% 
%         subplot(2,2,3)
%         plot(old_u_f(i,:), z_u_dl{1,1});
%         hold on
%         plot(u_f{1,1}(i,:), z_u_dl{1,1});
%         hold off
%         xlim([0,5])
%         ylim([0 1])
%         xlabel("$u_f$",'Interpreter','latex');
%         ylabel("$z$",'Interpreter','latex');
% 
%         subplot(2,2,4)
%         plot(old_u_p(i,:), z_u_dl{1,1});
%         hold on
%         plot(u_p{1,1}(i,:), z_u_dl{1,1});
%         hold off
%         xlim([0,5])
%         ylim([0 1])
%         xlabel("$u_p$",'Interpreter','latex');
%         ylabel("$z$",'Interpreter','latex');
%         drawnow();
%         writeVideo(v,getframe(gcf))
%     end
%     close all
%     close(v)
% end

function beta_val=beta_fn(phihat)
    phi_c=0.585;
    beta_val = (phi_c + phihat).^2./((1-phi_c-phihat).^3);
end

function mu_val = mu_Iv_fn(Iv)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    reg_param = 10^8;
    phi_c=0.585;
    mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
end

function diff_mat = smooth_diff(X,gap,t_scale)
    diff_mat = zeros(size(X,1)-gap,size(X,2));
    for i=gap+1:size(X,1)
        diff_mat(i,:)=(X(i,:)-X(i-gap,:))/(gap*t_scale);
    end
end