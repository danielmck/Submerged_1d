names = ["Rauter_6_9_deep.txt"];%"rau_5_deep_small_creep.txt"]; %,"rau_5_alt_creep_small.txt","rau_5.txt"]; %,"rau_8_5_no_creep.txt"];
sim_num=size(names,2);
sim_list = cell(1,sim_num);
sim_title = ["$8.5$ degree Slope", "$5$ degree Slope", "Linear Law, $S_0=5$", "Constant $\phi_m$"];

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

% mu1_Iv = [0.32,0.32];
% mu2_Iv = [0.7,0.7];
% Iv_0 = [0.005,0.005];
% 
% reg_param = [10^8,10^8];

% t_step = [10,10];
% buoyancy = -(rho_p-rho_f).*g.*phi_c.*cosd(theta);

v_scale = sqrt(g.*h);
% p_scale = rho_f.*g.*h;
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
        phi_m{i,1} = phi_c(i); %+(phi_rcp(i)-phi_c(i)).*(abs(dupdz{i,1})<shear_lim_dl(i)).*(shear_lim_dl(i)-abs(dupdz{i,1})).^2/shear_lim_dl(i)^2;
    else
        phi_m{i,1} = phi_c(i)+(phi_rcp(i)-phi_c(i)).*(abs(dupdz{i,1})<shear_lim_dl(i)).*(shear_lim_dl(i)-abs(dupdz{i,1})).^creep_type(i)/shear_lim_dl(i)^creep_type(i);
    end
    p_i{i,1} = eta_f_dl(i).*abs(dupdz{i,1})./((phi_m{i,1}./(phi{i,1}(:,1:end-1))-1).^2);

    p_p{i,1} = horzcat(p_c{i,1}(:,1:end-1) + p_i{i,1},zeros(n_times(i),1));
    p_e{i,1} = p_b_dl{i,1}.*ones(n_times(i),1)-p_p{i,1}; 

    Iv{i,1} = eta_f_dl(i).*abs(dupdz{i,1})./(p_p{i,1}(:,1:end-1));
    phi_Iv{i,1} = phi_c(i)./(1+sqrt(abs(Iv{i,1})));
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
    % old_pe_ave = depth_average(old_p_e',N,n_times);

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
% Iv_approx = eta_f_dl.*abs(dupdz)./(p_b_dl(1:end-1));
% phi_approx = (p_p{i,1}*phi_rcp+phi_rlp*a_dl)./(a_dl+p_p{i,1});
% dupdz_approx = (sind(theta)*(density_ratio*phi_c+(1-phi_c))*(1-z_u_dl)-0.32*p_p)*0.005/(eta_f_dl*0.38);
% A = sind(theta)*(density_ratio*phi_c+(1-phi_c))*(1-z_u_dl)./p_p;
% sqrt_Iv_approx = Iv_0/2*(-5/2*phi_c+sqrt((5/2*phi_c)^2-4*(mu1_Iv-A)*(mu2_Iv-mu1_Iv)/Iv_0))/(mu2_Iv-mu1_Iv);
% c1=(a_dl./p_p(:,1:end-1)-Iv+1);
% c2=(-a_dl./p_p(:,1:end-1).*(2*phi_m+phi_rlp)-2*phi_m-phi_rcp+Iv*phi_rcp);
% c3=(2*phi_m*phi_rcp+2*phi_rlp*phi_m.*a_dl./p_p(:,1:end-1)+phi_m.^2.*a_dl./p_p(:,1:end-1)+phi_m.^2);
% c4=-(phi_m.^2.*a_dl./p_p(:,1:end-1)*phi_rlp+phi_rcp*phi_m.^2);
% pressure_cube = phi(:,1:end-1).^3.*c1+phi(:,1:end-1).^2.*c2 + phi(:,1:end-1).*c3 +c4;

% press_lhs = (phi_m-phi(:,1:end-1)).^2.*(phi_rcp-phi(:,1:end-1));
% press_lhs_exp = -phi(:,1:end-1).^3+phi(:,1:end-1).^2.*(2*phi_m+phi_rcp)-phi(:,1:end-1).*(2*phi_m*phi_rcp+phi_m.^2)+phi_rcp*phi_m.^2;
% 
% press_rhs = (phi_m-phi(:,1:end-1)).^2.*(phi(:,1:end-1)-phi_rlp).*a_dl./p_p(:,1:end-1)+Iv.*phi(:,1:end-1).^2.*(phi_rcp-phi(:,1:end-1));
% press_rhs_exp = phi(:,1:end-1).^3.*(a_dl./p_p(:,1:end-1)-Iv)+phi(:,1:end-1).^2.*(Iv*phi_rcp-a_dl./p_p(:,1:end-1).*(2*phi_m+phi_rlp)) ... 
% + phi(:,1:end-1).*((2*phi_rlp*phi_m+phi_m.^2).*a_dl./p_p(:,1:end-1)) ...
% - (phi_m.^2.*a_dl./p_p(:,1:end-1)*phi_rlp);
% p = (3.*c1.*c3.^2-c2.^2)./(3.*(c1).^2); 
% q = (2.*c2.^3-9.*c1.*c2.*c3+27.*(c1).^2.*c4)./(27.*(c1).^3);
% delta = -(4*p.^3+27*q.^2);
% root = nthroot(real(-q/2-sqrt(q.^2/4+p.^3/27)),3)+nthroot(real(-q/2+sqrt(q.^2/4+p.^3/27)),3)-c2./(3*c1); %
% mag = 2*sqrt(-p/3);
% root1 = mag.*cos(acos(3*q./(p.*mag))/3)-c2./(3*c1);
% root2 = mag.*cos(acos(3*q./(p.*mag))/3-2*pi/3)-c2./(3*c1);
% root3 = mag.*cos(acos(3*q./(p.*mag))/3-4*pi/3)-c2./(3*c1);
% cube2 = root1.^3.*c1+root1.^2.*c2+root1.*c3+c4;
nfig = 1;
% axes('ColorOrder',brewermap(nfig+3,'*Purples'))
% subplot(2,1,1);
slowing_u = load("rau_5_deep_no_pe.txt");
hold on
% pe_approx = p_b_dl{i,1}(1:end-1)./mu_Iv{i,1}.*(mu_Iv{i,1} - (density_ratio(i).*phi_c(i)+(1-phi_c(i)))/((density_ratio(i)-1)*phi_c(i)).*tand(theta(i)));
% 
% SetPaperSize(12,12);
plot_times = [3000,300];
% shear_lim_zone = (sum(u_p{i,1}>1e-4,2)-sum(dupdz{i,1}>shear_lim_dl(i),2))/N(i);

for i= 1:sim_num
%     plot(shear_lim_dl(i)*ones(200,1),z_pe_dl{i,1}(1:end),'DisplayName',"$\frac{\partial p_i}{\partial \phi}$");
    for j=linspace(1,nfig,nfig)
        k = plot_times(i)+(j-1)*10;
        t_val = (k-1)*t_step;
%         shear_lim_zone = (sum(dupdz{i,1}<shear_lim_dl(i),2)-sum(u_p{i,1}<1e-4,2))/N(i);
%         plot(shear_lim_zone)
        plot(p_e{i,1}(k,1:end),z_pe_dl{i,1}(1:end),'DisplayName',"t="+num2str(t_val)); 
%         plot(p_b_dl{i,1}-p_crit{i,1},z_pe_dl{i,1}(1:end),'DisplayName',"$\frac{\partial p_i}{\partial \phi}$");
%         plot(dpidt_gd{i,1}(k,1:end),z_pe_dl{i,1}(1:end-1),'DisplayName',"$\frac{\partial p_i}{\partial \dot\gamma}$");
%         plot(dupdz{i,1}(k,1:end),z_pe_dl{i,1}(1:end-1),'DisplayName',sim_title(i)+ " at t="+num2str(t_val));
%         plot(p_b_dl{i,1}-p_crit{i,1},z_pe_dl{i,1}(1:end),'DisplayName',sim_title(i));
%         plot(slowing_u(250*(j-1)+1,201:end),z_pe_dl{i,1}(1:end),'DisplayName',"t="+num2str(0.1*250*(j-1)));
%         plot(p_b_dl{i,1}-p_crit{i,1}(1:end),z_pe_dl{i,1}(1:end),'DisplayName',"$\frac{\partial p_i}{\partial z}$");
%         plot(smooth_dpidt{i,1}(k,:),z_pe_dl{i,1}(1:end-1),'DisplayName',"$\frac{\partial p_i}{\partial z}$");
%         plot(smooth_dpcdt{i,1}(k,:),z_pe_dl{i,1}(1:end),'DisplayName',"$\frac{\partial p_c}{\partial t}$");
%         plot(-smooth_dpedt{i,1}(k,:),z_pe_dl{i,1}(1:end),'DisplayName',"$\frac{\partial p_p}{\partial t}$");
%         plot(-phi_c(i).*1./beta_pe{i,1}(k,:).*d2pdz2{i,1}(k,:),z_pe_dl{i,1}(1:end),'DisplayName',"From $\frac{\partial p_e}{\partial z}$ variation");
%         plot(-phi_c(i).*[diff(1./beta_u{i,1}(k,:)) 0]/dz_dl(i).*dpdz{i,1}(k,:),z_pe_dl{i,1}(1:end),'DisplayName',"From $\beta$ variation");
    %     plot(-(pi/2)^2*p_e(k,1).*cos(z_pe_dl*pi/2),z_pe_dl(1:end),'DisplayName',"Cosine approximation");

    %     plot(0.07*(2/pi)^2*cos(z_pe_dl*pi/2),z_pe_dl)
%         plot(linspace(0,t_step(i)*20,21),up_ave{i,1}(1:21),'DisplayName',"test")
%         plot(linspace(1,t_step(i)*500,500),p_e{i,1}(1:500,1)','DisplayName',sim_title(i))
    %     plot(linspace(2000,3500,1500),p_e(2001:3500,1)')
    %     plot(linspace(1,1500,1500),pi_ave(1:1500)')
    end
end
% annotation('arrow',[0.571120689655172 0.853448275862069],...
%     [0.535912751677852 0.740492170022371]);
% plot(shear_lim_dl(i)*ones(200,1),z_pe_dl{i,1}(1:end),'--','DisplayName',"$S_0$");
% plot(excess_p_c{i,1},z_pe_dl{i,1}(1:end),"--");
h_leg = legend('Location', "best",'UserData', 4);
% get(h_leg,'Position')
% HeightScaleFactor = 1.2;
% NewHeight = h_leg.Position(4) * HeightScaleFactor;
% h_leg.Position(1) = h_leg.Position(1)*3.5;
% h_leg.Position(2) = h_leg.Position(2) + (NewHeight - h_leg.Position(4));
% h_leg.Position(4) = NewHeight;
% get(h_leg,'Position')
ylabel("$z/h$");
xlabel('$\phi$');
% box on
% xlim([0, 5])
title("Increasing Volume Fraction on $8.5$ degree slope");
% ax = gca;
%      set(gca,'ytick',[])
%      ax.YAxis(1).Exponent = 0;
% ax.XAxis.Exponent = 0;
% xtickformat('%5.1e');
% PrintFig('Rau_8_5deg_phi_time_evo')
% subplot(2,1,2);
% hold on
% for j=linspace(1,nfig,nfig) 
%     plot(phi_c+phi_hat(1+50*(j-1),:),z_pe(1:end));
% end

% cd Results
% old_data = load("dil_8_5deg_deep_comp_new_phic.txt");
% cd ..
% old_p_e = old_data(:,1:200);
% old_phi_hat = old_data(:,201:400);
% old_u_f = old_data(:,401:600);
% old_u_p = old_data(:,601:end);

make_avi=false;
if (make_avi) 
    nframes=3001;
    t_vals = ((0:nframes))*10;
    Frames=moviein(nframes);
    figure
    v = VideoWriter('Rauter_comp_deep_8_5deg.avi');
    v.FrameRate=50;
    open(v)
    ax = gca();
    for i=1:1000
        sgtitle("t="+num2str(t_vals(i))); 
        subplot(2,2,1)
        plot(old_p_e(i,:), z_pe_dl{1,1});
        hold on
        plot(p_e{1,1}(i,:), z_pe_dl{1,1});
%         plot((-buoyancy(ani_num)-(density_ratio(ani_num)*phi_c(ani_num)+(1-phi_c(ani_num)))*sind(theta(ani_num))/0.342)*(1-z_pe{ani_num,1}), z_pe{1,1},"--");
        hold off
        legend("Iverson","Rauter","Location",'northeast','Interpreter','latex');
%         legend();
        xlim([0 0.25]);
        ylim([0 1]);
        xlabel("$p_e$",'Interpreter','latex');
        ylabel("$z$",'Interpreter','latex');

        subplot(2,2,2)
        plot(phi_c(1)+old_phi_hat(i,:), z_pe_dl{1,1});
        hold on
        plot(phi{1,1}(i,:), z_pe_dl{1,1});
        hold off
        xlim([0.53,0.63])
        ylim([0 1])
        xlabel('$\phi$','Interpreter','latex');
        ylabel("$z$",'Interpreter','latex');

        subplot(2,2,3)
        plot(old_u_f(i,:), z_u_dl{1,1});
        hold on
        plot(u_f{1,1}(i,:), z_u_dl{1,1});
        hold off
        xlim([0,5])
        ylim([0 1])
        xlabel("$u_f$",'Interpreter','latex');
        ylabel("$z$",'Interpreter','latex');

        subplot(2,2,4)
        plot(old_u_p(i,:), z_u_dl{1,1});
        hold on
        plot(u_p{1,1}(i,:), z_u_dl{1,1});
        hold off
        xlim([0,5])
        ylim([0 1])
        xlabel("$u_p$",'Interpreter','latex');
        ylabel("$z$",'Interpreter','latex');
        drawnow();
        writeVideo(v,getframe(gcf))
    end
    close all
    close(v)
end

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