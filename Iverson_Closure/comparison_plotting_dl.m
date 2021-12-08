% opengl firmware
% names = ["Ive_comp_4_deep_9_2_start_early.txt"]; %,"Ive_comp_4_deep_9_1_start_early.txt","Ive_comp_4_deep_9_2_start_early.txt"]; %"Ive_comp_4_deep.txt"];
names = ["Ive_comp_5_deep_custom_time.txt"];
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
    sim_type(k) = record.Type(in_table);
    N(k) = record.N(in_table);
    h(k) = record.h(in_table);
    d(k) = record.d(in_table);
    dz(k) = h/(N-0.5);
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
    elseif (sim_type(i) == "sin")
        p_e{i,1} = vec(:,1:N(i))';
        p_p{i,1} = (p_b{i,1}-p_e{i,1});
        phi{i,1} = vec(:,N(i)+1:2*N(i))';
        u_f{i,1} = vec(:,2*N(i)+1:end)';
        u_p{i,1} = vec(:,2*N(i)+1:end)';
    elseif (sim_type(i) == "ucon")
        p_e{i,1} = vec(:,1:N(i))';
        p_p{i,1} = (p_b{i,1}-p_e{i,1});
        phi{i,1} = vec(:,N(i)+1:end)';
        initial_u = load("no_phi_no_p.txt");
        u_f{i,1} = initial_u(7,1:N(i))'.*(ones(1,11));
        u_p{i,1} = initial_u(7,N(i)+1:end)'.*(ones(1,11));
    end
    
    if (sim_type(i)~="ucon")
        beta_pe{i,1} = 150*(phi_c(i) + phi{i,1}).^2.*eta_f_dl(i)./((1-phi_c(i)-phi{i,1}).^3.*d_dl(i)^2);
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
            phi_Iv{i,1} = -sqrt(abs(Iv{i,1}))*phi_c(i)./(1+sqrt(abs(Iv{i,1})));
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
    mu1=0.342;
    mu2 = 0.557;
    I0 = 0.069;
    M = 3.1159*2.*d_dl(i).*sqrt(density_ratio(i));
%     pp_const = (density_ratio(i)-1)*cosd(theta(i)).*phi_c(i);
%     Fg = (phi_c(i)*density_ratio(i)+(1-phi_c(i)))*sind(theta(i));
%     ss_A = Fg-mu1*pp_const;
%     ss_B = sqrt(pp_const)*M;
%     ss_dupdz = ss_A*(1-z_u{i,1})/(ss_B*sqrt(1-z_u{i,1})+eta_f_dl(i));
%     ss_u_p = ss_A/(3*ss_B^4)*(ss_B*(-2*ss_B^2*(1-z_u{i,1}).^(3/2)-3*ss_B*eta_f_dl(i)*z_u{i,1}-6*eta_f_dl(i)^2*(1-z_u{i,1}).^(1/2))+eta_f_dl(i)^3*log(ss_B*(1-z_u{i,1}).^(1/2)+eta_f_dl));
% %     plot(ss_A.*(1-z_u{i,1}),z_u{i,1})
%     a=0.0107;
%     b=-0.0101;
%     p=0.9169;%0.9078;
%     c=0.01128;
% %     plot(-crit_grad*(1-z_u{i,1}),z_u{i,1})
%     approx_up = (-0.25*p*a^2*(1-z_u{i,1}).^4-2/3*a*b*p*(1-z_u{i,1}).^3-b^2*p/2*(1-z_u{i,1}).^2)/(eta_f_dl(i)*phi_c(i)^2)+c;
% %     comp_vec = p_p{i,1}.*(phi{i,1}).^2/(phi_c(i)^2*eta_f_dl(i));
%     approx_dupdz = (a*(1-z_u{i,1})+b).^2.*p.*(1-z_u{i,1})./(phi_c(i)^2*eta_f_dl(i));
%     plot(approx_up,z_pe{i,1}(1:end))
%     test_vec = a*(1-z_u{i,1})+b;
    crit_grad = -(density_ratio(i)*phi_c(i)+(1-phi_c(i)))*sind(theta(i))/0.32;
    crit_pe = p_b{i,1}+crit_grad*(1-z_pe{i,1});
    
%     approx_f = (eta_f_dl(i)+M*sqrt(-crit_grad*(1-z_u{i,1}))).*approx_dupdz/(0.342);

%     erosion_point = get_erosion_time(u_p{i,1}',5e-6);
%     pe_grad = zeros(6001,1);
%     for k=1:6001
%         pe_grad(k,1) = dpdz{i,1}(round(max(erosion_point(k)*N(i)-3,1)),k);
%     end
%     ep_diff = erosion_point(430:2030)-erosion_point(400:2000);
%     ep_grad = zeros(1601,1);
%     ep_grad(1)=(ep_diff(1)+ep_diff(2)+ep_diff(3))/3;
%     ep_grad(2)=(ep_diff(1)+ep_diff(2)+ep_diff(3)+ep_diff(4))/4;
%     for l=3:1599record.Type(in_table)
%         ep_grad(l)=(ep_diff(l-2)+ep_diff(l-1)+ep_diff(l)+ep_diff(l+1)+ep_diff(l+2))/5;
%     end
%     ep_grad(1600)=(ep_diff(1598)+ep_diff(1599)+ep_diff(1600)+ep_diff(1601))/4;
%     ep_grad(1601)=(ep_diff(1599)+ep_diff(1600)+ep_diff(1601))/3;
%     plot_vec = -0.342*((crit_grad+8e-5).*(1-z_pe{i,1}.*ones(1,6001))+p_p{i,1})./(eta_f_dl(i)+M*sqrt(p_p{i,1}));%vertcat(zeros(1,n_times(i)),diff(mu{i,1},1,1))./dz_dl(i).*p_p{i,1};
    
%     approx_pe = p_b{i,1}+(crit_grad+4e-5*(2.5*0.6+0.4)).*(1-z_pe{i,1})+approx_f;
    % Plots initial profile  
%     subplot(sim_num,1,i);
    % Plots profile at other times
%     
    
%     plot(linspace(40,200,1601),dupdt{i,1}(155,400:2000))

%     p_p_deriv = (vertcat(zeros(1,n_times(i)),diff(p_p{i,1},1,1))./dz_dl(i).*mu{i,1});
%     mu_deriv = (vertcat(zeros(1,n_times(i)),diff(mu{i,1},1,1))./dz_dl(i).*p_p{i,1});
%     crit_phi = get_critical_value(u_p{i,1}',p_e{i,1}',false,1e-5);
    
    
%     pp_slowing = 1/(2*mu1^2)*(2*Fg*mu1*(1-z_pe{i,1})-2*M*dupdz{i,1}.*sqrt(Fg*mu1.*(1-z_pe{i,1})));
    

%     SetPaperSize(10,8);
%     X = phi{i,1}./(phi{i,1}+phi_c(i));
%     A = 2*d_dl*sqrt(density_ratio(i))/eta_f_dl;
%     eps = A.*X.^2;
%     B = sind(theta(i))*rho(i)*(1-z_pe{i,1});
%     cube = A*X.^2.*(mu2.*p_p{i,1}.^(3/2)-B.*sqrt(p_p{i,1}))+I0.*(mu1.*p_p{i,1}-B);
%     c1 = mu2;
%     c2 = I0.*mu1;
%     c3 = -B;
%     c4 = -I0.*B;
%     p = (3.*c1.*c3.*eps.^2-c2.^2)./(3.*(c1.*eps).^2); %(-3*A^2*X.^4.*B*mu2-mu1^2*I0^2)./(3*A^2*X.^4*mu2^2);
%     q = (2.*c2.^3-9.*c1.*eps.*c2.*c3.*eps+27.*(c1.*eps).^2.*c4)./(27.*(c1.*eps).^3); %(2*mu1^2*I0^2+9*A^2*X.^4*mu2*mu1*I0.*B-27*A^2*X.^4*mu2^2*I0.*B)./(27*A^3*X.^6*mu2^3);
%     delta = -(4*p.^3+27*q.^2);
% %     root = nthroot(real(-q/2-sqrt(q.^2/4+p.^3/27)),3)+nthroot(real(-q/2+sqrt(q.^2/4+p.^3/27)),3)-c2./(3*c1); %
%     mag = 2*sqrt(-p/3);
%     root1 = mag.*cos(acos(3*q./(p.*mag))/3)-c2./(3*c1.*eps);
%     root2 = mag.*cos(acos(3*q./(p.*mag))/3-2*pi/3)-c2./(3*c1.*eps);
%     root3 = mag.*cos(acos(3*q./(p.*mag))/3-4*pi/3)-c2./(3*c1.*eps);
%     cube2 = A*X.^2.*c1.*sqrt(p_p{i,1}).^3+c2.*p_p{i,1}+A*X.^2.*c3.*sqrt(p_p{i,1})+c4;
%     
%     ss_val = sqrt(sind(theta(i))*rho(i)*(1-z_pe{i,1})/mu1);
%     pert_val1 = A*X.^2.*B*(-mu2/mu1+1)/(2*c2);
%     pert_val2 = A.^2*X.^4.*B.^(3/2)/(2*mu1^(3/2)*I0^2)*(-mu2/mu1+1)*(-3/2*mu2/mu1+1);
%     pert_val1_v2 = A*X.^2.*(B*(-mu2/mu1+1)/(2*c2)+sqrt(B)*eta_f_dl(i)/(2*mu1^(3/2)));
%     pert_val2_v2 = A.^2*X.^4.*(B.^(3/2)/(2*mu1^(3/2)*I0^2)*(-mu2/mu1+1)*(-3/2*mu2/mu1+1)+eta_f_dl(i)/(2*I0^2*mu1^3)*(I0*mu1-mu1-mu2/2)*B+3/8*sqrt(B)*eta_f_dl(i)^2/sqrt(mu1^5));
%     pert_val = pert_val1 + pert_val2;
%     
%     pp_expand = ss_val.^2 + A*X.^2.*(-mu2/mu1+1).*(B/mu1).^(3/2)/I0+A.^2.*X.^4.*(-mu2/mu1+1).*(3-4*mu2/mu1).*(B/(2*mu1*I0)).^2;
%     pe_expand = p_b{i,1}-pp_expand;
%     pe_simple = p_b{i,1}-(ss_val.^2 + A*(phi{i,1}/phi_c(i)).^2.*(-mu2/mu1+1).*(B/mu1).^(3/2)/I0);
%     
%     c1_v2 = mu2+X.^2.*eta_f_dl(i)^2;
%     c2_v2 = I0.*mu1+X.^2.*eta_f_dl(i)^2*I0;
%     p_v2 = (3.*c1_v2.*c3.*eps.^2-c2_v2.^2)./(3.*(c1_v2.*eps).^2); %(-3*A^2*X.^4.*B*mu2-mu1^2*I0^2)./(3*A^2*X.^4*mu2^2);
%     q_v2 = (2.*c2_v2.^3-9.*c1_v2.*eps.*c2_v2.*c3.*eps+27.*(c1_v2.*eps).^2.*c4)./(27.*(c1_v2.*eps).^3);
%     root1_v2 = mag.*cos(acos(3*q_v2./(p_v2.*mag))/3)-c2_v2./(3*c1_v2.*eps);
    
%     phi_ave = depth_average(phi{1,1},N(1),6000);
    
%     phi_only = load("Results/phi_only_evo_slowing.txt");
%     X_phi_only = phi_only'./(phi_only'+phi_c(i));
%     pe_phi_only = p_b{i,1}-ss_val.^2 - A*X_phi_only.^2.*(-mu2/mu1+1).*(B/mu1).^(3/2)/I0-A.^2.*X_phi_only.^4.*(-mu2/mu1+1).*(3-4*mu2/mu1).*(B/(2*mu1*I0)).^2;
    
%     phi_simple = load("Results/phi_only_evo_slowing_simple.txt");
    
%     depo_pos = 6000*(1-get_erosion_time(u_p{i,1},1e-5));
    
%     plot_vec = eta_f_dl(i)*dupdz{i,1}./(phi{i,1}./(phi_c(i)+phi{i,1})).^2;
    t_max = 0;
    plot_vec = u_p{i,1};
    %,sum(t_vals{2,1}<t_max),sum(t_vals{3,1}<t_max),5001];
    
    nfig=1;
%     
%     plot(depo_pos,linspace(1,200,200)/200)
%     SetPaperSize(10,10);
%     subplot(1,2,1);
    C = brewermap(nfig+1, 'PuRd');
    
    plot_names = ["ICs from 9.0 degree slope","ICs from 9.1 degree slope","ICs from 9.2 degree slope","ICs from 13 degree slope"];
    for j=linspace(1,nfig,nfig)
%         t_vals{j,1}(t1(j))
        t_max = 0+(j-1)*0.1;
        t_val = (k-1)*t_step(i);
        t1=[sum(t_vals{1,1}<t_max)+1];
        hold on
        plot(plot_vec(1:end,t1),z_pe{i,1}(1:end),'DisplayName',"Full Profile")
%         plot(crit_pe,z_pe{i,1}(1:end),'DisplayName',"$\phi(I_v)$");
%         legend('Location', "best",'UserData', 14);
        ylabel("$z/h$");
        xlabel('$u_p$');
        box on
%         xlim([0,0.05])
%         plot(diffusion_term{i,1}(1:end,j),z_pe{i,1}(1:end),'DisplayName',"Diffusion Term")
    end
%     plot((density_ratio(i)-1)*phi_c(i)*sind(theta(i))/tand(13)*(1-z_pe{i,1}),z_pe{i,1}(1:end),'DisplayName',"Analytic Min Value")
%     plot(crit_pe,z_pe{i,1},"k--",'DisplayName',"Critical $p_e$")
%     legend('Location', "best",'UserData', 12);
%     ylabel("z/h");
%     xlabel('$\frac{\partial p_e}{\partial t}$');
% % % 
% % % %     
    da_vals = load("../Iverson_DA/DA_Results/Ive_da_5_deep_times.txt");
    da_times = da_vals(:,1);
    da_phi = da_vals(:,3);
    da_u = da_vals(:,4);
    da_pe = da_vals(:,5)-cosd(4);
    da_Iv = 3.*eta_f_dl.*da_u./((rho-1).*cosd(theta)-da_pe);
    da_phi_Iv = phi_c./(1+sqrt(da_Iv));
%     plot(da_pe(sum(da_times<300)),0,"rx",'DisplayName',"DA Model Basal Value")
    
%     title('Excess Pressure after Medium Timescale on 5 degree slope')

    u_ave1 = depth_average(phi{1,1},N(1),n_times);
%      SetPaperSize(10,10);
%      [peak_pos, peak_ind] = min(dilatancy{i,1},[],1);
%      yyaxis left;
     tmax = 5000;
%      [~,da_peak] = max(da_pe);
%      [~,full_peak] = max(p_e{i,1}(1,:));
     
     plot(t_vals{1,1}(t_vals{1,1}<tmax),p_e{i,1}(1,t_vals{1,1}<tmax),'DisplayName', "Full Model")
     plot(da_times(da_times<tmax),da_pe(da_times<tmax),'DisplayName', "Depth Averaged Model")
%      plot(p_e{i,1}(1,t_vals{1,1}<tmax),u_ave1(t_vals{1,1}<tmax),'DisplayName', "Full Model")
%      plot(da_pe(da_times<tmax),da_phi(da_times<tmax),'DisplayName', "Depth Averaged Model")
%      ylabel("$\Bar{u}$");
%      yyaxis right;
%      plot(linspace(1,5000,5000+1),p_e{i,1}(1,:),'DisplayName', "Peak Magnitude")
     ylabel("Basal $p_e$");
     ax = gca;
%      set(gca,'ytick',[])
%      ax.YAxis(1).Exponent = 0;
%      ax.XAxis.Exponent = 0;
     xlabel('t');
     legend('Location', "best",'UserData', 15);
     title("Decrease of Basal $p_e$ from $9.2$ degree slope ICs")
     
%      xtickformat('%5.1e');
% %     hold off
    pdfname = 'Ive_pe_full_vs_da_9_2_start';    
    
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
make_avi=false;
if (make_avi) 
    nframes=6000;
    t_vals = ((0:nframes))*10;
    Frames=moviein(nframes);
    figure
    v = VideoWriter('dil9_5_dilatant_depo.avi');
    v.FrameRate=50;
    open(v)
    ax = gca();
    for i=2000:6000
        subplot(2,2,1)
        plot(p_e{ani_num,1}(:,i), z_pe{ani_num,1});
        hold on
        plot((-buoyancy(ani_num)-(density_ratio(ani_num)*phi_c(ani_num)+(1-phi_c(ani_num)))*sind(theta(ani_num))/0.342)*(1-z_pe{ani_num,1}), z_pe{1,1},"--");
        hold off
        legend("t="+num2str(t_vals(i-1999),'%.1f'),"Critical $p_e$","Location",'northeast','Interpreter','latex');
%         legend();
        xlim([0 0.06]);
        ylim([0 1]);
        xlabel("$p_e$",'Interpreter','latex');
        ylabel("$z$",'Interpreter','latex');

        subplot(2,2,2)
        plot(phi_c(1)+phi{ani_num,1}(:,i), z_u{ani_num,1});
%         hold on
%         plot(phi{2,1}(:,i+29), z_u{2,1});
%         hold off
        xlim([0.55,0.6])
        ylim([0 1])
        xlabel('$\phi$','Interpreter','latex');
        ylabel("$z$",'Interpreter','latex');

        subplot(2,2,3)
        plot(u_f{ani_num,1}(:,i), z_u{ani_num,1});
%         hold on
%         plot(u_f{2,1}(:,i+29), z_u{2,1});
%         hold off
        xlim([0,0.5])
        ylim([0 1])
        xlabel("$u_f$",'Interpreter','latex');
        ylabel("$z$",'Interpreter','latex');

        subplot(2,2,4)
        plot(u_p{ani_num,1}(:,i), z_u{ani_num,1});
%         hold on
%         plot(u_p{2,1}(:,i+29), z_u{2,1});
%         hold off
        xlim([0,0.5])
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

function mu_val = mu_Iv_fn(Iv)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    phi_c=0.6;
    mu_val = tanh(1e8*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
end