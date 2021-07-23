function slowing_sep_var
N = 200; % number of discretisation points in z for each of tau, u
h = 4e-3; % layer height (m)
d=1.43e-5; % grain diameter (m)

mu1_I=0.342; % \frac{u+u_d}{\rho_f \phi_c}
mu2_I=0.557; % 
I_0 = 0.069;

mu1_Iv = 0.32;
mu2_Iv = 0.7;
Iv_0 = 0.005;
reg_param = 10^8;

phi_c=0.6; % Volume fraction
eta_f = 0.0010016; % Pa s
g=9.81; % m/s^2
rho_f = 1000; % kg/m^3
rho_p = 2500; % kg/m^3
theta = 9; % deg

alpha = 0.001; % 1/Pa
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
t_step = 10;
lambda = 0.00001;
A = 2*d_dl*sqrt(density_ratio)/eta_f_dl;
ss_val = sqrt(sind(theta)*rho*(1-z_pe_dl)/mu1_I);
B = sind(theta)*rho*(1-z_pe_dl);

% syms y(z)
% [V] = odeToVectorField(diff(y, 2) == -diff(y)^2/y+3*diff(y)/(1-z)-3/8*y/((1-z)^2)+100/2/((1-z)^(3/2)));
% M = matlabFunction(V,'vars', {'z','Y'});
% SetPaperSize(10,10);
% [~,y_sol] = ode45(@sep_var_z,z_pe_dl(1:end-1),[0.5 1]);
% plot(y_sol(:,1)',z_pe_dl(1:end-1))
% ylabel("z/h")
% xlabel("$Z(z)$")
% xlim([0,3]);
% title("Solution to the Separated ODE for $z$ with $\lambda = 10^{-5}$")
% PrintFig('SepVar_z')
cd Results
full_sim = load("dil9deg_slowing_short.txt");
vec = full_sim(1500,201:400)';
% vec = ones(200,1)*-0.05;
time_vals = (0:3500)*t_step;
opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');

[~,vec]=ode15s(@only_phi_slowing,time_vals,vec,opts);

save("phi_only_evo_slowing.txt", 'vec','-ascii')
cd ../

    function dydz = sep_var_z(z,y)
        dy1 = y(2);
        dy2 = -y(2)^2/y(1)+3*y(2)/(1-z)-3/8*y(1)/((1-z)^2)+lambda/2/((1-z)^(3/2));
        dydz = [dy1,dy2]';
    end
    
    function dphidt = only_phi_slowing(t,phi)
        beta = 150*(phi_c+phi').^2.*eta_f_dl./((1-phi_c+phi').^3.*d_dl^2);
        X = phi'./(phi_c+phi');
        p_p = ss_val.^2 + A*X.^2.*(-mu2_I/mu1_I+1).*(B/mu1_I).^(3/2)/I_0+A.^2.*X.^4.*(-mu2_I/mu1_I+1).*(3-4*mu2_I/mu1_I).*(B/(2*mu2_I*I_0)).^2;
        dpedz = [0 diff(p_b_dl-p_p)/dz_dl];
        dphidt = - phi_c*diff(1./beta.*dpedz)/dz_dl;
        dphidt = interp1(z_pe(1:end-1),dphidt,z_pe,'linear','extrap');
        dphidt = dphidt';
    end
end

