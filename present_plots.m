SetPaperSize(10,10);
points = logspace(-6,2,1000);
semilogx(points,mu_I_fn(points))
ylabel('$\mu(I)$')
xlabel('$I$')
% ax = gca;
% ax.YAxis.Exponent = 0;
% ax.XAxis.Exponent = 0;
% ytickformat('%.1e');
PrintFig('mu_profile')


function mu_val = mu_I_fn(I)
    mu1_I=0.342; % \frac{u+u_d}{\rho_f \phi_c}
    mu2_I=0.557; % 
    I_0 = 0.069;
    reg_param = 10^6;
    mu_val = tanh(reg_param*I).*(mu1_I+(mu2_I-mu1_I)./(1+I_0./abs(I)));
end
    
function beta_val=beta_fn(phihat)
        phi_c=0.6;
        eta_f = 0.0010016;
        d=1.43e-5;
        beta_val = 150*(phi_c + phihat).^2.*eta_f./((1-phi_c-phihat).^3.*d^2);
end
    