function mu_val = mu_Iv_fn(Iv)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    reg_param = 1*10^7;
    phi_c=0.585;
    mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
end