function mu_val = mu_I_fn(I)
    mu1_I=0.342; 
    mu2_I=0.557;  
    I_0 = 0.069;
    reg_param = 1*10^7;
    mu_val = tanh(reg_param*I).*(mu1_I+(mu2_I-mu1_I)./(1+I_0./abs(I)));
end