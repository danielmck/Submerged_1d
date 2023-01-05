function dmudIv = dmudIv_fn(Iv)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    reg_param = 1*10^7;
    phi_c=0.585;
    dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
    end