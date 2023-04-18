function [xi_final,y_final] = better_to_hopital(customit,reverse,params,provideit,master_xi,master_y)
% Does not work, an attempt to convert the no pe wave with no viscosity
% into the full wave
% Likely doesn't work as there is no way to find the flow height gradient
% at the critical point.
    [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
    P = (rho-rho_f)/rho;
    if ~exist("customit","var")
        customit = false;
    end
    if ~exist("reverse","var")
        reverse = false;
    end
    if ~exist("provideit","var")
        provideit = false;
    end
    % If there are custom initial parameters, the function may have to call
    % viscous_Iv_bvp_from_master to make the input
    master_name = "non_vis_convert.txt";
    master_file = load("../Results/"+master_name);
    master_xi = master_file(1,:);
    master_y = master_file(2:end,:);
    record = readtable('../Results/wave_record.csv');

    in_table = strcmp(record.Name, master_name);
    wave_type = record.wave_type(in_table);
    theta = record.theta(in_table);
    lambda_in = record.lambda(in_table);
    Fr = record.Fr(in_table);
    tau0 = record.tau0(in_table);
    alpha = record.alpha(in_table);
    d = record.d(in_table);
    
    
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
    chi = (rho_f+3*rho)/(4*rho);
    
    [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0);
    u_eq = Fr*sqrt(g*cosd(theta)*h0);
    phi_eq = phi_c/(1+sqrt(crit_Iv));
    
    crit_pb = rho_f*g*cosd(theta)*h0;
    
    p_tot = rho*g*cosd(theta);
    
    z_scale = h0;
    v_scale = u_eq;
    p_scale = crit_pb;
    t_scale = z_scale/v_scale;

    tau0_dl = tau0/p_scale;
    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    g_dl = g*t_scale/v_scale; 
    
    u_eq_dl = u_eq/v_scale;
    
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale;
    d_dl = d/z_scale;
    
    p_tot_grad_dl = p_tot/p_scale*z_scale;
    rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
    
    pres_h = 0;
    set_h_min = false;
    h_b_pert = 1e-3;
    
    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);

    master_h = master_y(3,:);
    [~,max_ind] = max(master_h);
    
    master_y = master_y(:,1:max_ind);
    master_y = vertcat(master_y(1:3,:),master_y(5:end,:));
    
    master_xi = master_xi(1:max_ind);
    lambda_in = master_xi(end);
    xi = master_xi/lambda_in;
    
    y_in = vertcat(master_y(2:end,:),master_y(1,:),lambda_in*ones(1,max_ind));
%     extr = [diff(master_xi)>0, true];
%     master_y = master_y(:,extr);
%     master_xi = master_xi(:,extr);
%     crit_ind = sum(master_xi<1);
%     crit_ratio = master_y(3,1)/master_y(2,1);
%     xi = horzcat(master_xi(1:crit_ind)*crit_ratio,crit_ratio+(1-crit_ratio)*(master_xi(crit_ind+1:end)-1));
%     y_in = vertcat(master_y(4:8,:),master_y(1:2,:));
    
    tol=1e-4;
    denom_eps = 5e-4;
    h_min = y_in(2,1);
%     lambda_in = y_in(7,1);
    rf=y_in(3,end);

    solInit1=bvpinit(xi,@bvp_guess);
    opts = bvpset('RelTol',tol,'NMax',3000, 'Stats', true);
    solOut = bvp5c(@viscous_syst, @bc_vals,solInit1,opts);
    
    out_vec = vertcat(xi,y);
    filename = "no_vis_better_master.txt";
    save("Results/"+filename,"out_vec","-ascii")
    write_record("Results/wave_record.csv",filename,{"full","water",Fr,theta,alpha,d,tau0})
        
    function guess = bvp_guess(xi)
        % Initial guess function from the ode soln
        guessd = sum(xi<xi)+1;
        if guessd == max(size(xi))
            guess = y_in(:,end);
        else
            gap = xi(guessd+1)-xi(guessd);
            if abs(gap)<1e-8
                guess = y_in(:,guessd);
            else
                guess = y_in(:,guessd)*(xi-xi(guessd))/gap + y_in(:,guessd+1)*(xi(guessd+1)-xi)/gap;
            end
        end  
    end

    function dydxi = viscous_syst(xi,y)
            % The system from the viroulet paper with a different friction
            % law. Need to solve for u_w and Q1 too.
            u_w = y(6);
            lambda = y(7);
            
            Q1 = y(1);
            h = y(2);
            u = (-Q1 + h.*u_w)./h;
            m = y(3);
            phi = y(4)./Q1;
            pb = y(5) + g_dl*cosd(theta)*rho_dl*chi.*h;
            

            zeta = 3/(2*alpha_dl*h) + P/4;
            p_p = p_tot_grad_dl*h-pb;
            Iv = 3*eta_f_dl*abs(u)/h/p_p;
            D = -2/beta_dl/h*(pb-h);

            denom = (h.^3/Fr^2-Q1.^2);
            fric_coeff = p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv);
            fb_val = tand(theta)-fric_coeff-tau0_dl*rho_f/rho/h;
            dhdxi = 1/Fr^2.*h.^3.*(fb_val+P*D*h)./denom;

            R_w3 = -phi*rho_f_dl/rho_dl*D;
            R_w4 = (-P*chi+zeta)*D - 2*3/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));

            dQdxi = -P*D;
            dmdxi = h/(lambda)*u^(1-pres_h);
            dy6dxi = -R_w3;
            dy7dxi = R_w4/(u-u_w);

            dydxi = [dQdxi,dhdxi,dmdxi,dy6dxi,dy7dxi,0,0];
            if abs(denom)<denom_eps
                A_pp_coeff = -mu_Iv_fn(Iv)./(p_tot_grad_dl.*h)+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./p_p;
                A_pb_coeff = -P.*2./beta_dl-A_pp_coeff;
                A_h_coeff = p_p./(p_tot_grad_dl.*h.^2).*mu_Iv_fn(Iv)+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*3*eta_f_dl.*abs(u)./h.^2./p_p+tau0_dl.*rho_f./rho./h.^2-P.*2./beta_dl+A_pp_coeff.*rho_dl.*g_dl.*cosd(theta);

                full_h_coeff = (h.^3./Fr^2.*(g_dl*cosd(theta)*rho_dl.*chi.*A_pb_coeff+A_h_coeff-p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./u.*Q1./h.^2)+3.*h.^2./Fr^2.*(fb_val+P.*D.*h));
                numer_h_term = full_h_coeff.*dhdxi;
                numer_other_term = h.^3./Fr^2.*(dy7dxi.*A_pb_coeff+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./u.*dQdxi./h.^2);
                              
                h_quad_roots = roots([3*h.^2./Fr^2, -full_h_coeff-2.*Q1.*dQdxi, -numer_other_term]);
                if imag(h_quad_roots(1))==0
                    dydxi(2) = min(h_quad_roots(h_quad_roots>0));
                else
                    dydxi(2) = 0;
                end
            end
            
%             dpbdxi = dy7dxi + g_dl*cosd(theta)*rho_dl*chi.*dhdxi;
%             dppdxi = p_tot_grad_dl*dhdxi-dpbdxi;
%             dudxi = -dQdxi./h+Q1/h^2*dhdxi;
%             dIvdxi = 3*eta_f_dl*dudxi/h/p_p-3*eta_f_dl*abs(u)/h^2/p_p*dhdxi-3*eta_f_dl*abs(u)/h/p_p^2*dppdxi;
%             ddenomdxi = (3*h.^2/Fr^2.*dhdxi-2*Q1.*dQdxi);
%             dforcebaldxi = (-dppdxi/(p_tot_grad_dl*h)+p_p/(p_tot_grad_dl*h^2)*dhdxi)*mu_Iv_fn(Iv)-p_p/(p_tot_grad_dl*h)*dmudIv_fn(Iv)*dIvdxi...
%                 +tau0_dl*rho_f/rho/h^2*dhdxi+P*2/beta_dl*(dpbdxi-dhdxi);
%             dnumerdxi = 1/Fr^2.*h.^3.*dforcebaldxi + 3/Fr^2.*h.^2.*(fb_val+P*D*h);
            dydxi = dydxi'*lambda;
        end

        function resid = bc_vals(ya,yb)
            if set_h_min
                h_start = max(ya(1)/ya(6)+h_b_pert,h_min);
                h_stop = (-h_start + sqrt(h_start.^2+8*ya(1)^2./(h_start/Fr^2)))/2;
                len_resid = [ya(2)-h_start, yb(2)-h_stop];
            else
                h_stop = (-ya(2) + sqrt(ya(2).^2+8*ya(1)^2./(ya(2)/Fr^2)))/2;
                len_resid = [ya(7)-lambda_in, yb(2)-h_stop];
            end
            resid = [ya(1)-yb(1), ya(3), yb(3)-rf,...
                (ya(4)-yb(4)), (ya(5)-yb(5)), len_resid];
        end
end