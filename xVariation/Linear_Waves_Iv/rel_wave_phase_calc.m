function rel_wave_phase_calc
% Finds the relative phase of the eigenvectors of the A matrix for a range
% of Fr and d or alpha
    str_end = "alpha_Fr_9deg_water_big_part";
    d = 1e-3;
    theta = 9;
    phase = 'Water';
    single_k = 1;
    k=1;
    
    rho_p = 2500;
    
    if strcmp(phase,'Air')
        rho_f = 1;
        eta_f = 1.18e-5; % Pa s
    else
        rho_f = 1000; % kg/m^3
        eta_f = 0.0010016;
    end

    Iv_eq = newt_solve_crit_Iv(theta,rho_p,rho_f);
    
    if (single_k == 0)
        k_unstab = load("Results/k_"+str_end+".txt");
        n_pts = size(k_unstab,1);
        k_txt = '';
    else
        k_txt = "_k_"+num2str(k);
        n_pts = 500;
    end
    
    Fr_list = linspace(0.005,5,n_pts);
    alpha_list = logspace(-6,log10(5e-4),n_pts);
    
    h_u_phase = zeros(n_pts);
    h_p_phase = zeros(n_pts);
    h_phi_phase = zeros(n_pts);
    h_u_phase2 = zeros(n_pts);
    h_p_phase2 = zeros(n_pts);
    h_phi_phase2 = zeros(n_pts);
    
    for i=1:n_pts
        alpha = alpha_list(i);
        for j = 1:n_pts
            Fr = Fr_list(j);
            if (single_k == 0)
                k = k_unstab(i,j);
            end
            A_mat = make_A_mat(k,rho_p,rho_f,theta,eta_f,d,alpha,Fr,Iv_eq);
            [A_eig, phase_mat] = phase_from_A(A_mat,k);
            [~,ind]=maxk(imag(A_eig),2);
            
            h_u_phase(i,j) = mod(phase_mat(ind(1),1)-phase_mat(ind(1),2),2*pi);
            h_p_phase(i,j) = mod(phase_mat(ind(1),1)-phase_mat(ind(1),3),2*pi);
            h_phi_phase(i,j) = mod(phase_mat(ind(1),1)-phase_mat(ind(1),4),2*pi);
            h_u_phase2(i,j) = mod(phase_mat(ind(2),1)-phase_mat(ind(2),2),2*pi);
            h_p_phase2(i,j) = mod(phase_mat(ind(2),1)-phase_mat(ind(2),3),2*pi);
            h_phi_phase2(i,j) = mod(phase_mat(ind(2),1)-phase_mat(ind(2),4),2*pi);
        end
    end
    save(strcat('Results/h_u_phase1_',str_end,k_txt,'.txt'),'h_u_phase','-ascii');
    save(strcat('Results/h_p_phase1_',str_end,k_txt,'.txt'),'h_p_phase','-ascii');
    save(strcat('Results/h_phi_phase1_',str_end,k_txt,'.txt'),'h_phi_phase','-ascii');
    save(strcat('Results/h_u_phase2_',str_end,k_txt,'.txt'),'h_u_phase2','-ascii');
    save(strcat('Results/h_p_phase2_',str_end,k_txt,'.txt'),'h_p_phase2','-ascii');
    save(strcat('Results/h_phi_phase2_',str_end,k_txt,'.txt'),'h_phi_phase2','-ascii');
end

function phase = get_phase(wave, x_val,k)
    found_rt = 0;
    if (size(wave,2) ~= size(x_val,2))
        "Vectors need to be same size"
    else
        for i = 1:(size(wave,2)-1)
            if (wave(i)*wave(i+1)<0)
                grad = (wave(i+1)-wave(i))/(x_val(i+1)-x_val(i));
                root = x_val(i) + wave(i)/grad;
                phase = mod(root*k + (grad > 0)*pi,2*pi);
                found_rt = 1;
                break
            end
        end
    end
    if (~found_rt)
        phase = NaN;
%         "No root found"
    end
end

function [A_eig, phase_mat] = phase_from_A(A_mat,k)
    [A_vec, A_val] = eig(A_mat);
    A_eig = diag(A_val);
    x_val = linspace(0,2*pi/k,100);
    phase_mat = zeros(4);
    for l = 1:4
        wave = real(exp(imag(A_eig(l))).*A_vec(:,l).*exp(1i.*k.*x_val));
%         peak = max(wave,[],2);
%         wave = wave./peak;
        
        for m = 1:4
            phase_mat(l,m) = get_phase(wave(m,:),x_val,k);
        end
    end
end