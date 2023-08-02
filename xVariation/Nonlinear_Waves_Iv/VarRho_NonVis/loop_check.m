N=30;
out_vec = zeros(1,N); %zeros(N,N,N,N);
Fr_num = 1.5; %linspace(0.6,5,N);
uw_num = 3.625; %linspace(2.5,3.0,N);
h_crit = linspace(1.05,1.5,N);

% for i=1:N
    Fr = Fr_num;
%     for j=1:N
        uw = uw_num;
%         for k = 1:N
            hc = h_crit;
            Q1 = -0.6492; %hc*(hc^2-uw);
            hm_crit = -Q1/uw;
            hmin_num = linspace(0.5,0.99,N);
            for l = 1:N
                hmin = hmin_num(l);
                umin = Q1/hmin+uw;
                hmax = (-hmin+sqrt(hmin^2+8*Q1^2*Fr^2/hmin))/2;
                umax = Q1/hmax+uw;
                sq_diff = umax/hmax^2-umin/hmin^2;
                out_vec(1,l) = sq_diff;
            end
%         end
%     end
% end
plot(hmin_num,out_vec)