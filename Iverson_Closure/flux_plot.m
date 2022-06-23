flux_contours = load('flux_contours.txt');
flux_contours_DA = load('../Iverson_DA/DA_flux_contours.txt');
init_flux = 0.1*linspace(1,50,50);

SetPaperSize(10,12);
% semilogy(init_flux,flux_contours_4(1,:), 'DisplayName', '$4$ degree slope')
C = brewermap(4, 'RdYlGn');

semilogy(init_flux,flux_contours(1,:), 'color' ,C(4,:), 'DisplayName', 'No Flux - Full Model')

hold on
semilogy(init_flux,flux_contours_DA(1,:), "--", 'color', C(4,:), 'DisplayName', 'No Flux - DA Model')
% semilogy(init_flux,flux_contours(3,:), 'DisplayName', '10\% of Initial Flux')
semilogy(init_flux,flux_contours(6,:), 'color', C(3,:), 'DisplayName', '25\% of Initial Flux - Full')
semilogy(init_flux,flux_contours_DA(6,:), "--", 'color', C(3,:), 'DisplayName', '25\% of Initial Flux - DA')
% semilogy(init_flux,flux_contours(7,:), 'DisplayName', '30\% of Initial Flux')
% semilogy(init_flux,flux_contours(9,:), 'DisplayName', '40\% of Initial Flux')
semilogy(init_flux,flux_contours(11,:), 'color', C(2,:),'DisplayName', '50\% of Initial Flux - Full')
semilogy(init_flux,flux_contours_DA(11,:), '--', 'color', C(2,:),'DisplayName', '50\% of Initial Flux - DA')
semilogy(init_flux,flux_contours(16,:), 'color', C(1,:),'DisplayName', '75\% of Initial Flux - Full')
semilogy(init_flux,flux_contours_DA(16,:), '--', 'color', C(1,:),'DisplayName', '75\% of Initial Flux - DA')
title('Decrease of Flux on 4 degree slope for Full and DA Models')
xlabel('Initial Flux')
ylabel('Time')
legend('Location', 'southeast','UserData', 15)
PrintFig('Ive_init_flux_time_4deg_full_DA')