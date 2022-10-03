flux_contours = load('flux_contours.txt');
flux_contours_DA = load('../Iverson_DA/DA_flux_contours.txt');
init_flux = 0.1*linspace(1,50,50);


% semilogy(init_flux,flux_contours_4(1,:), 'DisplayName', '$4$ degree slope')
SetPaperSize(13,8)
C = viridis(5);

semilogy(init_flux,flux_contours(1,:), 'color' ,C(4,:), 'DisplayName', 'No Flux - Full Model')

hold on
semilogy(init_flux,flux_contours_DA(1,:), "--", 'color', C(4,:), 'DisplayName', 'No Flux - DA Model')
% semilogy(init_flux,flux_contours(3,:), 'DisplayName', '10\% of Initial Flux')
semilogy(init_flux,flux_contours(6,:), 'color', C(3,:), 'DisplayName', '25\% - Full')
semilogy(init_flux,flux_contours_DA(6,:), "--", 'color', C(3,:), 'DisplayName', '25\% - DA')
% semilogy(init_flux,flux_contours(7,:), 'DisplayName', '30\% of Initial Flux')
% semilogy(init_flux,flux_contours(9,:), 'DisplayName', '40\% of Initial Flux')
semilogy(init_flux,flux_contours(11,:), 'color', C(2,:),'DisplayName', '50\% - Full')
semilogy(init_flux,flux_contours_DA(11,:), '--', 'color', C(2,:),'DisplayName', '50\% - DA')
semilogy(init_flux,flux_contours(16,:), 'color', C(1,:),'DisplayName', '75\% - Full')
semilogy(init_flux,flux_contours_DA(16,:), '--', 'color', C(1,:),'DisplayName', '75\% - DA')
% title('Decrease of Flux on 4 degree slope for Full and DA Models')
xlabel('Initial Flux')
ylabel('Time')
xlim([0,5])
legend('Location', 'east outside')
figname = "Ive_low_flux_4deg_slow comp.pdf";
exp_graph(gcf,figname)
movefile(figname,"Figures/SecondYearReport")

