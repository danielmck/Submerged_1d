function phase_plot
    str_end = "alpha_Fr_9deg_water_big_part";
    d = 1e-3;
    theta = 9;
    phase = "Water";
    comp_val = "u";
    mode_num = 2;
    k_val = 1;

    num_unstab = load("Results/nu_"+str_end+".txt");
    phase_mat = load("Results/h_"+comp_val+"_phase"+num2str(mode_num)+"_"+str_end+"_k_"+num2str(k_val)+".txt");
%     phase_mat_u = load("Results/h_u_phase1_"+str_end+"_k_16.txt");
%     phase_mat = mod(phase_mat_p - phase_mat_u,2*pi);
    n_pts = size(num_unstab,2);
    f = figure;
    width = 10;
    height = 10;
    set(f, 'PaperUnits', 'centimeters');
    set(f, 'PaperSize', [width height]);
    Fr_list = linspace(0.005,5,n_pts);
    alpha_list = logspace(-6,log10(5e-4),n_pts);
    %         set(f,'Position',[15 15 width height])
    % contourf(Fr_list,d_list,num_unstab,2)
    ax = gca;
    ax.YAxis.Exponent = 0;
    ytickformat('%5.1e');
%     xlim([0 1])
%     ylim([1e-5 5e-4])
    %         SetPaperSize(10,10)
    colormap(winter)
    xlabel('Froude Number')
    ylabel('Particle Compressibility')

    title(strcat("$\theta = ",num2str(theta),"$, $ d =",num2str(d,'%5.1e'),"$, Fluid = Water"))
    %         legend()
    set(gca,'YScale', 'log')

    hold on
    [Fr1,alpha1] = meshgrid(Fr_list,alpha_list);
    contourf(Fr1,alpha1,phase_mat,80,'LineColor','none','HandleVisibility','off');
    [x1 y1] = contour(Fr1,alpha1,num_unstab,[0.5 0.5],'*k','LineWidth',1.5);
    [x2 y2] = contour(Fr1,alpha1,num_unstab,[1.5 1.5],'*r','LineWidth',1.5);
%     [x3 y3] = contour(Fr1,alpha1,phase_mat,[pi pi],'*g','LineWidth',1.5);
    legend("1st Mode Stability Boundary","2nd Mode Stability Boundary",'location','southeast');
    c = colorbar;
    c.Ticks = [0,pi/2,pi,3*pi/2,2*pi];
    c.Limits = [0 2*pi];
    tlab = {"0","$\frac{\pi}{2}$","$\pi$","$\frac{3\pi}{2}$","$2\pi$"};
    xlabel('Froude Number')
    ylabel('Particle Compressibility')
    title(strcat("$\theta = ",num2str(theta),"$, $ d =",num2str(d,'%5.1e'),"$, Fluid = ",phase,", $k="+num2str(k_val)+"$"))
%     for j=1:(size(c.Ticks,2)-1)
%         tlab{j+1,1} = "$10^{"+num2str(c.Ticks(j))+"}$";
%     end
    c.TickLabels = tlab;
    if (mode_num == 1)
        rank = "1st";
    else
        rank = "2nd";
    end
    c.Label.String = "Relative phase of $h$ and $"+comp_val+"$ perturbations in "+rank+" Eigenmode";
    fig_name = "phase_h_"+comp_val+"_mode"+num2str(mode_num)+"_"+str_end+".png";
    exp_graph(f,fig_name)
    movefile(fig_name, 'Figures/StabilityPlots/ModePhasePlots');
    % contour(Fr1,alpha1,num_unstab,1.5,'*r');
end