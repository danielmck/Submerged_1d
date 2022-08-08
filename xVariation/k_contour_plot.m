function k_contour_plot
    str_end = "alpha_Fr_18deg_air_big_part";
    str_end2 = str_end+"";
    d = 1e-3;
    theta = 18;
    phase = "Air";

    num_unstab = load("Results/nu_"+str_end+".txt");
    k_unstab = load("Results/k_"+str_end2+".txt");
    n_pts = size(k_unstab,2);
    f = figure;
    width = 10;
    height = 10;
    set(f, 'PaperUnits', 'centimeters');
    set(f, 'PaperSize', [width height]);
    Fr_list = linspace(0.005,5,n_pts);
    d_list = logspace(-6,log10(5e-4),n_pts);
    %         set(f,'Position',[15 15 width height])
    % contourf(Fr_list,d_list,num_unstab,2)
    ax = gca;
    ax.YAxis.Exponent = 0;
    ytickformat('%5.1e');
    %         SetPaperSize(10,10)
    colormap(winter)
    xlabel('Froude Number')
    ylabel('Particle Compressibility')

    title(strcat("$\theta = ",num2str(theta),"$, $ d =",num2str(d,'%5.1e'),"$, Fluid = Water"))
    %         legend()
    set(gca,'YScale', 'log')

    %         set(gca,'TickLabelInterpreter','latex')
    % fig_name = 'UnstabMode_9deg_water_big_part.png';
    %         full_fig = strcat(fig_name,'.pdf');
    %         exportgraphics(f,full_fig,'Resolution',300)
    %         PrintFig(fig_name)
    % exp_graph(f,fig_name)
    % movefile(fig_name, 'Figures/StabilityPlots');
    hold on
    [Fr1,alpha1] = meshgrid(Fr_list,d_list);
    contourf(Fr1,alpha1,log10(k_unstab),20,'LineColor','none','HandleVisibility','off');
    [x1 y1] = contour(Fr1,alpha1,num_unstab,[0.5 0.5],'*k','LineWidth',1.5);
    [x2 y2] = contour(Fr1,alpha1,num_unstab,[1.5 1.5],'*r','LineWidth',1.5);
    legend("1st Mode Stability Boundary","2nd Mode Stability Boundary",'location','southeast');
    c = colorbar;
    c.Ticks = -5:3;
    c.Limits = [-5 3];
    tlab = {};
    xlabel('Froude Number')
    ylabel('Particle Compressibility')
    title(strcat("$\theta = ",num2str(theta),"$, $ d =",num2str(d,'%5.1e'),"$, Fluid = ",phase))
    for j=1:size(c.Ticks,2)
        tlab{j,1} = "$10^{"+num2str(c.Ticks(j))+"}$";
    end
    c.TickLabels = tlab;
    c.Label.String = 'Most Unstable Wavenumber of 1st Mode';
    fig_name = "Unstab_k_"+str_end2+".png";
    exp_graph(f,fig_name)
    movefile(fig_name, 'Figures/StabilityPlots');
    % contour(Fr1,alpha1,num_unstab,1.5,'*r');
end

