function exp_graph(fig, fname)
    a=findobj(fig,'type','axe');
    xlab=get(get(a,'xlabel'),'string');
    ylab=get(get(a,'ylabel'),'string');
    figtitle = get(get(a,'title'),'string');
    fs = 10;
    H = get(a,'Children');
    if size(xlab,2)
        xlabel(xlab,'Interpreter','latex','FontSize', fs)
    end
    if size(ylab,2)
        ylabel(ylab,'Interpreter','latex','FontSize', fs)
    end
    if size(figtitle,2)
        title(figtitle,'Interpreter','latex','FontSize', fs)
    end
    set(gca,'TickLabelInterpreter','latex','FontSize', fs)
    if (sum(size(get(a,'legend')))> 0)
        leg_loc= get(get(a,'legend'),'location');
        legend('location',leg_loc,'Interpreter','latex','FontSize', fs)
    end
    if (sum(size(get(a,'colorbar')))> 0)
        for c = get(a,'colorbar')
            c.Label.Interpreter = 'latex';
            c.Label.FontSize = fs;
            c.TickLabelInterpreter = 'latex';
        end
    end
%     if ~strcmp(fname(end-4:end),'.png')
%         fname = strcat(fname,".png");
%     end
    exportgraphics(fig,fname,'Resolution',300)
end
    