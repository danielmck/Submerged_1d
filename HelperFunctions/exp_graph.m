function exp_graph(fig, fname)
    a=findobj(fig,'type','axe');
    m = size(a,1);
    sub_ind = 1:m;
    posx = zeros(1,m);
    posy = zeros(1,m);
    for n = 1:m
        posx(n) = a(n).Position(1);
        posy(n) = a(n).Position(2);
    end
    Ncols = numel(unique(posx));
    Nrows= numel(unique(posy));
    for k=1:m
        ind = k;
        while (ind>1)
            pos1 = get(a(k),'Position');
            pos2 = get(a(sub_ind(ind-1)),'Position');
            if (pos1(2)<pos2(2)-0.01)
                temp = sub_ind(ind-1);
                sub_ind(ind-1) = k;
                sub_ind(ind) = temp;
                ind=ind-1;
            elseif (pos1(2)<pos2(2)+0.01)
                if (pos1(1)<pos2(1))
                    temp = sub_ind(ind-1);
                    sub_ind(ind-1) = k;
                    sub_ind(ind) = temp;
                    ind=ind-1;
                else
                    break
                end
            else
                break
            end
        end
    end
    for j=1:m
        if (m>1)
            subplot(Nrows,Ncols,sub_ind(j))
        end
        a_now = a(j,1);
        xlab=get(get(a_now,'xlabel'),'string');
        figtitle = get(get(a_now,'title'),'string');
        fs = 12;
        H = get(a_now,'Children');
        if size(xlab,2)
            xCol = get(a_now,'XColor');
            xlabel(xlab,'Interpreter','latex','FontSize', fs)
            set(gca,'XColor', xCol)
        end
        Yax = get(gca,'yaxis');
        if (size(Yax,2)==2)
            yyaxis right
            box on
            ylab=get(get(a_now,'ylabel'),'string');
            if size(ylab,2)
                yCol = get(a_now,'YColor');
                ylabel(ylab,'Interpreter','latex','FontSize', fs)
                set(gca,'YColor', yCol)
            else
                set(gca,'YColor', 'none')
            end
            yyaxis left
        end
        box on
        ylab=get(get(a_now,'ylabel'),'string');
        yCol = get(a_now,'YColor');
        ylabel(ylab,'Interpreter','latex','FontSize', fs)
        set(gca,'YColor', yCol)
        if size(figtitle,2)
            title(figtitle,'Interpreter','latex','FontSize', fs)
        end
        set(gca,'TickLabelInterpreter','latex','FontSize', fs)
        if (sum(size(get(a_now,'legend')))> 0)
            leg_loc= get(get(a_now,'legend'),'location');
            legend('location',leg_loc,'Interpreter','latex','FontSize', fs)
        end
        if (sum(size(get(a_now,'colorbar')))> 0)
            for c = get(a_now,'colorbar')
                c.Label.Interpreter = 'latex';
                c.Label.FontSize = fs;
                c.TickLabelInterpreter = 'latex';
            end
        end
    end
    if ~isempty(findobj(fig, 'Type','subplottext'))
        sglab = findobj(fig, 'Type','subplottext');
        sgtitle(sglab.String,'FontSize',fs, 'Interpreter','latex');
    end
%     if ~strcmp(fname(end-4:end),'.png')
%         fname = strcat(fname,".png");
%     end
    exportgraphics(fig,fname,'Resolution',300)
end
    