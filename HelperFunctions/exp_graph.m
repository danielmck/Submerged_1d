function exp_graph(fig, fname)
    a=findobj(fig,'type','axe');
    temp_a=a;
    m = size(a,1);
    sub_ind = ones(1,m);
    posx = zeros(1,m);
    posy = zeros(1,m);
    for n = 1:m
        posx(n) = a(n).Position(1);
        posy(n) = a(n).Position(2);
    end
    Ncols = numel(unique(posx));
    Nrows= numel(unique(posy));
    xlabs = strings(1,m);
    figtitles = strings(1,m);
    ylabs = strings(1,m);
    xCols = zeros(m,3);
    yCols = zeros(m,3);
    for k=1:m
        a_now = a(k);
        if size(get(get(a_now,'xlabel'),'string'),1)<1
            xlabs(k)="";
        else
            xlabs(k)=get(get(a_now,'xlabel'),'string');
        end
        if size(get(get(a_now,'title'),'string'),1)<1
            figtitles(k)="";
        else
            figtitles(k)=get(get(a_now,'title'),'string');
        end
        if size(get(get(a_now,'ylabel'),'string'),1)<1
            ylabs(k)="";
        else
            ylabs(k)=get(get(a_now,'ylabel'),'string');
        end
%         figtitles(k) = get(get(a_now,'title'),'string');
        xCols(k,:) = get(a_now,'XColor');
%         ylabs(k)=get(get(a_now,'ylabel'),'string');
        yCols(k,:) = get(a_now,'YColor');
        for l=1:m
            pos1 = get(a(k),'Position');
            pos2 = get(a(l),'Position');
            if (l~=k)
                if (pos1(2)<pos2(2)-0.01)
%                     temp = sub_ind(ind-1);
%                     sub_ind(ind-1) = k;
                    sub_ind(k) = sub_ind(k)+1;
%                     ind=ind-1;
                elseif (pos1(2)<pos2(2)+0.01)
                    if (pos1(1)>pos2(1))
%                         temp = sub_ind(ind-1);
%                         sub_ind(ind-1) = k;
%                         sub_ind(ind) = temp;
%                         ind=ind-1;
                        sub_ind(k) = sub_ind(k)+1;
                    end
                end
            end
        end
    end
    for j=1:m
        if (m>1)
            subplot(Nrows,Ncols,sub_ind(j))
        end
        a_now = a(j);
        xlab=get(get(a_now,'xlabel'),'string');
        figtitle = get(get(a_now,'title'),'string');
        fs = 12;
        if size(xlabs(k),2)
            xCol = get(a_now,'XColor');
            xlabel(xlabs(j),'Interpreter','latex','FontSize', fs)
            set(gca,'XColor', xCols(j,:))
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
%         ylab=get(get(a_now,'ylabel'),'string');
%         yCol = get(a_now,'YColor');
        
%         temp_anow = temp_a(7);
%         
        ylabel(ylabs(j),'Interpreter','latex','FontSize', fs)
        set(gca,'YColor', yCols(j,:))
        if size(figtitles(j),2)
            title(figtitles(j),'Interpreter','latex','FontSize', fs)
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

        sgtitle(sglab.String,'FontSize',fs+1, 'Interpreter','latex');
    end
%     if ~strcmp(fname(end-4:end),'.png')
%         fname = strcat(fname,".png");
%     end
    exportgraphics(fig,fname,'Resolution',300)
end
    