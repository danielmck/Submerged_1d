function SetPaperSize(width,height)
    set(gcf,'PaperUnits','centimeters');
    left = 0; bottom = 0; 
    myfiguresize = [left, bottom, width, height];
    set(gcf, 'PaperPosition', myfiguresize);

    unis = get(gcf,'units');
    ppos = get(gcf,'paperposition');
    set(gcf,'units',get(gcf,'paperunits'));
    pos = get(gcf,'position');
    pos(3:4) = ppos(3:4);
    set(gcf,'position',pos);
    set(gcf,'units',unis);
end