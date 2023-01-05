function f = PlotSetup(width,height)
    f = figure;
    set(f, 'PaperUnits', 'centimeters');
    set(f, 'PaperSize', [width height]);
end