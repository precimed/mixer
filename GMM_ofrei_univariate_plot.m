f = figure;
ax = axes('Parent',f,'position',[0.13 0.49  0.77 0.44]);

vo = -60;
slider1 = uicontrol('Parent',f,'Style','slider','Position', [81,100-vo,419,23], 'value',0.9, 'min',0, 'max',1);
slider1text = uicontrol('Parent',f,'Style','text','Position',[10,100-vo,70,23], 'String','h2','BackgroundColor',f.Color);

vo = -30; % vertical offset
slider2 = uicontrol('Parent',f,'Style','slider','Position', [81,100-vo,419,23], 'value',1, 'min',0.5, 'max',1.5);
slider2text = uicontrol('Parent',f,'Style','text','Position',[10,100-vo,70,23], 'String','sigma0','BackgroundColor',f.Color);

vo = 00; 
slider3 = uicontrol('Parent',f,'Style','slider','Position', [81,100-vo,419,23], 'value',-1, 'min',-7, 'max',0);
slider3text = uicontrol('Parent',f,'Style','text','Position',[10,100-vo,70,23], 'String','log10(pi)','BackgroundColor',f.Color);

vo = 30; 
slider4 = uicontrol('Parent',f,'Style','slider','Position', [81,100-vo,419,23], 'value', 10, 'min',0, 'max',500);
slider4text = uicontrol('Parent',f,'Style','text','Position',[10,100-vo,70,23], 'String','mean TLD','BackgroundColor',f.Color);

vo = 60; 
slider5 = uicontrol('Parent',f,'Style','slider','Position', [81,100-vo,419,23], 'value', 1.0, 'min',0, 'max',1, 'SliderStep', [1 1] * 0.02);
slider5text = uicontrol('Parent',f,'Style','text','Position',[10,100-vo,70,23], 'String','mean ld r2','BackgroundColor',f.Color);

vo = 90; 
slider6 = uicontrol('Parent',f,'Style','slider','Position', [81,100-vo,419,23], 'value', 5, 'min',1, 'max',7);
slider6text = uicontrol('Parent',f,'Style','text','Position',[10,100-vo,70,23], 'String','GWAS N','BackgroundColor',f.Color);

vo = -90;  show_checkboxes = 'on';
cbx1 = uicontrol('Parent', f, 'Style', 'checkbox', 'Visible', show_checkboxes, 'Position',[81,100-vo,23,23], 'value', false);
cbx4text = uicontrol('Parent',f,'Style','text','Visible', show_checkboxes,'Position',[10,100-vo,70,23], 'String','log','BackgroundColor',f.Color);
cbx2 = uicontrol('Parent', f, 'Style', 'checkbox', 'Visible', show_checkboxes,'Position',[81+100,100-vo,23,23], 'value', true);
cbx4text = uicontrol('Parent',f,'Style','text','Visible', show_checkboxes,'Position',[10+100,100-vo,70,23], 'String','convolution','BackgroundColor',f.Color);

nsnp = 1e7;
meanhet = 0.2165;

cb = @(es,ed) GMM_ofrei_univariate_plot_callback(struct(...
    'opts_struct', struct('sigma_beta', slider1.Value / sqrt(nsnp * meanhet * 10^slider3.Value), 'sigma0', slider2.Value, 'pivec', 10^slider3.Value), ...
    'tld', slider4.Value, ...
    'N', 10^slider6.Value, ...
    'use_logscale', cbx1.Value, ...
    'convolution', cbx2.Value, ...
    'mean_ld_r2', slider5.Value));
slider1.Callback = cb;
slider2.Callback = cb;
slider3.Callback = cb;
slider4.Callback = cb;
slider5.Callback = cb;
slider6.Callback = cb;
cbx1.Callback = cb;
cbx2.Callback = cb;
cb([],[]);