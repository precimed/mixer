f = figure;
%subplot(2,2,[3 4]);
%margin = 0.03;
%ax1 = axes('Parent',f,'position',[margin     0.49  0.5-2*margin 0.44]);
%ax2 = axes('Parent',f,'position',[0.5+margin 0.49  0.5-2*margin 0.44]);

slider1 = uicontrol('Parent',f,'Style','slider','Position', [81,100,419,23], 'value',0.9, 'min',0, 'max',1);
slider1text = uicontrol('Parent',f,'Style','text','Position',[10,100,70,23], 'String','h2','BackgroundColor',f.Color);

vo = 30; % vertical offset
slider2 = uicontrol('Parent',f,'Style','slider','Position', [81,100-vo,419,23], 'value',1, 'min',0.5, 'max',1.5);
slider2text = uicontrol('Parent',f,'Style','text','Position',[10,100-vo,70,23], 'String','(sigma0)^2','BackgroundColor',f.Color);

vo = -30;  
slider3 = uicontrol('Parent',f,'Style','slider','Position', [581,100-vo,419,23], 'value',-1, 'min',-7, 'max',0);
slider3text = uicontrol('Parent',f,'Style','text','Position',[510,100-vo,70,23], 'String','lg(pi)','BackgroundColor',f.Color);

vo = 90; 
slider4 = uicontrol('Parent',f,'Style','slider','Position', [81,100-vo,419,23], 'value', 10, 'min',1, 'max',100);
slider4text = uicontrol('Parent',f,'Style','text','Position',[10,100-vo,70,23], 'String','LD','BackgroundColor',f.Color);

vo = -30; 
slider5 = uicontrol('Parent',f,'Style','slider','Position', [81,100-vo,419,23], 'value', 5, 'min',1, 'max',7);
slider5text = uicontrol('Parent',f,'Style','text','Position',[10,100-vo,70,23], 'String','lg(N)','BackgroundColor',f.Color);

vo = -60; 
cbx1 = uicontrol('Parent', f, 'Style', 'checkbox', 'Position',[81,100-vo,23,23], 'value', false);
cbx1text = uicontrol('Parent',f,'Style','text','Position',[10,100-vo,70,23], 'String','log','BackgroundColor',f.Color);
cbx2 = uicontrol('Parent', f, 'Style', 'checkbox', 'Position',[81+100,100-vo,23,23], 'value', true);
cbx2text = uicontrol('Parent',f,'Style','text','Position',[10+100,100-vo,70,23], 'String','convolution','BackgroundColor',f.Color);

nsnp = 1e7;
meanhet = 0.2165;

vo = 60; 
slider6 = uicontrol('Parent',f,'Style','slider','Position', [81,100-vo,419,23], 'value', 0, 'min',-1, 'max',1);
slider6text = uicontrol('Parent',f,'Style','text','Position',[10,100-vo,70,23], 'String','rho0','BackgroundColor',f.Color);

vo = 0; 
slider7 = uicontrol('Parent',f,'Style','slider','Position', [581,100-vo,419,23], 'value',  0, 'min',-1, 'max',1);
slider7text = uicontrol('Parent',f,'Style','text','Position',[510,100-vo,70,23], 'String','rho_beta','BackgroundColor',f.Color);

vo = 30; 
cbx8 = uicontrol('Parent',f,'Style','checkbox','Position', [581,100-vo,419,23], 'value', false);
cbx8text = uicontrol('Parent',f,'Style','text','Position',[510,100-vo,70,23], 'String','Genetic Overlap','BackgroundColor',f.Color);


cb = @(es,ed) GMM_ofrei_bivariate_plot_callback(struct(...
    'opts_struct', struct('sigma_beta', sqrt(slider1.Value / (nsnp * meanhet * 10^slider3.Value)), 'sigma0', slider2.Value, 'pivec', 10^slider3.Value), ...
    'opts_struct2', struct('sigma_beta', [1 1] * sqrt(slider1.Value / (nsnp * meanhet * 10^slider3.Value)), 'sigma0', [1 1] * slider2.Value, 'rho0', -slider6.Value, 'rho_beta', -slider7.Value, 'pivec', [1-cbx8.Value 1-cbx8.Value cbx8.Value] * 10^slider3.Value), ...
    'tld', slider4.Value, ...
    'N', 10^slider5.Value, ...
    'use_logscale', cbx1.Value, ...
    'convolution', cbx2.Value, ...
    'ld_bins', 30));
slider1.Callback = cb;
slider2.Callback = cb;
slider3.Callback = cb;
slider4.Callback = cb;
slider5.Callback = cb;
slider6.Callback = cb;
slider7.Callback = cb;
cbx8.Callback = cb;
cbx1.Callback = cb;
cbx2.Callback = cb;
cb([],[]);