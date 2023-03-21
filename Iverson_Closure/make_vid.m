name = "Ive_5deg_13init.txt"; %"Ive_comp_4_deep_1_flux.txt"];%,"Ive_comp_4_deep_9_2_start.txt""Ive_comp_4_deep_9_1_start.txt","
% names = ["Ive_comp_5_deep_custom_time_IC.txt"]; %,"Ive_comp_5_deep_custom_time_IC.txt"];%
record = readtable('Results/result_record.csv');


in_table = strcmp(record.Name, name);
cd Results
sim = load(name);
cd ../
sim_type = record.sim_type(in_table);
N = record.N(in_table);
h = record.h(in_table);
d = record.d(in_table);
dz = h/(N-0.5);
phi_c = record.phi_c(in_table);
density_ratio = record.rho_r(in_table);
eta_f_dl = record.eta_f(in_table);
theta = record.theta(in_table);
alpha_dl = record.alpha(in_table);
t_step = record.t_step(in_table);
n_times = size(sim,1);


buoyancy = -(density_ratio-1)*phi_c*cosd(theta);
z_pe = linspace(1/(2*N),1,N)';
z_u = linspace(0,1-1/(2*N),N)';
p_b = (density_ratio-1)*phi_c*cosd(theta)*(1-z_pe);


t_vals = sim(:,1);
vec = sim(:,2:end);

    
    % Processes the data differently depending on the simulation type
if (sim_type == "dil")
    p_e = vec(:,1:N)';
    p_p = p_b-p_e;
    phi = phi_c+vec(:,N+1:2*N)';
    u_f = vec(:,2*N+1:3*N)';
    u_p = vec(:,3*N+1:end)';
elseif (sim_type == "pcon")
    p_e = s_frac*p_b;
    p_p = (1-s_frac)*p_b;
    % If looking at initial values for dilatancy sim, need to set phi
    % to zero here
    phi = phi_c+alpha_dl*p_p;
    u_f = vec(:,1:N)';
    u_p = vec(:,N+1:end)';
elseif (sim_type == "pdriv")
    p_e = vec(:,1:N)';
    p_p = (p_b-p_e);
    phi = phi_c+alpha_dl*p_p;
    u_f = vec(:,N+1:2*N)';
    u_p = vec(:,2*N+1:end)';
elseif (sim_type == "sin")
    p_e = vec(:,1:N)';
    p_p = (p_b-p_e);
    phi = phi_c+vec(:,N+1:2*N)';
    u_f = vec(:,2*N+1:end)';
    u_p = vec(:,2*N+1:end)';
elseif (sim_type == "ucon")
    p_e = vec(:,1:N)';
    p_p = (p_b-p_e);
    phi = phi_c+vec(:,N+1:end)';
    initial_u = load("no_phi_no_p.txt");
    u_f = initial_u(7,1:N)'.*(ones(1,11));
    u_p = initial_u(7,N+1:end)'.*(ones(1,11));
end

C = viridis(4);

nframes=1000;

ax = gca();
workingDir = 'fluxvids';
cd fluxvids/images
width = 12;
height = 7;
width_p = ceil(width*37.7952755906);
height_p = ceil(height*37.7952755906);
f = figure;
figsave(f,"img.jpg",[1400 700])
delete img.jpg
for i=1:1000
    fig_name = "img"+num2str(i,'%04d')+".jpg";
    
    % Set a size if desired
    set(f,'Position',[15 15 1000 500])
%     hold on
%         SetPaperSize(12,7);
    if i<500
        t_max = 0.0+(i-1);
        rate = 1;
    else
        t_max = 500.0+(i-500)*10;
        rate = 10;
    end

    subplot(1,2,2)
    t_index = sum(t_vals<t_max)+1;
    plot(p_e(:,t_index), z_pe,'color',C(1,:));
    xlim([0 0.6]);
    ylim([0 1]);
    xlabel("$p_e$",'Interpreter','latex','FontSize',12);
%         ylabel("$z$",'Interpreter','latex');

    subplot(1,2,1)
    t_index = sum(t_vals<t_max)+1;
    plot(u_p(:,t_index), z_pe,'color',C(1,:));
    xlim([0 25])
    ylim([0 1])
    xlabel("$u$",'Interpreter','latex','FontSize',12);
    ylabel("$z$",'Interpreter','latex','FontSize',12);

    sgtitle("$\theta = "+num2str(theta)+"^{\circ}$, $t="+num2str(t_max)+"$, Rate = $\times"+num2str(rate)+"$",'Interpreter','latex')      

%     set(f,'Position',[15 15 1200 700])
    im_len = 0;
    im_hgt = 0;
    if (i == 1)
        exportgraphics(f,fig_name,"Resolution",300)
        im = imread(fig_name);
        req_len = size(im,1);
        req_hgt = size(im,2);
    else
        while (im_len~=req_len) || (im_hgt~=req_hgt)
            exportgraphics(f,fig_name,"Resolution",300)
            im = imread(fig_name);
            im_len = size(im,1);
            im_hgt = size(im,2);
        end
    end
%     figsave(f,"img"+num2str(i,'%04d')+".jpg",[1200 700])
    clf;
end
close(f)
%     close all
%     close(v)
cd ../../

%%
imageNames = dir(fullfile(workingDir,'images','*.jpg'));
imageNames = {imageNames.name}';

outputVideo = VideoWriter(fullfile(workingDir,'shuttle_out.avi'));
outputVideo.FrameRate = 30;
open(outputVideo)
for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir,'images',imageNames{ii}));
   writeVideo(outputVideo,img)
end
close(outputVideo)

function figsave(fig,file,rez,txt,bak)
    %Save figure as image, with custom resolutions and text scaling.
    % figsave(fig,file,rez,txt,bak)
    %
    %Example:
    % clf,text(0.1,0.5,{'This text should be';'50 pixels high';'and the image';'900W x 600H pix'},'FontSize',50)
    % figsave(gcf,'Figure.jpg',[900 600])
    if nargin<1 || isempty(fig),  fig  = gcf; end          %figure handle
    if nargin<2 || isempty(file), file = 'Figure.jpg'; end %output file name
    if nargin<3 || isempty(rez),  rez  = [900 600]; end    %resolution [horizontal vertical]
    if nargin<4 || isempty(txt),  txt  = 1; end            %text scale factor
    if nargin<5 || isempty(bak),  bak  = 1; end            %preserve background colour
    set(fig,'PaperPosition',[0 0 rez/(100*txt)],'PaperUnits','inches'); %set paper size (does not affect display)
    if bak
        set(fig,'InvertHardcopy','off'); %preserve background colour
    end
    imwrite(print(gcf,'-RGBImage',['-r' num2str(100*txt,'%f')]),file) %print RGB image to file (slow)
end