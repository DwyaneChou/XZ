clc
clear

ncfile    = '..\run\output_xz_3.nc';
pic_path  = '.\';
varname   = 'theta';

time_start = 1;
time_end   = 73;

it = time_end;

% x_min = 0;
% x_max = 20000;
% z_min = 0;
% z_max = 6000;

x_min = 0;
x_max = 19200;
z_min = 0;
z_max = 4800;

R2D    = 180/pi;
radius = 6371229;
g      = 9.80616;

LevelList = -12:1:-0.5;

% var = ncread(ncfile,varname);
x  = ncread(ncfile,'x');
z  = ncread(ncfile,'z');
nt = time_end - time_start + 1;

var = ncread(ncfile,varname,[1,1,it],[Inf,Inf,1]);

disp(['Plotting time ',num2str(it),'/',num2str(nt)])

figure%('visible','off')
var_p = var - 300;
% theta
plt = contour(x,z,var_p,'LevelList',LevelList,'LineStyle','-');
xlim([x_min,x_max])
ylim([z_min,z_max])
colormap(jet)

% pcolor(x,z,var)
% shading interp
% xlim([min(min(x)),max(max(x))])
% ylim([min(min(z)),max(max(z))])
% set(gca,'Clim',[300,302])
% colormap(jet)
% colorbar
% 
% % output picture
% title([varname,' at ',num2str((it-1)*history_interval),' second(s)'])
% print(gcf,'-r600','-dpng',[pic_path,'\',varname,'_',num2str(it-1,'%.4d'),'.png']);
