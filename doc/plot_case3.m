clc
clear

ncfile    = '..\run\output_xz_3.nc';
pic_path  = '.\';
varname   = 'theta';

time_start = 1;
time_end   = 151;

it = time_end;

R2D    = 180/pi;
radius = 6371229;
g      = 9.80616;

LevelList = -9:1:0;

% var = ncread(ncfile,varname);
x  = ncread(ncfile,'x');
z  = ncread(ncfile,'z');
nt = time_end - time_start + 1;

var = ncread(ncfile,varname,[1,1,it],[Inf,Inf,1]);

disp(['Plotting time ',num2str(it),'/',num2str(nt)])
% figure('visible','off')

var_p = var - 300;
% theta
plt = contour(x,z,var_p,'LevelList',LevelList,'LineStyle','-');
xlim([0,19200])
ylim([0,4800])
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
