clc
clear

% ncfile    = '..\run\output_xz_2.nc';
% ncfile    = '..\run\output_xz_2_0p28_0p18.nc';
ncfile    = '..\run\output_xz_2_gaussian_src.nc';
pic_path  = '.\';
varname   = 'w';

time_start = 1;
time_end   = 501;

it = time_end;

R2D    = 180/pi;
radius = 6371229;
g      = 9.80616;

if strcmp(varname,'u')
    max_value = 12;
    min_value = 8;
    dc = 0.2;
elseif strcmp(varname,'w')
    max_value = 2;
    min_value = -2;
    dc = 0.05;
end
TextList = min_value:dc:max_value;

% var = ncread(ncfile,varname);
x  = ncread(ncfile,'x');
z  = ncread(ncfile,'z');
nt = time_end - time_start + 1;

var = ncread(ncfile,varname,[1,1,it],[Inf,Inf,1]);

disp(['Plotting time ',num2str(it),'/',num2str(nt)])
% figure('visible','off')
figure
if strcmp(varname,'u')||strcmp(varname,'w')
    % plt = contour(x,z,var,'LevelStep',0.05,'ShowText','off','LineStyle','-');
    plt = contour(x,z,var,'LevelStep',dc,'LineStyle','-','LineWidth',2);
    % plt = contour(x,z,var,'LevelStep',0.05,'ShowText','on','TextList',TextList,'LineStyle','-','LineWidth',2);
    xlim([-10000,10000])
    ylim([0,10000])
    set(gca,'Clim',[min_value,max_value])
end
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
