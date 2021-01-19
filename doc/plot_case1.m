clc
clear

ncfile    = '..\run\output_xz_1.nc';
pic_path  = '.\';
varname   = 'w';

time_start = 1;
time_end   = 201;

it = time_end;

history_interval = 5;

R2D    = 180/pi;
radius = 6371229;
g      = 9.80616;

if strcmp(varname,'u')
    max_value = 10;
    min_value = -10;
    TextList  = [-10,-5,-3,0,3,5,10];
elseif strcmp(varname,'w')
    max_value = 15;
    min_value = -15;
    TextList  = [-10,-9,-5,-3,0,3,5,9,10];
elseif strcmp(varname,'theta')
    max_value = 302;
    min_value = 300;
    TextList  = [0,0.03,0.2,0.4];
    LevelList = zeros(1,12);
    LevelList(1   ) = 0.01;
    LevelList(2   ) = 0.03;
    LevelList(3:12) = 0.2:0.2:2;
end

% var = ncread(ncfile,varname);
x  = ncread(ncfile,'x');
z  = ncread(ncfile,'z');
nt = time_end - time_start + 1;

var = ncread(ncfile,varname,[1,1,it],[Inf,Inf,1]);

disp(['Plotting time ',num2str(it),'/',num2str(nt)])
% figure('visible','off')

if strcmp(varname,'u')||strcmp(varname,'w')
    var_p = var;
    var_n = var;
    var_p(var_p<0)=0;
    var_n(var_n>0)=0;
    
    % Positive
    plt = contour(x,z,var_p,'LevelStep',1,'ShowText','on','TextList',TextList,'LineStyle','-');
    xlim([5000,15000])
    % xlim([min(min(x)),max(max(x))])
    ylim([min(min(z)),max(max(z))])
    hold on
    % Negative
    plt = contour(x,z,var_n,'LevelStep',1,'ShowText','on','TextList',TextList,'LineStyle','--');
    xlim([5000,15000])
    % xlim([min(min(x)),max(max(x))])
    ylim([min(min(z)),max(max(z))])
elseif strcmp(varname,'theta')
    var_p = var - 300;
    % theta
    plt = contour(x,z,var_p,'LevelList',LevelList,'ShowText','on','TextList',TextList,'LineStyle','-');
    xlim([5000,15000])
    % xlim([min(min(x)),max(max(x))])
    ylim([min(min(z)),max(max(z))])
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
