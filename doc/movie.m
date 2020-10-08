clc
clear

ncfile    = '..\run\output_xz_1.nc';
pic_path  = '.\picture';
varname   = 'theta';

time_start = 1;
time_end   = 201;

history_interval = 5;

R2D    = 180/pi;
radius = 6371229;
g      = 9.80616;

% var = ncread(ncfile,varname);
x  = ncread(ncfile,'x');
z  = ncread(ncfile,'z');
nt = time_end - time_start + 1;

parfor it = time_start:time_end
    var = ncread(ncfile,varname,[1,1,it],[Inf,Inf,1]);
    
    disp(['Plotting time ',num2str(it),'/',num2str(nt)])
    figure('visible','off')
    
    pcolor(x,z,var)
    shading interp
    xlim([min(min(x)),max(max(x))])
    ylim([min(min(z)),max(max(z))])
    set(gca,'Clim',[300,302])
    colormap(jet)
    colorbar
    
    % output picture
    title([varname,' at ',num2str((it-1)*history_interval),' second(s)'])
    print(gcf,'-r600','-dpng',[pic_path,'\',varname,'_',num2str(it-1,'%.4d'),'.png']);
    
end