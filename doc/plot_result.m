clc
clear

ncfile    = '..\run\output_xz_1.nc';
pic_path  = '.\';
varname   = 'rho';

x_min = -75000;
x_max = 75000;
y_min = 0;
y_max = 15000;

time_start = 1;
% time_end   = 300;
time_end   = 50;

it = time_end;

R2D    = 180/pi;
radius = 6371229;
g      = 9.80616;

% var = ncread(ncfile,varname);
x     = ncread(ncfile,'x');
z     = ncread(ncfile,'z');
sqrtG = ncread(ncfile,'sqrtG');
nt = time_end - time_start + 1;

var0 = ncread(ncfile,varname,[1,1,1 ],[Inf,Inf,1]);
var  = ncread(ncfile,varname,[1,1,it],[Inf,Inf,1]);

disp(['Plotting time ',num2str(it),'/',num2str(nt)])
figure%('visible','off')
var_p = var;
% theta
plt = contour(x,z,var_p,'LineStyle','-');
xlim([x_min,x_max])
ylim([y_min,y_max])
colormap(jet)

% Norm error
f_diff = ( var - var0 ).^2;
fr     = var0.^2;
S1     = s_function(f_diff,sqrtG);
S2     = s_function(fr    ,sqrtG);
L2     = sqrt(S1/S2);

% dx=1000,dz=500,L2=0.078258622350089
% dx=500,dz=250,L2=0.075650091340192

% % output picture
% title([varname,' at ',num2str((it-1)*history_interval),' second(s)'])
print(gcf,'-r600','-dpng',[pic_path,'\',varname,'_',num2str(it-1,'%.4d'),'.png']);

function S = s_function(f,A)
S = sum(f.*A) / sum(A);
end