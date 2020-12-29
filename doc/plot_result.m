clc
clear

ncfile    = '..\run\output_xz_1_250m.nc';
pic_path  = '.\';
varname   = 'rho';

x_min = -75000;
x_max = 75000;
y_min = 0;
y_max = 25000;

time_start = 1;
time_end   = 301;
% time_end   = 51;

it = time_end;

R2D    = 180/pi;
radius = 6371229;
g      = 9.80616;

% var = ncread(ncfile,varname);
x     = ncread(ncfile,'x');
z     = ncread(ncfile,'z');
sqrtG = ncread(ncfile,'sqrtG');
nx = size(x,1);
nz = size(z,1);
nt = time_end - time_start + 1;

var0 = ncread(ncfile,varname,[1,1,1 ],[Inf,Inf,1]);
var  = ncread(ncfile,varname,[1,1,it],[Inf,Inf,1]);

disp(['Plotting time ',num2str(it),'/',num2str(nt)])
figure%('visible','off')
% plt = contour(x,z,var,'LineStyle','-');
plt = pcolor(x,z,var);
xlim([x_min,x_max])
ylim([y_min,y_max])
colormap(jet)

% Norm error
f_diff = ( var - var0 ).^2;
fr     = var0.^2;
S1     = s_function(f_diff,sqrtG);
S2     = s_function(fr    ,sqrtG);
L2     = sqrt(S1/S2);

% dx=1000,dz=500,L2= 0.031939202568543
% dx=500,dz=250,L2= 0.003351975927782
% dx=250,dz=125,L2= 7.369259322361821e-04
% dx=125,dz=125,L2=

diff=sum(sum(abs(var-var0)))/(nx*nz);
% dx = 1000, diff= 5.867174223902470e-05
% dx = 500 , diff= 5.339639190593575e-06
% dx = 250 , diff= 8.786019802160002e-07
% dx = 125 , diff= 

% % output picture
% title([varname,' at ',num2str((it-1)*history_interval),' second(s)'])
print(gcf,'-r600','-dpng',[pic_path,'\',varname,'_',num2str(it-1,'%.4d'),'.png']);

function S = s_function(f,A)
S = sum(sum(f.*A)) ./ sum(sum(A));
end