clc
clear

z_max = 42000;
nc    = 126;  % number of cells on vertical distribution
m     = 2;    % stretch coefficient on x-dir

nl    = nc+1; % number of levels(cell boundaries) on vertical distribution
deta  = 1/nc;
eta   = 0:deta:1;
eta   = eta';

x = m * pi * 2 * ( eta - 0.5 );
k = z_max * 2; %  stretch coefficient on y-dir
f = k * ( atan(x) + pi/2 ) / pi;

F_n = size(1,nl);
for i = 2:nl
    F_n(i) = F_n(i-1) + f(i) * deta;
end

F_n = sum(f) * deta;

x_min = -pi * m;
x_max = x;
F_a = k / ( 2 * pi^2 * m ) * ( ( x_max .* atan(x_max) - log( sqrt( 1 + x_max.^2 ) ) ) - ( x_min .* atan(x_min) - log( sqrt( 1 + x_min.^2 ) ) ) ) + k / 2 * eta ;

figure
plot(eta,f,'r')
hold on
plot(eta,F_a,'b')

level = zeros(2,nl);
for i = 1:nl
    level(:,i) = F_a(i);
end

figure
plot(level,'b')