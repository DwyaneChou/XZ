clc
clear

x_min = -2.5;
x_max = 2.5;
z_min = -2.5;
z_max = 1.5;

dx = 1;
dz = 1;

x = x_min:dx:x_max;
z = z_min:dz:z_max;

nx = size(x,2)-1;
nz = size(z,2)-1;

N = nx*nz;

A = zeros(N,N);
j = 0;
for k = 1:nz
    for i = 1:nx
        j = j + 1;
        x1 = x(i  );
        x2 = x(i+1);
        z1 = z(k  );
        z2 = z(k+1);
        A(j,1:N) = calc_rectangle_poly_integration(nx,nz,x1,x2,z1,z2);
    end
end

invA = inv(A);
detA = det(A);

function c = calc_rectangle_poly_integration(nx,ny,x_min,x_max,y_min,y_max)
k = 0;
N = nx*ny;
c = zeros(1,N);
for j = 0:ny-1
    for i = 0:nx-1
        k = k + 1;
        c(k) = ( x_max^(i+1) - x_min^(i+1) ) * ( y_max^(j+1) - y_min^(j+1) ) / ( ( i + 1 ) * ( j + 1 ) );
    end
end
end