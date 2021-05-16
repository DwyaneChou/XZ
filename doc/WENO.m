clc
clear

run_seconds      = 2;
dx               = 1/60;
dt               = 2/1200;
history_interval = 2;
output_path      = 'picture\';
integral_scheme  = 'RK4'; % Choose from 'RK4', 'IRK2';

% For IRK2 only
IRK_residual   = 1.e-11;
max_outer_iter = 500;
max_inner_iter = 1000;

x_min = -1;
x_max = 1;
x_res = dx;

x = x_min:x_res:x_max;

u = ones(size(x));
f = zeros(size(x));

% Square wave
f(x>-0.4&x<0.4) = 1;

% % Sine wave
% f = sin( x/(x_max-x_min)*2*pi );

n   = length(f);

nt  = run_seconds/dt;
ht  = history_interval/dt;

stat.u      = u;
stat.f      = f;
stat.x      = x;
stat.dx     = x_res;
stat.n      = n;
stat.ht     = ht;
stat.nt     = nt;
stat.dt     = dt;

plot_result(stat,0,output_path)

total_mass0 = sum(stat.f(1:end-1));

iht = 1;
for it = 1:nt
    if strcmp(integral_scheme,'RK4')
        stat_new = RK4(stat,dt);
    elseif strcmp(integral_scheme,'IRK2')
        stat_new = IRK2(stat,dt,IRK_residual,max_outer_iter,max_inner_iter);
    end
    
    stat = cp_stat(stat_new);
    
    if mod(it,ht)==0
        total_mass = sum(stat.f(1:end-1));
        MCR = (total_mass - total_mass0)/total_mass0;
        disp(['Output result at ',num2str(iht*stat.ht*stat.dt),' second(s)',' MCR = ',num2str(MCR,'%e'),...
              ' Max/Min = ',num2str(max(stat.f),'%e'),' ',num2str(min(stat.f),'%e')])
        plot_result(stat,iht,output_path)
        iht = iht + 1;
    end
end

L2 = sqrt(sum((f-stat.f).^2)/sum(f.^2));

function tend = spatial_discrete(stat)
u      = stat.u;
f      = stat.f;
dx     = stat.dx;
n      = stat.n;

eps = 1.E-16;
nm1 = n-1;

% Calculate derivative on cell center by WENO
stencil_width = 3;
fim = zeros(stencil_width,nm1);
fip = zeros(stencil_width,nm1);
fi  = f(1:nm1);
for i = 1:stencil_width
    fim(i,i+1:nm1    ) = fi(1:nm1-i);
    fim(i,1:i        ) = fi(nm1-i+1:nm1);
    fip(i,1:nm1-i    ) = fi(i+1:nm1);
    fip(i,nm1-i+1:nm1) = fi(1:i);
end

fwR(1,:) = fim(2,:)/3. - 7./6. * fim(1,:) + 11./6. * fi;
fwR(2,:) =-fim(1,:)/6. + 5./6. * fi       + 1./3.  * fip(1,:);
fwR(3,:) = fi      /3. + 5./6. * fip(1,:) - 1./6.  * fip(2,:);

IS(1,:) = 0.25    * (      fim(2,:) - 4. * fim(1,:) + 3. * fi       ).^2 ...
        + 13./12. * (      fim(2,:) - 2. * fim(1,:) +      fi       ).^2;
IS(2,:) = 0.25    * (      fim(1,:)                 -      fip(1,:) ).^2 ...
        + 13./12. * (      fim(1,:) - 2. * fi       +      fip(1,:) ).^2;
IS(3,:) = 0.25    * ( 3. * fi       - 4. * fip(1,:) +      fip(2,:) ).^2 ...
        + 13./12. * (      fi       - 2. * fip(1,:) +      fip(2,:) ).^2;
    
C(1) = 0.1;
C(2) = 0.6;
C(3) = 0.3;

% % origin WENO weight
% alpha = zeros(stencil_width,nm1); % init alpha
% for i = 1:stencil_width
%     alpha(i,:) = C(i) ./ ( eps + IS(i,:) ).^2;
% end
% 
% omega = alpha ./ sum(alpha,1);
% % origin WENO weight

% WENO_Z weight
tau40 = abs( IS(1,:) - IS(2,:) );
tau41 = abs( IS(2,:) - IS(3,:) );
tau5  = abs( IS(3,:) - IS(1,:) );

alpha = zeros(stencil_width,nm1); % init alpha
for i = 1:stencil_width
    alpha(i,:) = C(i) .* ( 1. + tau5 ./ ( eps + IS(i,:) ) );
end

omega = alpha ./ sum(alpha,1);

% idx = find(tau40<=min(IS,[],1) & tau41>min(IS,[],1));
% omega(:,idx) = repmat( [1./3.,2./3.,0.]',1,size( idx, 2 ) );
% 
% idx = find(tau40>min(IS,[],1) & tau41<=min(IS,[],1));
% omega(:,idx) = repmat( [0.,2./3.,1./3.]',1,size( idx, 2 ) );
% WENO_Z weight

% Calculate tend
fR_WENO        = dot(omega,fwR,1);
fL_WENO(2:nm1) = fR_WENO(1:nm1-1);
fL_WENO(1    ) = fR_WENO(nm1);

dfdx_WENO = ( fR_WENO - fL_WENO ) / dx;
% END of WENO

tend.f    = -dfdx_WENO;
tend.f(n) = tend.f(1);
tend.f    = u .* tend.f;

% % True solution
% dfdx=cos(stat.x/100.*pi) * pi / 100.;

end

function stat = RK4(stat,dt)
weights = [1/6,1/3,1/3,1/6];

tend1 = spatial_discrete(stat);
stat2 = update_stat(stat,tend1,0.5*dt);

tend2 = spatial_discrete(stat2);
stat3 = update_stat(stat,tend2,0.5*dt);

tend3 = spatial_discrete(stat3);
stat4 = update_stat(stat,tend3,dt);

tend4 = spatial_discrete(stat4);

tend.f  = weights(1).*tend1.f  + weights(2).*tend2.f  + weights(3).*tend3.f  + weights(4).*tend4.f;

stat = update_stat(stat,tend,dt);

end

function stat = IRK2(stat,dt,IRK_residual,max_outer_iter,max_inner_iter)
eps = 1.e-14;

n = stat.n;

tend1 = spatial_discrete(stat);
stat3 = update_stat(stat,tend1,0.5*dt);
tend2 = spatial_discrete(stat3);
stat1 = update_stat(stat,tend1,dt);
stat2 = update_stat(stat,tend2,dt);
stat3 = cp_stat(stat2);

q_old = stat.f;

for out_loop = 1:max_outer_iter
    dq    = stat2.f - stat1.f;
    q     = stat2.f;
    
    x0 = reshape(dq,[],1);
    
    delta = sqrt( ( 1 + norm(q) ) * eps ) / norm(dq);
    
    F = calc_F_IRK2(q,q_old,stat3,dt);
    
    q_pert = q + dq * delta;
    
    F_pert = calc_F_IRK2(q_pert,q_old,stat3,dt);
    
    Adq = ( F_pert - F ) ./ delta;
    
    r0 = -F - Adq;
    
    beta = norm(r0);
    
    disp(['norm2 r0 = ',num2str(beta),' after loop ',num2str(out_loop)]);
    
    if(beta<IRK_residual)
        break
    end
    
    V = zeros(n,max_inner_iter+1);
    R = zeros(n,max_inner_iter);
    H = zeros(max_inner_iter+1,max_inner_iter);
    
    V(:,1) = r0 / beta;
    
    delta = sqrt( ( 1 + norm(q) ) * eps );
    
    % GMRES(m)
    ava_iter = 0; % available iter
    for jter = 1:max_inner_iter
        q_pert = q + V(:,jter)' * delta;
        F_pert = calc_F_IRK2(q_pert,q_old,stat3,dt);
        Adq = ( F_pert - F ) ./ delta;
        R(:,jter) = Adq;
        for iter = 1:jter
            H(iter,jter) = R(:,jter)' * V(:,iter);
            R(:,jter) = R(:,jter) - H(iter,jter) * V(:,iter);
        end
        H(jter+1,jter) = norm( R(:,jter) );
        ava_iter = ava_iter + 1;
        
        if abs(H(jter+1,jter)) < 1e-13
            %sprintf('done without residual')
            break;
        else
            V(:,jter+1) = R(:,jter) / H(jter+1,jter);
        end
    end
    
    e1 = zeros(ava_iter+1,1);
    e1(1) = beta;
    
    [T,bk] = givens( H(1:ava_iter+1,1:ava_iter), e1 );
    
    newT = zeros(ava_iter,ava_iter);
    for i = 1:ava_iter
        for j = 1:i
            newT(j,i) = T(j,i);
        end
    end
    y = newT \ bk;
    
    dx = V(:,1:ava_iter) * y;
    xm = x0 + dx;
    
    stat1.f = stat2.f;
    
    stat2.f = stat1.f + xm';
end

stat.f = stat2.f;

end

function [T,bk] = givens( H,b )
%givens: 通过givens变换化上Hessenborg阵为上三角矩阵
%  化 Hx=b 为 Tx=c
    [~,n] = size(H);
    %提取
    Ht = H(n+1,:);
    H = H(1:n,:);
    b = b(1:n);
    %Rotate Matrix need to recrate every iteration
    %R = eye(n,n);%Rotate Matrix
    for k = 1:n-1
        R = eye(n,n);
        down = (H(k,k)^2+H(k+1,k)^2)^(1/2);
        s = H(k+1,k)/down;
        c = H(k,k)/down;
        R(k:k+1,k:k+1) = [c,s;-s,c];
        H = R*H;
        H(k+1,k) = 0;
        b = R*b;
    end
    T = [H;Ht];
    bk = b;
end

function stat = update_stat(stat,tend,dt)
stat.f  = stat.f + tend.f  * dt;
end

function stat_out = cp_stat(stat_in)
stat_out.u      = stat_in.u;
stat_out.f      = stat_in.f;
stat_out.x      = stat_in.x;
stat_out.dx     = stat_in.dx;
stat_out.n      = stat_in.n;
stat_out.ht     = stat_in.ht;
stat_out.nt     = stat_in.nt;
stat_out.dt     = stat_in.dt;
end

function F = calc_F_IRK2(q,q_old,stat3,dt)
stat3.f = 0.5 * ( q + q_old );
% stat3.f = q;
% stat3.f = q_old;

tend3 = spatial_discrete(stat3);

F = ( q - q_old ) / dt - tend3.f;
end

function plot_result(stat,time_idx,output_path)
x_min = min(stat.x);
x_max = max(stat.x);

figure('Visible','off')
plot(stat.x,stat.f)
xlim([x_min,x_max])
ylim([-1.5,1.5])
title(['WENO at ',num2str(time_idx*stat.ht*stat.dt),' second(s)'])
print(gcf,'-r400','-dpng',[output_path,'WENO_',num2str(time_idx,'%.4d'),'.png']);
end