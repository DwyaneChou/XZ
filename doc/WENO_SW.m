clc
clear

run_seconds      = 1;
dx               = 1/200;
dt               = 0.05;
history_interval = 0.1;
output_path      = 'picture\';
integral_scheme  = 'IRK2'; % Choose from 'RK4', 'IRK2';

nVar = 2;

% For IRK2 only
IRK_residual   = 1.e-7;
max_outer_iter = 500;
max_inner_iter = 5000;

x_min = -1;
x_max = 1;
x_res = dx;

x = x_min:x_res:x_max;

% % square wave
% u   = zeros(size(x));
% phi = 100. * ones(size(x));
% phi(x>-0.2&x<0.2) = 150;

% sine wave
u   = zeros(size(x));
phi = 50. * sin( x/(x_max-x_min)*2*pi ) + 150 * ones(size(x)); 

n   = length(x);
nt  = run_seconds/dt;
ht  = history_interval/dt;

stat.var(1,:) = phi;
stat.var(2,:) = phi .* u;

stat.nVar   = nVar;
stat.x      = x;
stat.dx     = x_res;
stat.n      = n;
stat.ht     = ht;
stat.nt     = nt;
stat.dt     = dt;

plot_result(stat,0,output_path)

phi = stat.var(1,:);
total_mass0 = sum(phi(1:end-1));

iht = 1;
for it = 1:nt
    if strcmp(integral_scheme,'RK4')
        stat_new = RK4(stat,dt);
    elseif strcmp(integral_scheme,'IRK2')
        stat_new = IRK2(stat,dt,IRK_residual,max_outer_iter,max_inner_iter);
    end
    
    stat = cp_stat(stat_new);
    
    if mod(it,ht)==0
        phi = stat.var(1,:);
        total_mass = sum(phi(1:end-1));
        MCR = (total_mass - total_mass0)/total_mass0;
        disp(['Output result at ',num2str(iht*stat.ht*stat.dt),' second(s)',' MCR = ',num2str(MCR,'%e'),...
              ' Max/Min = ',num2str(max(phi),'%e'),' ',num2str(min(phi),'%e')])
        plot_result(stat,iht,output_path)
        iht = iht + 1;
    end
end

function tend = spatial_discrete(stat)
nVar = stat.nVar;
dx   = stat.dx;
n    = stat.n;

q  = stat.var;
qR = zeros(nVar,n-1);
qL = zeros(nVar,n-1);
for iVar = 1:nVar
    qR(iVar,:) = WENO(q(iVar,:), 1);
    qL(iVar,:) = WENO(q(iVar,:),-1);
end

qR(:,2:n) = qR;
qR(:,1)   = qR(:,n);
qL(:,n)   = qL(:,1);

fR = calc_F(qR);
fL = calc_F(qL);

uL      = qL(2,:) ./ qL(1,:);
phiL    = qL(1,:);
lambdaL = abs( uL ) + sqrt( phiL );

uR      = qR(2,:) ./ qR(1,:);
phiR    = qR(1,:);
lambdaR = abs( uR ) + sqrt( phiR );

f = 0.5 * ( fR + fL - max(lambdaL,lambdaR) .* ( qL - qR ) );

source_term = 0.;

tend.var = - ( f(:,2:n) - f(:,1:n-1) ) / dx + source_term;
tend.var(:,n) = tend.var(:,1);
end

function F = calc_F(q)
phi  = q(1,:);
phiu = q(2,:);
u    = phiu ./ phi;

F(1,:) = phiu;
F(2,:) = phiu.*u + 0.5 * phi.^2;
end

function f_WENO = WENO(f,dir)
n   = length(f);
eps = 1.E-15;
nm1 = n-1;
nm2 = n-2;

% Calculate derivative on cell center by WENO
stencil_width = 2;
nStencil      = 3;
fim = zeros(stencil_width,nm1);
fip = zeros(stencil_width,nm1);
fi  = f(1:nm1);
for i = 1:stencil_width
    fim(i,i+1:nm1    ) = fi(1:nm1-i);
    fim(i,1:i        ) = fi(nm1-i+1:nm1);
    fip(i,1:nm1-i    ) = fi(i+1:nm1);
    fip(i,nm1-i+1:nm1) = fi(1:i);
end

if dir > 0
    fR(1,:) = fim(2,:)/3. - 7./6. * fim(1,:) + 11./6. * fi;
    fR(2,:) =-fim(1,:)/6. + 5./6. * fi       + 1./3.  * fip(1,:);
    fR(3,:) = fi      /3. + 5./6. * fip(1,:) - 1./6.  * fip(2,:);
    
    IS(1,:) = 0.25    * (      fim(2,:) - 4. * fim(1,:) + 3. * fi       ).^2 ...
            + 13./12. * (      fim(2,:) - 2. * fim(1,:) +      fi       ).^2;
    IS(2,:) = 0.25    * (      fim(1,:)                 -      fip(1,:) ).^2 ...
            + 13./12. * (      fim(1,:) - 2. * fi       +      fip(1,:) ).^2;
    IS(3,:) = 0.25    * ( 3. * fi       - 4. * fip(1,:) +      fip(2,:) ).^2 ...
            + 13./12. * (      fi       - 2. * fip(1,:) +      fip(2,:) ).^2;
else
    fL(1,:) = fip(2,:)/3. - 7./6. * fip(1,:) + 11./6. * fi;
    fL(2,:) =-fip(1,:)/6. + 5./6. * fi       + 1./3.  * fim(1,:);
    fL(3,:) = fi      /3. + 5./6. * fim(1,:) - 1./6.  * fim(2,:);
    
    IS(1,:) = 0.25    * (      fip(2,:) - 4. * fip(1,:) + 3. * fi       ).^2 ...
            + 13./12. * (      fip(2,:) - 2. * fip(1,:) +      fi       ).^2;
    IS(2,:) = 0.25    * (      fip(1,:)                 -      fim(1,:) ).^2 ...
            + 13./12. * (      fip(1,:) - 2. * fi       +      fim(1,:) ).^2;
    IS(3,:) = 0.25    * ( 3. * fi       - 4. * fim(1,:) +      fim(2,:) ).^2 ...
            + 13./12. * (      fi       - 2. * fim(1,:) +      fim(2,:) ).^2;
end
    
C(1) = 0.1;
C(2) = 0.6;
C(3) = 0.3;

% % origin WENO weight
% alpha = zeros(nStencil,nm1); % init alpha
% for i = 1:nStencil
%     alpha(i,:) = C(i) ./ ( eps + IS(i,:) ).^2;
% end
% 
% omega = alpha ./ sum(alpha,1);

% WENO_Z weight
tau40 = abs( IS(1,:) - IS(2,:) );
tau41 = abs( IS(2,:) - IS(3,:) );
tau5  = abs( IS(3,:) - IS(1,:) );

alpha = zeros(nStencil,nm1); % init alpha
for i = 1:nStencil
    alpha(i,:) = C(i) .* ( 1. + tau5 ./ ( eps + IS(i,:) ) );
end

omega = alpha ./ sum(alpha,1);

% idx = find(tau40<=min(IS,[],1) & tau41>min(IS,[],1));
% omega(:,idx) = repmat( [1./3.,2./3.,0.]',1,size( idx, 2 ) );
% 
% idx = find(tau40>min(IS,[],1) & tau41<=min(IS,[],1));
% omega(:,idx) = repmat( [0.,2./3.,1./3.]',1,size( idx, 2 ) );

% Calculate tend
if dir > 0
    f_WENO = dot(omega,fR,1);
else
    f_WENO = dot(omega,fL,1);
end

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

tend.var = weights(1).*tend1.var  + weights(2).*tend2.var  + weights(3).*tend3.var  + weights(4).*tend4.var;

stat = update_stat(stat,tend,dt);

end

function stat = IRK2(stat,dt,IRK_residual,max_outer_iter,max_inner_iter)
eps = 1.e-14;

nVar = stat.nVar;
n    = stat.n;
ne   = nVar * n; % Number of element

tend1 = spatial_discrete(stat);
stat1 = update_stat(stat,tend1,0.5*dt);
tend2 = spatial_discrete(stat1);
stat2 = update_stat(stat,tend2,0.5*dt);
stat3 = cp_stat(stat2);

q_old = stat.var;
q_old = reshape(q_old,[],1);
for out_loop = 1:max_outer_iter
    dq = stat2.var - stat1.var;
    dq = reshape(dq,[],1);
    q  = stat2.var;
    q  = reshape(q,[],1);
    
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
    
    V = zeros(ne,max_inner_iter+1);
    R = zeros(ne,max_inner_iter);
    H = zeros(max_inner_iter+1,max_inner_iter);
    
    V(:,1) = r0 / beta;
    
    delta = sqrt( ( 1 + norm(q) ) * eps );
    
    % GMRES(m)
    ava_iter = 0; % available iter
    for jter = 1:max_inner_iter
        q_pert = q + V(:,jter) * delta;
        F_pert = calc_F_IRK2(q_pert,q_old,stat3,dt);
        Adq = ( F_pert - F ) ./ delta;
        R(:,jter) = Adq;
        for iter = 1:jter
            H(iter,jter) = R(:,jter)' * V(:,iter);
            R(:,jter) = R(:,jter) - H(iter,jter) * V(:,iter);
        end
        H(jter+1,jter) = norm( R(:,jter) );
        ava_iter = ava_iter + 1;
        
        if abs(H(jter+1,jter)) < 1e-10
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
    
    stat1.var = stat2.var;
    
    dvar = reshape(xm,2,[]);
    stat2.var = stat1.var + dvar;
end

tend = spatial_discrete(stat2);
stat = update_stat(stat,tend,dt);

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

function F = calc_F_IRK2(q,q_old,stat3,dt)
% stat3.f = 0.5 * ( q + q_old );
% stat3.var = q;
% stat3.f = q_old;

stat3.var = reshape(q,2,[]);

tend3 = spatial_discrete(stat3);
dq = reshape(tend3.var,[],1);
F = ( q - q_old ) / ( dt / 2. ) - dq;
end

function stat = update_stat(stat,tend,dt)
stat.var = stat.var + tend.var * dt;
end

function stat_out = cp_stat(stat_in)
stat_out.var    = stat_in.var;
stat_out.nVar   = stat_in.nVar;
stat_out.x      = stat_in.x;
stat_out.dx     = stat_in.dx;
stat_out.n      = stat_in.n;
stat_out.ht     = stat_in.ht;
stat_out.nt     = stat_in.nt;
stat_out.dt     = stat_in.dt;
end

function plot_result(stat,time_idx,output_path)
x     = stat.x;
x_min = min(x);
x_max = max(x);

var_plot = stat.var(1,:);

figure('Visible','off')
plot(x,var_plot)
xlim([x_min,x_max])
ylim([50,200])
title(['WENO\_SW at ',num2str(time_idx*stat.ht*stat.dt),' second(s)'])
print(gcf,'-r400','-dpng',[output_path,'WENO_SW_',num2str(time_idx,'%.4d'),'.png']);
end