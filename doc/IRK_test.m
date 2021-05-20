clc
clear

run_seconds      = 2;
dx               = 0.025;
dt               = 0.2;
history_interval = 2;
output_path      = 'picture\';
integral_scheme  = 'IRK'; % Choose from 'RK4', 'IRK2', 'IRK';

% For IRK2 only
IRK_stage      = 4;
IRK_residual   = 1.e-12;
max_outer_iter = 500;
max_inner_iter = 500;

% Set DIRK Coef
if IRK_stage==1
    A(1,1) = 0.5;
    c(1,1) = 0.5;
    b(1,1) = 1;
elseif IRK_stage==2
    alpha = (3 + sqrt(3)) / 6;
    A = zeros(2,2);
    A(1,1) = alpha;
    A(2,1) = 1-2*alpha; A(2,2) = alpha;
    c = [alpha, 1-alpha];
    b = [0.5, 0.5];
elseif IRK_stage==3
    syms x
    x=solve(1/6-3/2*x+3*x^2-x^3==0,x);
    lambda = double(x(2));
    A = zeros(3,3);
    A(1,1) = lambda;
    A(2,1) = (1-lambda)/2; A(2,2) = lambda;
    A(3,1) = (-6*lambda^2+16*lambda-1)/4; A(3,2) = (6*lambda^2-20*lambda+5)/4; A(3,3) = lambda;
    c = [lambda,(1+lambda)/2,1];
    b = A(3,:);
% elseif IRK_stage==3
%     A = zeros(3,3);
%     A(1,1) = 1/3;
%     A(2,1) = 1/2; A(2,2) = 1/2;
%     A(3,1) = 3/4; A(3,2) =-1/4; A(3,3) = 1/2;
%     c = [1/3,1,1];
%     b = [3/4,-1/4,1/2];
% elseif IRK_stage==3
%     alpha = 2 * cos(pi/18) / sqrt(3);
%     A = zeros(3,3);
%     A(1,1) = (1+alpha)/2;
%     A(2,1) = -alpha/2; A(2,2) = (1+alpha)/2;
%     A(3,1) = 1+alpha; A(3,2) = -(1+2*alpha); A(3,3) = (1+alpha)/2;
%     c = [(1+alpha)/2,0.5,(1-alpha)/2];
%     b = [1/(6*alpha^2),1-1/(3*alpha^2),1/(6*alpha^2)];
elseif IRK_stage==4
    alpha = 0.2416942608;
    beta  = 0.0604235652;
    eta   = 0.1291528696;
    sigma = 0.5 - beta - eta - alpha;
    A = zeros(4,4);
    A(1,1) =  alpha;
    A(2,1) = -alpha; A(2,2) =  alpha;
    A(3,1) =      0; A(3,2) =1-alpha; A(3,3) = alpha;
    A(4,1) =   beta; A(4,2) =    eta; A(4,3) = sigma; A(4,4) = alpha;
    c = [alpha,0,1,0.5];
    b = [0,1/6,1/6,2/3];
else
    error(['Unknown IRK_stage']);
end

x_min = -1;
x_max = 1;
x_res = dx;

xi = x_min:x_res:x_max;

n = length(xi)-1;

x = xi(1:n) + dx/2;

u = ones(size(x));
f = zeros(size(x));

% % Square wave
% f(x>-0.4&x<0.4) = 1;

% Sine wave
f = sin( x/(x_max-x_min)*2*pi );

nt  = run_seconds/dt;
ht  = history_interval/dt;

jab = zeros(n,n);
for i = 2:n
    jab(i,i-1) = -u(i-1);
end
jab(1,n) = -u(n);
for i = 1:n-1
    jab(i,i+1) = u(i+1);
end
jab(n,1) = u(1);
jab = - jab / (2*dx);

I = eye(n,n);
for iStage = 1:IRK_stage
    J(:,:,iStage) = I - A(iStage,iStage) * dt * jab;
    iJ(:,:,iStage) = inv(squeeze(J(:,:,iStage)));
end

stat.u      = u;
stat.f      = f;
stat.x      = x;
stat.dx     = x_res;
stat.n      = n;
stat.ht     = ht;
stat.nt     = nt;
stat.dt     = dt;
stat.A      = A;
stat.c      = c;
stat.b      = b;
stat.iJ     = iJ;

plot_result(stat,0,output_path)

total_mass0 = sum(stat.f);

iht = 1;
for it = 1:nt
    if strcmp(integral_scheme,'RK4')
        stat_new = RK4(stat,dt);
    elseif strcmp(integral_scheme,'IRK2')
        stat_new = DIRK(stat,dt,IRK_stage,IRK_residual,max_outer_iter,max_inner_iter);
    elseif strcmp(integral_scheme,'IRK')
        stat_new = IRK(stat,IRK_stage,dt,IRK_residual,max_outer_iter);
    end
    
    stat = stat_new;
    
    if mod(it,ht)==0
        total_mass = sum(stat.f);
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

eps = 1.E-15;

% Calculate derivative on cell center by WENO
stencil_width = 3;
fim = zeros(stencil_width,n);
fip = zeros(stencil_width,n);
fi  = f(1:n);
for i = 1:stencil_width
    fim(i,i+1:n  ) = fi(1:n-i);
    fim(i,1:i    ) = fi(n-i+1:n);
    fip(i,1:n-i  ) = fi(i+1:n);
    fip(i,n-i+1:n) = fi(1:i);
end

tend.f = -( fip(1,:) - fim(1,:) ) / ( 2*dx );

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

function stat = DIRK(stat,dt,IRK_stage,IRK_residual,max_outer_iter,max_inner_iter)
eps = 1.e-15;

A = stat.A;
c = stat.c;
b = stat.b;
n = stat.n;

% Initial condition by Explicit RK
tend1 = spatial_discrete(stat);
stat1 = update_stat(stat,tend1,c(1,1)*dt);
tend2 = spatial_discrete(stat1);
stat2 = update_stat(stat,tend2,c(1,1)*dt);

f1 = stat1.f;
f2 = stat2.f;

% Prepare stat for iteration
q_old = stat.f;

for iStage = 1:IRK_stage
    stat_IRK(iStage) = stat;
    tend_IRK(iStage) = tend1;
end

for iStage = 1:IRK_stage
    qt = 0;
    for i = 1:iStage-1
        qt = qt + dt * A(iStage,i) * tend_IRK(i).f;
    end
    for out_loop = 1:max_outer_iter
        dq = f2 - f1;
        q  = f2;
        
        x0 = reshape(dq,[],1);
        
        delta = sqrt( ( 1 + norm(q) ) * eps ) / norm(dq);
        
        F = calc_F_DIRK(q,q_old,qt,stat2,dt,iStage);
        
        q_pert = q + dq * delta;
        
        F_pert = calc_F_DIRK(q_pert,q_old,qt,stat2,dt,iStage);
        
        Adq = ( F_pert - F ) ./ delta;
        
        r0 = -F - Adq;
        
        beta = norm(r0);
        
        disp(['norm2 r0 = ',num2str(beta),' after loop ',num2str(out_loop),' Stage ',num2str(iStage)]);
        
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
        for j = 1:max_inner_iter
            q_pert = q + V(:,j)' * delta;
            F_pert = calc_F_DIRK(q_pert,q_old,qt,stat2,dt,iStage);
            Adq = ( F_pert - F ) ./ delta;
            R(:,j) = Adq;
            for i = 1:j
                H(i,j) = R(:,j)' * V(:,i);
                R(:,j) = R(:,j) - H(i,j) * V(:,i);
            end
            H(j+1,j) = norm( R(:,j) );
            
            if abs(H(j+1,j)) < 1e-10
                %sprintf('done without residual')
                break;
            else
                V(:,j+1) = R(:,j) ./ H(j+1,j);
                ava_iter = ava_iter + 1;
            end
        end
        
        e1 = zeros(ava_iter,1);
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
        
        f1 = f2;
        
        f2 = f1 + xm';
    end
    stat_IRK(iStage).f = f2;
    tend_IRK(iStage) = spatial_discrete(stat_IRK(iStage));
end

for iStage = 1:IRK_stage
    stat.f = stat.f + dt * b(iStage) * tend_IRK(iStage).f;
end
end

function stat = IRK(stat,IRK_stage,dt,IRK_residual,max_outer_iter)
A           = stat.A;
b           = stat.b;
iJ          = stat.iJ;
q_old       = stat.f;
y           = q_old;
tend        = spatial_discrete(stat);

for iStage = 1:IRK_stage
    stat_IRK(iStage) = stat;
    tend_IRK(iStage) = tend;
end

for iStage = 1:IRK_stage
    y = q_old;
    for i = 1:iStage-1
        y = y + A(iStage,i) * dt * tend_IRK(i).f;
    end
    zm = y + A(iStage,iStage) * dt * tend_IRK(iStage).f;
    
    stat_IRK(iStage).f = zm;
    tend_IRK(iStage) = spatial_discrete(stat_IRK(iStage));
    for iter = 1:max_outer_iter
        dz = squeeze(iJ(:,:,iStage)) * ( -zm + y + A(iStage,iStage) * dt * tend_IRK(iStage).f )';
        dz = dz';
        
        norm_dz = norm(dz);
        if norm_dz<IRK_residual
            break
        end
        stat_IRK(iStage).f = zm + dz;
        tend_IRK(iStage) = spatial_discrete(stat_IRK(iStage));
        zm = stat_IRK(iStage).f;
    end
end

for iStage = 1:IRK_stage
    stat.f = stat.f + b(iStage) * dt * tend_IRK(iStage).f;
end
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

function F = calc_F_DIRK(q,q_old,qt,stat,dt,iStage)
A = stat.A;
stat.f = q;
tend = spatial_discrete(stat);

F = q - q_old - qt - dt * A(iStage,iStage) * tend.f;
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