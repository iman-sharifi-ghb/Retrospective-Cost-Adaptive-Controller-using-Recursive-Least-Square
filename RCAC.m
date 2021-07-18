clc;
clear;
close all

global beta0 T
beta0 = 0;
T     = 0.01;

LR = 0.1;       % Learning Rate
x0 = [1;0];   % Initial Condition
Tf = 100;
t  = 0:T:Tf;
N = Tf/T;

% Initialization
u  = zeros(1, N+1);
e1 = zeros(1, N+1);
e2 = zeros(1, N+1);

ym = zeros(N+1, 1);
x1 = x0(1)*ones(1, N+1);
x2 = x0(2)*ones(1, N+1);
y  = x2;
Outputs = x0;

yd = [2*ones(1, N/5), 1*ones(1,N/5), 0*ones(1, N/5), 1*ones(1,N/5), 2*ones(1,N/5+1)];
% yd = ones(1,N+1);

ep = 0;
ei = 0;
ed = 0;

Nc = 3;
P1  = 100*eye(Nc);
P2  = 100*eye(Nc);
teta1 = repmat([0;0;0],1,N);
teta2 = repmat([0;0;0],1,N);

phi_k = zeros(N, 3);
sigma = +1;
u = zeros(N,1);
E = zeros(3,1);
e = zeros(N+1,1);

% Main Loop
for k = 1:N
   
    Outputs = NonLinDynamic(Outputs, u(k));
    Outputs = Outputs + 0.01*norm(Outputs)*randn(size(Outputs));
    
    x1(k+1) = Outputs(1);
    x2(k+1) = Outputs(2);
    
    y(k+1)  = Outputs(2);
    
    % System Identification
    phi = [x1(k), x2(k), u(k)]';
    [teta1(: ,k+1) , P1] = RLS1(phi, y(k+1), teta1(: ,k), P1, Nc) ;
    ym(k+1) = phi'*teta1(:,k+1);
    
    % Update PID Control Gains        
    phi_k(k, :) = E';
    [teta2(:, k+1), P2] = RLS2(phi_k(k, :), u(k), e(k), teta2(:, k), P2, sigma);    

    % Update u(t)
    e(k+1) = yd(k+1) - ym(k+1);
    ep      = e(k+1) - e(k);
    ei      = T/2*(e(k+1) + e(k));
    if k > 1
        ed  = 1/T*(e(k+1)-2*e(k)+e(k-1));
    end
    E = [ep;ei;ed];
    
    u(k+1)  = u(k) + teta2(:, k+1)'*E;
    u(k+1)
    
end
%% plot Outpots
figure;
plot(t, yd, '--g', 'LineWidth', 1.5)
hold on, grid on
plot(t, y, 'LineWidth', 2)
hold on, 
plot(t, ym, 'LineWidth', 2)
xlabel('Time'), ylabel('Outputs')
legend('desired', 'y', 'ym')
% axis([0 10 -5 5])

figure;
plot(t, teta2(1,:), 'LineWidth', 2)
hold on
plot(t, teta2(2,:), 'LineWidth', 2)
hold on
plot(t, teta2(3,:), 'g', 'LineWidth', 2)
xlabel('Time'), ylabel('Amp')
title('PID Gain Tuning')
legend('K_p', 'K_i', 'K_d')
grid on


