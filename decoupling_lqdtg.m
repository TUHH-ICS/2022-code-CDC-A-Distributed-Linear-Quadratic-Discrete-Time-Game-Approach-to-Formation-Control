clc;clear all;close all;
%% code for paper: Decoupling Approach for Solving Linear Quadratic
%  Discrete-Time Games with Application to Consensus Problem 
%  prima.aditya@tuhh.de

n = 2;     %dimensional plane
I = eye(n); 
N = 4;     %amount of agents
r = n*N;   %adjust the size
M = 6;     %number of edges
h = n*M;   %adjust the size

tf   = 10;       %finite time for open-loop horizon
Np   = 0;       %prediction horizon
step = tf*10;
dt   = tf/step; %step size/ sampling time

eta = 5;

D = [-1 -1  0  0 -1  0;
      1  0 -1  0  0 -1;
      0  1  1  1  0  0;
      0  0  0 -1  1  1];
miu = 1;
%% Problem 1
A  = [zeros(r), eye(r);zeros(r), zeros(r)];
F  = eye(2*r) + dt*A;
mu = [1  1  0  0  1  0;
      1  0  1  0  0  1;
      0  1  1  1  0  0;
      0  0  0  1  1  1];
for l=1:N %for amount of robots 
    %state weights
    W(:,:,l)    = diag(mu(l,:));
    L(:,:,l)    = kron((D*W(:,:,l)*D'),I);
    Q(:,:,l)    = dt*[L(:,:,l),zeros(r);
                   zeros(r),L(:,:,l)]; 
    %terminal weights
    Wf(:,:,l)   = eta*W(:,:,l);
    Lf(:,:,l)   = kron((D*Wf(:,:,l)*D'),I);  
    Qf(:,:,l)   = dt*[Lf(:,:,l),zeros(r);
                   zeros(r),Lf(:,:,l)]; 
    %control input weighting matirx
    R(:,:,l)    = dt*eye(n);
    %creating B matrix for each agent
    b(:,:,l)         = zeros(r,n);
    b((l-1)*n+1,1,l) = 1;
    b((l-1)*n+2,2,l) = 1;
    B(:,:,l)         = [zeros(r,n); b(:,:,l)];
    G(:,:,l)         = dt*B(:,:,l) + (dt^2)/2*A*B(:,:,l); 
end
%initialize the position and velocity for each agent
p1  = [1; 5];    p3  = [5; 4];
p2  = [0; 0.5];  p4  = [4; 0];    
v1  = [0; 1];    v3  = [0; 1];
v2  = [0; 1];    v4  = [0; 1];    
po  = [p1;p2;p3;p4];
vo  = [v1;v2;v3;v4];

%% Problem 2
Wn   = miu*eye(M);
Ln   = D*Wn*D';
Qn   = dt*[kron(Ln,I), zeros(r); zeros(r), kron(Ln,I)];
Rn   = dt*eye(r);
Qfn  = eta*Qn;
%matrix coefficients for vs
Bn   = [zeros(r);eye(r)];
%exact discretization
Gn   = dt*Bn + (dt^2)/2*A*Bn;

%% Problem 3
%matrix coefficients for es
At    = [zeros(h), eye(h);zeros(h), zeros(h)];
Ft    = eye(2*h) + dt*At; %Euler method
mut   = N*eye(M); %to define positive weights of each agent -> W
Qt    = dt*kron(kron(eye(2),mut),eye(n));
Rt    = dt*eye(n*M);
Qft   = eta*Qt;
%matrix coefficients for vs
Bt    = [zeros(h);eye(h)];
%exact discretization
Gt    = dt*Bt + (dt^2)/2*At*Bt;
Phi   = -D';
Phia  = kron(Phi,eye(n));
Phiai = pinv(Phia);

for idx=1:Np+1%indexing for prediction horizon
    if idx==1
        z(:,idx)= [po(1:2)-po(3:4); po(1:2)-po(5:6); po(3:4)-po(5:6); po(7:8)-po(5:6); po(1:2)-po(7:8); po(3:4)-po(7:8); ...
                   vo(1:2)-vo(3:4); vo(1:2)-vo(5:6); vo(3:4)-vo(5:6); vo(7:8)-vo(5:6); vo(1:2)-vo(7:8); vo(3:4)-vo(7:8)];
        %have to initialize the position and velocity vector
        x(:,idx)     = [po;vo];
        xn(:,idx)    = [po;vo];
        xhat(:,idx)  = [po;vo];
    else
        idxn = ((idx-1) * 10) + 1;
        z(:,idxn) = [xf{idx-1}(1:2,idxn)- xf{idx-1}(3:4,idxn);
                     xf{idx-1}(1:2,idxn)- xf{idx-1}(5:6,idxn);
                     xf{idx-1}(3:4,idxn)- xf{idx-1}(5:6,idxn);
                     xf{idx-1}(7:8,idxn)- xf{idx-1}(5:6,idxn);
                     xf{idx-1}(1:2,idxn)- xf{idx-1}(7:8,idxn); 
                     xf{idx-1}(3:4,idxn)- xf{idx-1}(7:8,idxn);
                     xf{idx-1}(9:10,idxn)- xf{idx-1}(11:12,idxn); 
                     xf{idx-1}(9:10,idxn)- xf{idx-1}(13:14,idxn);
                     xf{idx-1}(11:12,idxn)-xf{idx-1}(13:14,idxn);
                     xf{idx-1}(15:16,idxn)-xf{idx-1}(13:14,idxn);
                     xf{idx-1}(9:10,idxn)- xf{idx-1}(15:16,idxn);
                     xf{idx-1}(11:12,idxn)-xf{idx-1}(15:16,idxn)];
        %if reced the horizon then x(l) be the initial state
        x(:,idxn)    = xd{idx-1}(:,idxn);
        xn(:,idxn)   = xe{idx-1}(:,idxn);
        xhat(:,idxn) = xf{idx-1}(:,idxn);
    end
    %define the moving finite time
    T = (idx-1)*10+ step;
    %initialize the Riccati solution to be solved backward 
    for l=1:N %for amount of robots 
        P(:,:,l,T+1) = Qf(:,:,l); %initial Riccati for Problem 1
    end
    Pn(:,:,T+1) = Qfn; %initial Riccati for Problem 2
    Pt(:,:,T+1) = Qft; %initial Riccati for Problem 3
    %solve Problem 1
    for k = 1:T
        %-------------------------------------------------------
        for j = T-1:-1:1%backward computation
            for v = 1:N % v stands for vertices
                S(:,:,v) = G(:,:,v)*(inv(R(:,:,v)))*G(:,:,v)';
            end % looping robot for S
            La(:,:,k) = eye(2*r);
            for v = 1:N
                La(:,:,k) = La(:,:,k) + S(:,:,v)*P(:,:,v,j+1);
            end
            for v = 1:N
                P(:,:,v,j) = Q(:,:,v) + (F'*P(:,:,v,j+1) * inv(La(:,:,k))*F);
            end
        end %looping j backward
        for v = 1:N
               u(:,:,v,k) = -inv(R(:,:,v))*G(:,:,v)' * P(:,:,v,k+1) * inv(La(:,:,k)) * F * x(:,k);
        end
        uall(:,k) = [u(:,:,1,k);u(:,:,2,k);u(:,:,3,k);u(:,:,4,k)];
        x(:,k+1)    = F * x(:,k) +  G(:,:,1)*u(:,:,1,k) + G(:,:,2)*u(:,:,2,k) + G(:,:,3)*u(:,:,3,k) + G(:,:,4)*u(:,:,4,k);%    x(:,k+1) = inv(La(:,:,k)) * F * x(:,k);

        %solve Problem 2
        for j = T-1:-1:1%backward computation
            Pn(:,:,j) = Qn + F'*Pn(:,:,j+1)*F - F'*Pn(:,:,j+1)*Gn*inv(Rn+Gn'*Pn(:,:,j+1)*Gn)*Gn'*Pn(:,:,j+1)*F;
        end
        Kn(:,:,k) = -inv(Rn+Gn'*Pn(:,:,k+1)*Gn)*Gn'*Pn(:,:,k+1)*F;
        un(:,k)   = Kn(:,:,k) * xn(:,k);
        xn(:,k+1) = F * xn(:,k) + Gn*un(:,k);

        %solve Problem 3
        for j = T-1:-1:1%backward computation
        Pt(:,:,j) = Qt + Ft'*Pt(:,:,j+1)*Ft - Ft'*Pt(:,:,j+1)*Gt*inv(Rt+Gt'*Pt(:,:,j+1)*Gt)*Gt'*Pt(:,:,j+1)*Ft;
        end
        Kt(:,:,k) = -inv(Rt+Gt'*Pt(:,:,k+1)*Gt)*Gt'*Pt(:,:,k+1)*Ft;
        a(:,k)    = Kt(:,:,k) * z(:,k);
        z(:,k+1)  = Ft * z(:,k) + Gt * a(:,k);
        for v = 1:N
             uhat(:,:,k) = Phiai * a(:,k);
        end
        xhat(:,k+1) = F * xhat(:,k) +G(:,1:2)*uhat(1:2,k) + G(:,3:4)*uhat(3:4,k) + G(:,5:6)*uhat(5:6,k) + G(:,7:8)*uhat(7:8,k) ;
    end
    xd{idx} = x;
    xe{idx} = xn;
    xf{idx} = xhat;
    zf{idx} = z; 
end

Tot = (tf*10)+(Np*10);

%% Regulator section - Figure 2
figure('Name', 'Relative state and control', 'NumberTitle', 'off')
subplot(2,3,1)
plot(1:Tot+1,z(1,:),1:Tot+1,z(3,:),1:Tot+1,z(5,:),1:Tot+1,z(7,:),1:Tot+1,z(9,:),1:Tot+1,z(11,:),'linewidth',1.8)
yline(0,'--k')
legend('$p_x^1 - p_x^2$','$p_x^1 - p_x^3$','$p_x^2 - p_x^3$','$p_x^4 - p_x^3$','$p_x^1 - p_x^4$','$p_x^2 - p_x^4$','fontsize',12,'interpreter','latex')
xlim([1 Tot+1])
xlabel('time-steps','fontsize',12)
ylabel('x-axis','fontsize',12)
title('Relative positions')
grid on
set(gca,'color',[0.9,0.9,0.9]);
subplot(2,3,2)
plot(1:Tot+1,z(13,:),1:Tot+1,z(15,:),1:Tot+1,z(17,:),1:Tot+1,z(19,:),1:Tot+1,z(21,:),1:Tot+1,z(23,:),'linewidth',1.8)
yline(0,'--k')
legend('$v_x^1 - v_x^2$','$v_x^1 - v_x^3$','$v_x^2 - v_x^3$','$v_x^4 - v_x^3$','$v_x^1 - v_x^4$','$v_x^2 - v_x^4$','fontsize',12,'interpreter','latex')
xlim([1 Tot+1])
xlabel('time-steps','fontsize',12)
ylabel('x-axis','fontsize',12)
title('Relative velocities')
grid on
set(gca,'color',[0.9,0.9,0.9]);
subplot(2,3,3)
plot(1:Tot,a(1,:),1:Tot,a(3,:),1:Tot,a(5,:),1:Tot,a(7,:),1:Tot,a(9,:),1:Tot,a(11,:),'linewidth',1.8)
yline(0,'--k')
legend('$u_x^1 - u_x^2$','$u_x^1 - u_x^3$','$u_x^2 - u_x^3$','$u_x^4 - u_x^3$','$u_x^1 - u_x^4$','$u_x^2 - u_x^4$','fontsize',12,'interpreter','latex')
xlim([1 Tot])
xlabel('time-steps','fontsize',12)
ylabel('x-axis','fontsize',12)
title('Relative controls')
grid on
xlim([1 Tot])
set(gca,'color',[0.9,0.9,0.9]);
subplot(2,3,4)
plot(1:Tot+1,z(2,:),1:Tot+1,z(4,:),1:Tot+1,z(6,:),1:Tot+1,z(8,:),1:Tot+1,z(10,:),1:Tot+1,z(12,:),'linewidth',1.8)
yline(0,'--k')
legend('$p_y^1 - p_y^2$','$p_y^1 - p_y^3$','$p_y^2 - p_y^3$','$p_y^4 - p_y^3$','$p_y^1 - p_y^4$','$p_y^2 - p_y^4$','fontsize',12,'interpreter','latex')
xlim([1 Tot+1])
xlabel('time-steps','fontsize',12)
ylabel('y-axis','fontsize',12)
title('Relative positions')
grid on
set(gca,'color',[0.9,0.9,0.9]);
subplot(2,3,5)
plot(1:Tot+1,z(14,:),1:Tot+1,z(16,:),1:Tot+1,z(18,:),1:Tot+1,z(20,:),1:Tot+1,z(22,:),1:Tot+1,z(24,:),'linewidth',1.8)
yline(0,'--k')
legend('$v_y^1 - v_y^2$','$v_y^1 - v_y^3$','$v_y^2 - v_y^3$','$v_y^4 - v_y^3$','$v_y^1 - v_y^4$','$v_y^2 - v_y^4$','fontsize',12,'interpreter','latex')
xlim([1 Tot+1])
xlabel('time-steps','fontsize',12)
ylabel('y-axis','fontsize',12)
title('Relative velocities')
grid on
set(gca,'color',[0.9,0.9,0.9]);
subplot(2,3,6)
plot(1:Tot,a(2,:),1:Tot,a(4,:),1:Tot,a(6,:),1:Tot,a(8,:),1:Tot,a(10,:),1:Tot,a(12,:),'linewidth',1.8)
yline(0,'--k')
legend('$u_y^1 - u_y^2$','$u_y^1 - u_y^3$','$u_y^2 - u_y^3$','$u_y^4 - u_y^3$','$u_y^1 - u_y^4$','$u_y^2 - u_y^4$','fontsize',12,'interpreter','latex')
xlim([1 Tot])
xlabel('time-steps','fontsize',12)
ylabel('y-axis','fontsize',12)
title('Relative controls')
grid on
xlim([1 Tot])
set(gca,'color',[0.9,0.9,0.9]);

% % %relative dynamics of control inputs
% figure('Name', 'Relative acceleration control input', 'NumberTitle', 'off')
% subplot(2,1,1)
% plot(1:Tot,a(1,:),1:Tot,a(3,:),1:Tot,a(5,:),1:Tot,a(7,:),1:Tot,a(9,:),1:Tot,a(11,:),'linewidth',1.8)
% yline(0,'--k')
% legend('$u_x^1 - u_x^2$','$u_x^1 - u_x^3$','$u_x^2 - u_x^3$','$u_x^4 - u_x^3$','$u_x^1 - u_x^4$','$u_x^2 - u_x^4$','fontsize',12,'interpreter','latex')
% xlim([1 Tot])
% xlabel('time-steps','fontsize',12)
% ylabel('Relative control','fontsize',12)
% % title('Relative acceleration control input in x direction')
% grid on
% xlim([1 Tot])
% set(gca,'color',[0.9,0.9,0.9]);
% subplot(2,1,2)
% plot(1:Tot,a(2,:),1:Tot,a(4,:),1:Tot,a(6,:),1:Tot,a(8,:),1:Tot,a(10,:),1:Tot,a(12,:),'linewidth',1.8)
% yline(0,'--k')
% legend('$u_y^1 - u_y^2$','$u_y^1 - u_y^3$','$u_y^2 - u_y^3$','$u_y^4 - u_y^3$','$u_y^1 - u_y^4$','$u_y^2 - u_y^4$','fontsize',12,'interpreter','latex')
% xlim([1 Tot])
% xlabel('time-steps','fontsize',12)
% ylabel('Relative control','fontsize',12)
% % title('Relative acceleration control input in y direction')
% grid on
% xlim([1 Tot])
% set(gca,'color',[0.9,0.9,0.9]);
% 
% %relative dynamics of position
% figure('Name', 'Relative position', 'NumberTitle', 'off')
% subplot(2,1,1)
% plot(1:Tot+1,z(1,:),1:Tot+1,z(3,:),1:Tot+1,z(5,:),1:Tot+1,z(7,:),1:Tot+1,z(9,:),1:Tot+1,z(11,:),'linewidth',1.8)
% yline(0,'--k')
% legend('$p_x^1 - p_x^2$','$p_x^1 - p_x^3$','$p_x^2 - p_x^3$','$p_x^4 - p_x^3$','$p_x^1 - p_x^4$','$p_x^2 - p_x^4$','fontsize',12,'interpreter','latex')
% xlim([1 Tot+1])
% xlabel('time-steps','fontsize',12)
% ylabel('Relative position','fontsize',12)
% % title('Relative dynamics of position in x direction')
% grid on
% set(gca,'color',[0.9,0.9,0.9]);
% subplot(2,1,2)
% plot(1:Tot+1,z(2,:),1:Tot+1,z(4,:),1:Tot+1,z(6,:),1:Tot+1,z(8,:),1:Tot+1,z(10,:),1:Tot+1,z(12,:),'linewidth',1.8)
% yline(0,'--k')
% legend('$p_y^1 - p_y^2$','$p_y^1 - p_y^3$','$p_y^2 - p_y^3$','$p_y^4 - p_y^3$','$p_y^1 - p_y^4$','$p_y^2 - p_y^4$','fontsize',12,'interpreter','latex')
% xlim([1 Tot+1])
% xlabel('time-steps','fontsize',12)
% ylabel('Relative position','fontsize',12)
% % title('Relative dynamics of position in y direction')
% grid on
% set(gca,'color',[0.9,0.9,0.9]);
% 
% % %relative dynamics of velocity
% figure('Name', 'Relative velocity', 'NumberTitle', 'off')
% subplot(2,1,1)
% plot(1:Tot+1,z(13,:),1:Tot+1,z(15,:),1:Tot+1,z(17,:),1:Tot+1,z(19,:),1:Tot+1,z(21,:),1:Tot+1,z(23,:),'linewidth',1.8)
% yline(0,'--k')
% legend('$v_x^1 - v_x^2$','$v_x^1 - v_x^3$','$v_x^2 - v_x^3$','$v_x^4 - v_x^3$','$v_x^1 - v_x^4$','$v_x^2 - v_x^4$','fontsize',12,'interpreter','latex')
% xlim([1 Tot+1])
% xlabel('time-steps','fontsize',12)
% ylabel('Relative velocity','fontsize',12)
% % title('Relative dynamics of velocity in x direction')
% grid on
% set(gca,'color',[0.9,0.9,0.9]);
% subplot(2,1,2)
% plot(1:Tot+1,z(14,:),1:Tot+1,z(16,:),1:Tot+1,z(18,:),1:Tot+1,z(20,:),1:Tot+1,z(22,:),1:Tot+1,z(24,:),'linewidth',1.8)
% yline(0,'--k')
% legend('$v_y^1 - v_y^2$','$v_y^1 - v_y^3$','$v_y^2 - v_y^3$','$v_y^4 - v_y^3$','$v_y^1 - v_y^4$','$v_y^2 - v_y^4$','fontsize',12,'interpreter','latex')
% xlim([1 Tot+1])
% xlabel('time-steps','fontsize',12)
% ylabel('Relative velocity','fontsize',12)
% % title('Relative dynamics of velocity in y direction')
% grid on
% set(gca,'color',[0.9,0.9,0.9]);

%% Comparison of Nash and optimal solution from decoupling framework
%to plot control signals - Figure 3
figure('Name', 'Comparison of optimal control signals', 'NumberTitle', 'off')
subplot(2,4,1)
plot(1:Tot,un(1,:),'-b',1:Tot,uhat(1,:),'--r','linewidth',1.5)
legend('Nash','Decoupling','fontsize',11,'location','best')
title('Agent 1 (x-axis)','fontweight','bold','fontsize',12)
xlim([1 Tot])
grid on
xlabel('time steps','fontsize',12)
ylabel('control inputs','fontsize',12)
set(gca,'color',[0.9,0.9,0.9]);
subplot(2,4,2)
plot(1:Tot,un(3,:),'-b',1:Tot,uhat(3,:),'--r','linewidth',1.5)
legend('Nash','Decoupling','fontsize',11,'location','best')
title('Agent 2 (x-axis)','fontweight','bold','fontsize',12)
xlim([1 Tot])
grid on
xlabel('time steps','fontsize',12)
set(gca,'color',[0.9,0.9,0.9]);
subplot(2,4,3)
plot(1:Tot,un(5,:),'-b',1:Tot,uhat(5,:),'--r','linewidth',1.5)
legend('Nash','Decoupling','fontsize',11,'location','best')
title('Agent 3 (x-axis)','fontweight','bold','fontsize',12)
xlim([1 Tot])
grid on
xlabel('time steps','fontsize',12)
set(gca,'color',[0.9,0.9,0.9]);
subplot(2,4,4)
plot(1:Tot,un(7,:),'-b',1:Tot,uhat(7,:),'--r','linewidth',1.5)
legend('Nash','Decoupling','fontsize',11,'location','best')
title('Agent 4 (x-axis)','fontweight','bold','fontsize',12)
xlim([1 Tot])
grid on
xlabel('time steps','fontsize',12)
set(gca,'color',[0.9,0.9,0.9]);

subplot(2,4,5)
plot(1:Tot,un(2,:),'-b',1:Tot,uhat(2,:),'--r','linewidth',1.5)
legend('Nash','Decoupling','fontsize',11,'location','best')
title('Agent 1 (y-axis)','fontweight','bold','fontsize',12)
xlim([1 (tf/dt)-1])
grid on
ylabel('control inputs','fontsize',12)
xlabel('time steps','fontsize',12)
set(gca,'color',[0.9,0.9,0.9]);
subplot(2,4,6)
plot(1:Tot,un(4,:),'-b',1:Tot,uhat(4,:),'--r','linewidth',1.5)
legend('Nash','Decoupling','fontsize',11,'location','best')
title('Agent 2 (y-axis)','fontweight','bold','fontsize',12)
xlim([1 (tf/dt)-1])
grid on
xlabel('time steps','fontsize',12)
set(gca,'color',[0.9,0.9,0.9]);
subplot(2,4,7)
plot(1:Tot,un(6,:),'-b',1:Tot,uhat(6,:),'--r','linewidth',1.5)
legend('Nash','Decoupling','fontsize',11,'location','best')
title('Agent 3 (y-axis)','fontweight','bold','fontsize',12)
xlim([1 (tf/dt)-1])
grid on
xlabel('time steps','fontsize',12)
set(gca,'color',[0.9,0.9,0.9]);
subplot(2,4,8)
plot(1:Tot,un(8,:),'-b',1:Tot,uhat(8,:),'--r','linewidth',1.5)
legend('Nash','Decoupling','fontsize',11,'location','best')
title('Agent 4 (y-axis)','fontweight','bold','fontsize',12)
xlim([1 (tf/dt)-1])
grid on
xlabel('time steps','fontsize',12)
set(gca,'color',[0.9,0.9,0.9]);


%% Consensus section
%to plot the position alignment - Figure 4

figure('Name', 'Consensus', 'NumberTitle', 'off')
plot1 = plot(xn(1,:),xn(2,:),'-b','linewidth',3.5);
hold on
plot2 = plot(xn(3,:),xn(4,:),'-b','linewidth',3.5);
hold on
plot3 = plot(xn(5,:),xn(6,:),'-b','linewidth',3.5);
hold on
plot4 = plot(xn(7,:),xn(8,:),'-b','linewidth',3.5);
hold on
plot5 = plot(xhat(1,:),xhat(2,:),':r','linewidth',3.5);
hold on
plot6 = plot(xhat(3,:),xhat(4,:),':r','linewidth',3.5);
hold on
plot7 = plot(xhat(5,:),xhat(6,:),':r','linewidth',3.5);
hold on
plot8 = plot(xhat(7,:),xhat(8,:),':r','linewidth',3.5);
set(get(get(plot2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
set(get(get(plot3, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
set(get(get(plot4, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
set(get(get(plot6, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
set(get(get(plot7, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
set(get(get(plot8, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
ag1 = {'Agent 1'}; ag2 = {'Agent 2'}; ag3 = {'Agent 3'}; ag4 = {'Agent 4'};
text(0.2,5,ag1,'Fontsize',12)
text(0.5,0.75,ag2,'Fontsize',12)
text(4.3,5.2,ag3,'Fontsize',12)
text(4.2,0.75,ag4,'Fontsize',12)
plot(xn(1,1),xn(2,1),'ob',xn(3,1),xn(4,1),'ob',xn(5,1),xn(6,1),'ob',xn(7,1),xn(8,1),'ob','linewidth',5)
plot(xn(1,end),xn(2,end),'^k',xn(3,end),xn(4,end),'^k',xn(5,end),xn(6,end),'^k',xn(7,end),xn(8,end),'^k','linewidth',5)
legend('Nash strategy','Decoupling approach','fontsize',12);
xlabel('x-axis','fontsize',12)
ylabel('y-axis','fontsize',12)
grid on
set(gca,'color',[0.9,0.9,0.9]);


