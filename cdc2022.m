clc;clear all;close all;
%% code for paper: A Distributed LQDTG to MAS Consensus
%  prima.aditya@tuhh.de

n    = 2;          %dimensional plane
I    = eye(n); 
N    = 4;          %amount of agents
r    = n*N;        %adjust the size
M    = 4;          %number of edges
h    = n*M;        %adjust the size
tf   = 10;         %try 5 for receding horizon
Np   = 0;          %try 5 prediction horizon
step = tf*10;
dt   = tf/step; %step size/ sampling time
eta  = 10;
%incidence matrix
D = [-1 0 0 -1;
    1 -1 0 0;
    0 1 -1 1;
    0 0 1 0];
Di = kron(D,eye(n));
miu = 1;
%% Problem 1 - Nash strategy
A  = [zeros(r), eye(r);zeros(r), zeros(r)];
F  = eye(2*r) + dt*A;
mu = abs(D);
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
p1  = [0.5; 5];  p3  = [5; 4];
p2  = [0; 0.5];  p4  = [4; 0];    
v1  = [0; 1];    v3  = [0; 1];
v2  = [0; 1];    v4  = [0; 1];    
po  = [p1;p2;p3;p4];
vo  = [v1;v2;v3;v4];

%% Arranging Problem 1 as a single (coupled) Riccati problem - 1a
Wn   = miu*eye(M);
Ln   = D*Wn*D';
Qn   = dt*[kron(Ln,I), zeros(r); zeros(r), kron(Ln,I)];
Rn   = dt*eye(r);
Qfn  = eta*Qn;
%matrix coefficients for vs
Bn   = [zeros(r);eye(r)];
%exact discretization
Gn   = dt*Bn + (dt^2)/2*A*Bn;

%% Arranging Problem 2 as a single Riccati problem
%matrix coefficients for es
At    = [zeros(h), eye(h);zeros(h), zeros(h)];
Ft    = eye(2*h) + dt*At; %Euler method
mut   = eye(M); %to define positive weights of each agent -> W
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

alpha   = 0.1;

Iter    = 50; 
% Iter    = 200; %for receding horizon
In      = eye(N*n);
u0      = zeros(r,1);
eps     = 1e-3;
udata   = [];
for idx=1:Np+1%indexing for prediction horizon
    if idx==1
        z(:,idx)= [po(1:2)-po(3:4); po(3:4)-po(5:6); po(5:6)-po(7:8); po(1:2)-po(5:6);...
                   vo(1:2)-vo(3:4); vo(3:4)-vo(5:6); vo(5:6)-vo(7:8); vo(1:2)-vo(5:6);];
        %have to initialize the position and velocity vector
        x(:,idx)     = [po;vo];
        xpinv(:,idx) = [po;vo];
        xhat(:,idx)  = [po;vo];
    else
        idxn = ((idx-1) * 10) + 1;
        z(:,idxn) = [xf{idx-1}(1:2,idxn)- xf{idx-1}(3:4,idxn); xf{idx-1}(3:4,idxn)- xf{idx-1}(5:6,idxn);xf{idx-1}(5:6,idxn)- xf{idx-1}(7:8,idxn);xf{idx-1}(1:2,idxn)- xf{idx-1}(5:6,idxn);xf{idx-1}(9:10,idxn)- xf{idx-1}(11:12,idxn); xf{idx-1}(11:12,idxn)- xf{idx-1}(13:14,idxn); xf{idx-1}(13:14,idxn)- xf{idx-1}(15:16,idxn); xf{idx-1}(9:10,idxn)- xf{idx-1}(13:14,idxn)];
        %if reced the horizon then x(l) be the initial state
        x(:,idxn)    = xd{idx-1}(:,idxn);
        xpinv(:,idxn)= xe{idx-1}(:,idxn);
        xhat(:,idxn) = xf{idx-1}(:,idxn);
    end
    %define the moving finite time horizon
    T = (idx-1)*10 + step;
    %initialize the Riccati solution to be solved backward 
    for l=1:N %for amount of robots 
        P(:,:,l,T+1) = Qf(:,:,l); %initial Riccati for Problem 1
    end
    Pn(:,:,T+1) = Qfn; %initial Riccati for Problem 1a
    [Qft, Kft, Eige] = idare(Ft,Gt,Qt,Rt);
    Pt(:,:,T+1) = Qft; %initial Riccati for Problem 2
    
    
    for k = 1:T
            %-------------------------------------------------------
        %solve Problem 1
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
        uNash(:,k) = [u(:,:,1,k);u(:,:,2,k);u(:,:,3,k);u(:,:,4,k)];
        x(:,k+1)    = F * x(:,k) +  G(:,:,1)*u(:,:,1,k) + G(:,:,2)*u(:,:,2,k) + G(:,:,3)*u(:,:,3,k) + G(:,:,4)*u(:,:,4,k);%    x(:,k+1) = inv(La(:,:,k)) * F * x(:,k);

%         %solve Problem 1a (as a single Riccati - not a game)
%         for j = T-1:-1:1%backward computation
%             Pn(:,:,j) = Qn + F'*Pn(:,:,j+1)*F - F'*Pn(:,:,j+1)*Gn*inv(Rn+Gn'*Pn(:,:,j+1)*Gn)*Gn'*Pn(:,:,j+1)*F;
%         end
%         Kn(:,:,k) = -inv(Rn+Gn'*Pn(:,:,k+1)*Gn)*Gn'*Pn(:,:,k+1)*F;
%         un(:,k)   = Kn(:,:,k) * xn(:,k);
%         xn(:,k+1) = F * xn(:,k) + Gn*un(:,k);

        %solve Problem 2 (distributed manner)
        for j = T-1:-1:1%backward computation
        Pt(:,:,j) = Qt + Ft'*Pt(:,:,j+1)*Ft - Ft'*Pt(:,:,j+1)*Gt*inv(Rt+Gt'*Pt(:,:,j+1)*Gt)*Gt'*Pt(:,:,j+1)*Ft;
        end
        Kt(:,:,k) = -inv(Rt+Gt'*Pt(:,:,k+1)*Gt)*Gt'*Pt(:,:,k+1)*Ft;
        a(:,k)    = Kt(:,:,k) * z(:,k);
        z(:,k+1)  = Ft * z(:,k) + Gt * a(:,k);
            %test comparison between uhat != upinv
            for v = 1:N
                 upinv(:,k) = Phiai * a(:,k);
            end
            xpinv(:,k+1) = F * xpinv(:,k) +G(:,1:2)*upinv(1:2,k) + G(:,3:4)*upinv(3:4,k) + G(:,5:6)*upinv(5:6,k) + G(:,7:8)*upinv(7:8,k) ;
       
            %iterative method to get back uhat from stored a
            %tic;
            t = 1;
            if k>1
                u0 = usol;
            end
%             while (norm(a(:,k)-Phia*u0(:,t)) > eps) && (t <= Iter)%for receding horizon
           while (t <= Iter)
                u0_temp = (In - 2*alpha*Phia'*Phia)*u0(:,t) + 2*alpha*Phia'*a(:,k); 
                u0 = [u0 u0_temp];
            t = t+1;
            end
            udata = [udata,u0];
            %toc;
            usol      = u0(:,end);
            uhat(:,k) = usol;   
        xhat(:,k+1) = F * xhat(:,k) +G(:,1:2)*uhat(1:2,k) + G(:,3:4)*uhat(3:4,k) + G(:,5:6)*uhat(5:6,k) + G(:,7:8)*uhat(7:8,k) ;
        %calculate the receding stability equation 
        term1 = (Ft + Gt*Kt(:,:,end))' * Qft * (Ft + Gt*Kt(:,:,end)) - Qft;
        term2 = -Qt - Kt(:,:,end)'*Rt*Kt(:,:,end);
    end
    xd{idx} = x;
    xe{idx} = xpinv;
    xf{idx} = xhat;
    zf{idx} = z; 
end

Tot = (tf*10)+(Np*10);
%% Regulator section - Figure 2
figure('Name', 'Relative in x', 'NumberTitle', 'off')
subplot(3,1,1)
plot(1:Tot+1,z(1,:),1:Tot+1,z(3,:),1:Tot+1,z(5,:),1:Tot+1,z(7,:),'linewidth',1.8)
yline(0,'--k');
legend('$q_x^1$','$q_x^2$','$q_x^3$','$q_x^4$','fontsize',12,'interpreter','latex')
xlim([1 Tot+1])
xlabel('time-steps','fontsize',12)
ylabel('x-axis','fontsize',12)
title('Relative positions')
grid on
subplot(3,1,2)
plot(1:Tot+1,z(9,:),1:Tot+1,z(11,:),1:Tot+1,z(13,:),1:Tot+1,z(15,:),'linewidth',1.8)
yline(0,'--k');
legend('$w_x^1$','$w_x^2$','$w_x^3$','$w_x^4$','fontsize',12,'interpreter','latex')
xlim([1 Tot+1])
xlabel('time-steps','fontsize',12)
ylabel('x-axis','fontsize',12)
title('Relative velocities')
grid on
subplot(3,1,3)
plot(1:Tot,a(1,:),1:Tot,a(3,:),1:Tot,a(5,:),1:Tot,a(7,:),'linewidth',1.8)
yline(0,'--k');
legend('$a_x^1$','$a_x^2$','$a_x^3$','$a_x^4$','fontsize',12,'interpreter','latex')
xlim([1 Tot])
xlabel('time-steps','fontsize',12)
ylabel('x-axis','fontsize',12)
title('Relative control inputs')
grid on
xlim([1 Tot])

figure('Name', 'Relative in y', 'NumberTitle', 'off')
subplot(3,1,1)
plot(1:Tot+1,z(2,:),1:Tot+1,z(4,:),1:Tot+1,z(6,:),1:Tot+1,z(8,:),'linewidth',1.8)
yline(0,'--k');
legend('$q_y^1$','$q_y^2$','$q_y^3$','$q_y^4$','fontsize',12,'interpreter','latex')
xlim([1 Tot+1])
xlabel('time-steps','fontsize',12)
ylabel('y-axis','fontsize',12)
title('Relative positions')
grid on
subplot(3,1,2)
plot(1:Tot+1,z(10,:),1:Tot+1,z(12,:),1:Tot+1,z(14,:),1:Tot+1,z(16,:),'linewidth',1.8)
yline(0,'--k');
legend('$w_y^1$','$w_y^2$','$w_y^3$','$w_y^4$','fontsize',12,'interpreter','latex')
xlim([1 Tot+1])
xlabel('time-steps','fontsize',12)
ylabel('y-axis','fontsize',12)
title('Relative velocities')
grid on
subplot(3,1,3)
plot(1:Tot,a(2,:),1:Tot,a(4,:),1:Tot,a(6,:),1:Tot,a(8,:),'linewidth',1.8)
yline(0,'--k');
legend('$a_y^1$','$a_y^2$','$a_y^3$','$a_y^4$','fontsize',12,'interpreter','latex')
xlim([1 Tot])
xlabel('time-steps','fontsize',12)
ylabel('y-axis','fontsize',12)
title('Relative control inputs')
grid on
xlim([1 Tot])
%%
% % Agent 1 for iteration Iter = 2;
% fig = figure(99);clf;
% ax  = axes;
% plot(ax,udata(1,:),'*b')
% hold on
% for i = 1:length(upinv(1,:))
%     plot(ax,(Iter+1)*i,upinv(1,i),'xr','linewidth',5);
% end
% grid on
% xlabel('Time/ iterations')
% ylabel('Control input of Agent 1')
% legend('Distributed approach','Centralized solution','fontsize',12,'location','northeast');
% xlim([0 200])
% % create a new axes object (you could also get the active axes object via 'gca')
% % define Name-Value pairs for the zoom_plot function:
% % Name-Value pairs for the axes:
% options.axes.Names = {'Position','XLim'};
% options.axes.Values = {[.5 .4 .35 .35],[20,40]};
% % Name-Value pairs for the rectangle:
% options.rectangle.Names = {};
% options.rectangle.Values = {};
% % Name-Value pairs for the arrows:
% options.arrows.Names = {'HeadLength','HeadWidth'};
% options.arrows.Values = {8,8};
% % call the function with options:
% [zoom_utils] = zoom_plot(ax,options);
% 
% % Agent 4 for iteration Iter = 2;
% fig = figure(98);clf;
% ax  = axes;
% plot(ax,udata(7,:),'*b')
% hold on
% for i = 1:length(upinv(7,:))
%     plot(ax,(Iter+1)*i,upinv(7,i),'xr','linewidth',5);
% end
% grid on
% xlabel('Time/ iterations')
% ylabel('Control input of Agent 4')
% legend('Distributed approach','Centralized solution','fontsize',12,'location','southeast');
% xlim([0 200])
% % create a new axes object (you could also get the active axes object via 'gca')
% % define Name-Value pairs for the zoom_plot function:
% % Name-Value pairs for the axes:
% options.axes.Names = {'Position','XLim'};
% options.axes.Values = {[.5 .4 .3 .3],[20,40]};
% % Name-Value pairs for the rectangle:
% options.rectangle.Names = {};
% options.rectangle.Values = {};
% % Name-Value pairs for the arrows:
% options.arrows.Names = {'HeadLength','HeadWidth'};
% options.arrows.Values = {8,8};
% % call the function with options:
% [zoom_utils] = zoom_plot(ax,options);
%%
% Agent 1 for iteration Iter = 50;
fig = figure(97);clf;
ax  = axes;
plot(ax,udata(1,:),'*b')
hold on
for i = 1:length(upinv(1,:))
    plot(ax,(Iter+1)*i,upinv(1,i),'xr','linewidth',5);
end
grid on
xlabel('Time/ iterations')
ylabel('Control input of Agent 1')
legend('Distributed approach','Centralized solution','fontsize',12,'location','southwest');
xlim([0 1000])
% create a new axes object (you could also get the active axes object via 'gca')
% define Name-Value pairs for the zoom_plot function:
% Name-Value pairs for the axes:
options.axes.Names = {'Position','XLim'};
options.axes.Values = {[.5 .5 .35 .35],[100,220]};
% Name-Value pairs for the rectangle:
options.rectangle.Names = {};
options.rectangle.Values = {};
% Name-Value pairs for the arrows:
options.arrows.Names = {'HeadLength','HeadWidth'};
options.arrows.Values = {8,8};
% call the function with options:
[zoom_utils] = zoom_plot(ax,options);

% Agent 4 for iteration Iter = 50;
fig = figure(96);clf;
ax  = axes;
plot(ax,udata(7,:),'*b')
hold on
for i = 1:length(upinv(7,:))
    plot(ax,(Iter+1)*i,upinv(7,i),'xr','linewidth',5);
end
grid on
xlabel('Time/ iterations')
ylabel('Control input of Agent 4')
legend('Distributed approach','Centralized solution','fontsize',12,'location','northwest');
xlim([0 1000])
% create a new axes object (you could also get the active axes object via 'gca')
% define Name-Value pairs for the zoom_plot function:
% Name-Value pairs for the axes:
options.axes.Names = {'Position','XLim'};
options.axes.Values = {[.5 .3 .35 .35],[100,220]};
% Name-Value pairs for the rectangle:
options.rectangle.Names = {};
options.rectangle.Values = {};
% Name-Value pairs for the arrows:
options.arrows.Names = {'HeadLength','HeadWidth'};
options.arrows.Values = {9,9};
% call the function with options:
[zoom_utils] = zoom_plot(ax,options);


%% Consensus section
%to plot the position alignment - Figure 4

figure('Name', 'Consensus', 'NumberTitle', 'off')

plot1 = plot(xhat(1,:),xhat(2,:),'-b','linewidth',3.5);
hold on
plot2 = plot(xhat(3,:),xhat(4,:),'-b','linewidth',3.5);
hold on
plot3 = plot(xhat(5,:),xhat(6,:),'-b','linewidth',3.5);
hold on
plot4 = plot(xhat(7,:),xhat(8,:),'-b','linewidth',3.5);
hold on
plot5 = plot(xpinv(1,:),xpinv(2,:),':r','linewidth',3.5);
hold on
plot6 = plot(xpinv(3,:),xpinv(4,:),':r','linewidth',3.5);
hold on
plot7 = plot(xpinv(5,:),xpinv(6,:),':r','linewidth',3.5);
hold on
plot8 = plot(xpinv(7,:),xpinv(8,:),':r','linewidth',3.5);
hold on
set(get(get(plot2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
set(get(get(plot3, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
set(get(get(plot4, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
set(get(get(plot6, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
ag1 = {'Agent 1'}; ag2 = {'Agent 2'}; ag3 = {'Agent 3'}; ag4 = {'Agent 4'};
text(0.1,6,ag1,'Fontsize',12)
text(0.5,0.75,ag2,'Fontsize',12)
text(4.3,5.2,ag3,'Fontsize',12)
text(4.2,0.75,ag4,'Fontsize',12)
plot(x(1,1),x(2,1),'ob',x(3,1),x(4,1),'ob',x(5,1),x(6,1),'ob',x(7,1),x(8,1),'ob','linewidth',5)
plot(x(1,end),x(2,end),'^k',x(3,end),x(4,end),'^k',x(5,end),x(6,end),'^k',x(7,end),x(8,end),'^k','linewidth',5)
legend('Distributed approach','Centralized solution','fontsize',12);
xlabel('x-axis','fontsize',12)
ylabel('y-axis','fontsize',12)
grid on
