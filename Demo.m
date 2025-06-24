clear all;
close all;
clc
%%
N=3000;T=0.01;

q=zeros(3,N);dq=zeros(3,N);ddq=zeros(3,N);t=zeros(1,N);tao=zeros(3,N);fmax=zeros(3,N);Au=zeros(3,3,N);theta=zeros(3,N);
q2=zeros(3,N);dq2=zeros(3,N);ddq2=zeros(3,N);tao2=zeros(3,N);fmax2=zeros(3,N);Au2=zeros(3,3,N);theta2=zeros(3,N);
q3=zeros(3,N);dq3=zeros(3,N);ddq3 =zeros(3,N);tao3=zeros(3,N);
bu2=zeros(3,N);
x1=zeros(3,N);x2=zeros(3,N);x3=zeros(3,N);dx1=zeros(3,N);dx2=zeros(3,N);dx3=zeros(3,N); 
Ax1=zeros(3,N);Ax2=zeros(3,N);Ax3=zeros(3,N);dAx1=zeros(3,N);dAx2=zeros(3,N);dAx3=zeros(3,N); 
e0=zeros(3,N);de0=zeros(3,N);
e=zeros(3,N);de=zeros(3,N);
e2=zeros(3,N);de2=zeros(3,N);
e3=zeros(3,N);de3=zeros(3,N);

dq(:,1)=[0; 0 ;0];dq2(:,1)=[0; 0 ;0];dq3(:,1)=[0; 0 ;0];
 q(:,1)=[0.8;-0.6;0.001];q2(:,1)=[0.8;-0.6;0.001];q3(:,1)=[0.8;-0.6;0.001];
% q(:,1)=[0.8;0.3;0.01];q2(:,1)=[0.8;0.3;0.01];q3(:,1)=[0.8;0.3;0.01];
x1(:,1)=[0.8;-0.6;0.001]; Ax1(:,1)=[0.8;-0.6;0.001]; 
d=zeros(3,N);
D=eye(3,3)/10;K=eye(3,3);
l=5;l2=8;
% w(:,1) =5*rand(6,1);
w(:,1) =zeros(6,1);w2(:,1) =zeros(6,1);
P = zeros(6,6);P2 = zeros(6,6);
Q = zeros(6,1);Q2 = zeros(6,1);
Tau=0.001;Tau2=0.003;

%SMC
kp1=1;kd1=1;ki1=1;
lenta=10;yita2=0.5;
%SMC-ESO
a1=10;a2=80; a3=15;asplon=0.1;
kp2=100;kd2=20;c=10;

  %%Reference
for i=1:3000
    t(i)=i*0.01;
     if i<500
     qd(1,i)=0.4;
     elseif 500<=i && i<700
     qd(1,i)=-0.001*i+0.9;
     
     elseif 700<=i && i<1000
     qd(1,i)=0.2;
   
     elseif 1000<=i && i<1500
     qd(1,i)=0.001*i-0.8;
     
     elseif 1500<=i && i<2000
     qd(1,i)=-0.002*i+3.7;
     
     elseif 2000<=i && i<2500
     qd(1,i)=-0.3;
    
     elseif 2500<=i && i<2700
     qd(1,i)=0.001*i-2.8;
   
     elseif 2700<=i && i<3001
     qd(1,i)=-0.1;
     
     end
     qd(2,i)=-qd(1,i)+0.1;

      qd(3,i)=0.1*cos(t(i));
end
%%


for i=1:N
  t(i)=i*T;
   



%   qd(:,i)=[cos(t(i))+sin(t(i)/2);sin(t(i))+0.5*cos(t(i)/2);0.2*cos(t(i))];
%   dqd(:,i)=[-sin(t(i))+0.5*cos(t(i)/2);cos(t(i))-0.25*sin(t(i)/2);-0.2*sin(t(i))];    
%   ddqd(:,i)=[-cos(t(i))-0.25*sin(t(i)/2);-sin(t(i))-0.125*cos(t(i)/2);-0.2*cos(t(i))];
% QD=load('qd.mat');
% qd=QD.qd;
dqd=zeros(3,N);
ddqd=zeros(3,N);
dqd(3,i)=-0.2*sin(t(i));ddqd(3,i)=-0.2*cos(t(i));

     d(:,i)=2*[(sin(t(i))+sin(0.5*t(i)))/2;(cos(t(i))+cos(0.5*t(i)))/2;(sin(t(i))+cos(0.5*t(i)))/4];
  %Model
  B(:,:,i)=qinputB(q(:,i));
  C(:,:,i)=qinputC(dq(:,i),q(:,i));

  B2(:,:,i)=qinputB(q2(:,i));
  C2(:,:,i)=qinputC(dq2(:,i),q2(:,i));

  B3(:,:,i)=qinputB(q3(:,i));
  C3(:,:,i)=qinputC(dq3(:,i),q3(:,i));

 %Errors
  e(:,i)=-qd(:,i)+q(:,i);
  de(:,i)=-dqd(:,i)+dq(:,i);%

  e2(:,i)=-qd(:,i)+q2(:,i);
  de2(:,i)=-dqd(:,i)+dq2(:,i);

  e3(:,i)=-qd(:,i)+q3(:,i);
  de3(:,i)=-dqd(:,i)+dq3(:,i); 

  eg2(:,i)=-qd(:,i)+x1(:,i);deg2(:,i)=-dqd(:,i)+x2(:,i);
  sg2(:,i)=c*eg2(:,i)+deg2(:,i);

 %ESO
    Ae0(:,i)=q(:,i)-Ax1(:,i);
    dAx1(:,i)=Ax2(:,i)+a1*Ae0(:,i)/asplon;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
    Ax1(:,i+1)=Ax1(:,i)+T*dAx1(:,i);
    dAx2(:,i)=tao(:,i)+Ax3(:,i)+a2*Ae0(:,i)/asplon^2;
    Ax2(:,i+1)=Ax2(:,i)+T*dAx2(:,i);
    dAx3(:,i)=a3*Ae0(:,i)/asplon^3;
    Ax3(:,i+1)=Ax3(:,i)+T*dAx3(:,i);

    e0(:,i)=q3(:,i)-x1(:,i);
    dx1(:,i)=x2(:,i)+a1*e0(:,i)/asplon;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
    x1(:,i+1)=x1(:,i)+T*dx1(:,i);
    dx2(:,i)=bu2(:,i)+x3(:,i)+a2*e0(:,i)/asplon^2;
    x2(:,i+1)=x2(:,i)+T*dx2(:,i);
    dx3(:,i)=a3*e0(:,i)/asplon^3;
    x3(:,i+1)=x3(:,i)+T*dx3(:,i);



rou(:,:,i)=[e(1,i)^2 e(1,i)*e(2,i) e(1,i)*e(3,i) e(2,i)^2  e(2,i)*e(3,i)  e(3,i)^2];
drou(:,:,i)=[2*e(1,i) e(2,i) e(3,i) 0 0 0;0 e(1,i) 0 2*e(2,i) e(3,i) 0; 0 0 e(1,i) 0 e(2,i) 2*e(3,i) ]; 

rou2(:,:,i)=[e2(1,i)^2 e2(1,i)*e2(2,i) e2(1,i)*e2(3,i) e2(2,i)^2  e2(2,i)*e2(3,i)  e2(3,i)^2];
drou2(:,:,i)=[2*e2(1,i) e2(2,i) e2(3,i) 0 0 0;0 e2(1,i) 0 2*e2(2,i) e2(3,i) 0; 0 0 e2(1,i) 0 e2(2,i) 2*e2(3,i) ]; 

tao(:,i) =-1/2 * B(:,:,i)* drou(:,:,i)* w(:,i)+1*(B(:,:,i)*ddqd(:,i)-25*e(:,i)-5*de(:,i))-Ax3(:,i);
tao2(:,i) =-1/2 * B2(:,:,i)* drou2(:,:,i)* w2(:,i)+1*(B2(:,:,i)*ddqd(:,i)-25*e2(:,i)-15*de2(:,i));

bu2(:,i)=ddqd(:,i)-kd2*deg2(:,i)-kp2*eg2(:,i)-x3(:,i)-yita2*sign(sg2(:,i));
tao3(:,i)=B3(:,:,i)*bu2(:,i);


 fmax(:,i)=[norm([e(1,i) de(1,i)]);norm([e(2,i) de(2,i)]);norm([e(3,i) de(3,i)])];
fmax2(:,i)=[norm([e2(1,i) de2(1,i)]);norm([e2(2,i) de2(2,i)]);norm([e2(3,i) de2(3,i)])];


Au(:,1,i)=-B(:,:,i)\((C(:,:,i)+D)*dq(:,i)+K*q(:,i)-tao(:,i));
Au(:,2,i)=dq(:,i);
Au(:,3,i)=q(:,i);
E(:,:,i) = drou(:,:,i)'*Au(:,:,i);

Au2(:,1,i)=-B2(:,:,i)\((C2(:,:,i)+D)*dq2(:,i)+K*q2(:,i)-tao2(:,i));
Au2(:,2,i)=dq2(:,i);
Au2(:,3,i)=q2(:,i);
E2(:,:,i) = drou2(:,:,i)'*Au2(:,:,i);


theta (:,i)= [ e(1,i)^2 + de(1,i)^2 + tao(1,i)^2 ; + e(2,i)^2 + de(2,i)^2 + tao(2,i)^2;e(3,i)^2 + de(3,i)^2 + tao(3,i)^2 ];
theta2(:,i)= [  e2(1,i)^2  + tao2(1,i)^2 + fmax2(1,i)^2 ;  e2(2,i)^2  + tao2(2,i)^2 + fmax2(2,i)^2;  e2(3,i)^2  + tao2(3,i)^2 + fmax2(3,i)^2];

 
dP(:,:,i) = -l*P(:,:,i)+E(:,:,i)*E(:,:,i)';
dQ (:,i)= -l*Q(:,i)+E(:,:,i)*theta(:,i);
P(:,:,i+1) = P(:,:,i)+dP(:,:,i)*T;
Q(:,:,i+1) = Q(:,i)+dQ(:,i)*T;
W(:,i) = P(:,:,i)*w(:,i)+Q(:,i); 
w(:,i+1) = w(:,i) - Tau*W(:,i);  

dP2(:,:,i) = -l2*P2(:,:,i)+E2(:,:,i)*E2(:,:,i)';
dQ2(:,i)= -l2*Q2(:,i)+E2(:,:,i)*theta2(:,i);
P2(:,:,i+1) = P2(:,:,i)+dP2(:,:,i)*T;
Q2(:,:,i+1) = Q2(:,i)+dQ2(:,i)*T;
W2(:,i) = P2(:,:,i)*w2(:,i)+Q2(:,i);
w2(:,i+1) = w2(:,i) - Tau2*W2(:,i);  





%Model updation
ddq(:,i)=-B(:,:,i)\((C(:,:,i)+D)*dq(:,i)+K*q(:,i)-tao(:,i))+d(:,i);
dq(:,i+1)=dq(:,i)+ddq(:,i)*T;
q(:,i+1)=q(:,i)+dq(:,i)*T;

ddq2(:,i)=-B2(:,:,i)\((C2(:,:,i)+D)*dq2(:,i)+K*q2(:,i)-tao2(:,i)+d(:,i));
dq2(:,i+1)=dq2(:,i)+ddq2(:,i)*T;
q2(:,i+1)=q2(:,i)+dq2(:,i)*T;

ddq3(:,i)=-B3(:,:,i)\((C3(:,:,i)+D)*dq3(:,i)+K*q3(:,i)-tao3(:,i)+d(:,i));
dq3(:,i+1)=dq3(:,i)+ddq3(:,i)*T;
q3(:,i+1)=q3(:,i)+dq3(:,i)*T;



posd(:,i)=Positionq(qd(:,i));
pos1(:,i)=Positionq(q(:,i));
pos2(:,i)=Positionq(q2(:,i));
pos3(:,i)=Positionq(q3(:,i));
V(i)=rou(:,:,i)*w(:,i);
k(i)=kappa(q(:,i));
p(i)=phi(q(:,i));
k2(i)=kappa(qd(:,i));
p2(i)=phi(qd(:,i));
end

q(:,N+1)=[];dq(:,N+1)=[];w(:,N+1)=[];
q2(:,N+1)=[];dq2(:,N+1)=[];w2(:,N+1)=[];
q3(:,N+1)=[];dq3(:,N+1)=[];

%Index


EE3=mse(e3(1,:)+e3(2,:)+e3(3,:));
EE2=mse(e2(1,:)+e2(2,:)+e2(3,:));
EE=mse(e(1,:)+e(2,:)+e(3,:));

% EEE3=abs(e3);EEE2=abs(e2);EEE=abs(e);
% EE3=mean(EEE3(1,:)+EEE3(2,:)+EEE3(3,:));
% EE2=mean(EEE2(1,:)+EEE2(2,:)+EEE2(3,:));
% EE=mean(EEE(1,:)+EEE(2,:)+EEE(3,:));

Tao3=abs(tao3);Tao2=abs(tao2);Tao=abs(tao);
UU3=mean(Tao3(1,:)+Tao3(2,:)+Tao3(3,:));
UU2=mean(Tao2(1,:)+Tao2(2,:)+Tao2(3,:));
UU=mean(Tao(1,:)+Tao(2,:)+Tao(3,:));

%Figures
figure(1)
plot(t,qd(1,:),'k--','LineWidth',1.5)
hold on 
plot(t,q3(1,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on   
plot(t,q2(1,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
hold on 
plot(t,q(1,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
legend('Reference','ADRC controller','Robust ADP controller','Proposed method','Fontangle','italic','Fontweight','bold','Fontsize',16);
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
ylabel('Δx [m]');
xlabel('time [s]');
set (gcf,'position',[500,500,600,300] )

% 
figure(2)
plot(t,qd(2,:),'k--','LineWidth',1.5)
hold on 
plot(t,q3(2,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on 
plot(t,q2(2,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
hold on 
plot(t,q(2,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
legend('Reference','ADRC controller','Robust ADP controller','Proposed method','Fontangle','italic','Fontweight','bold','Fontsize',16);
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
ylabel('Δy [m]');
xlabel('time [s]');
set (gcf,'position',[500,500,600,300] )

figure(3)
plot(t,qd(3,:),'k--','LineWidth',1.5)
hold on 
plot(t,q3(3,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on 
plot(t,q2(3,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
hold on 
plot(t,q(3,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
legend('Reference','ADRC controller','Robust ADP controller','Proposed method','Fontangle','italic','Fontweight','bold','Fontsize',16);
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
ylabel('ΔL [m]');
xlabel('time [s]');
set (gcf,'position',[500,500,600,300] )


figure(4)
plot(t,dqd(1,:),'k--','LineWidth',1.5)
hold on 
plot(t,dq3(1,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on 
plot(t,dq2(1,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
hold on 
plot(t,dq(1,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
legend('Reference','ADRC controller','Robust ADP controller','Proposed method','Fontangle','italic','Fontweight','bold','Fontsize',16);
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
ylabel('dΔx [m]');
xlabel('time [s]');
set (gcf,'position',[500,500,600,300] )

figure(5)
plot(t,dqd(2,:),'k--','LineWidth',1.5)
hold on 
plot(t,dq3(2,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on 
plot(t,dq2(2,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
hold on 
plot(t,dq(2,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
legend('Reference','ADRC controller','Robust ADP controller','Proposed method','Fontangle','italic','Fontweight','bold','Fontsize',16);
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
ylabel('dΔy [m]');
xlabel('time [s]');
set (gcf,'position',[500,500,600,300] )

figure(6)
plot(t,dqd(3,:),'k--','LineWidth',1.5)
hold on 
plot(t,dq3(3,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on 
plot(t,dq2(3,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
hold on 
plot(t,dq(3,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
legend('Reference','ADRC controller','Robust ADP controller','Proposed method','Fontangle','italic','Fontweight','bold','Fontsize',16);
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
ylabel('dΔL [m]');
xlabel('time [s]');
set (gcf,'position',[500,500,600,300] )
% 
figure(7)
plot(t,tao3(1,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on 
plot(t,tao2(1,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
hold on 
plot(t,tao(1,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
legend('ADRC controller','Robust ADP controller','Proposed method','Fontangle','italic','Fontweight','bold','Fontsize',16);
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
set (gcf,'position',[500,500,600,300] )
ylabel('\taux [N/Nm]','Fontangle','italic');
xlabel('time [s]');
axis([-inf inf -5 5])

figure(8)
plot(t,tao3(2,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on 
plot(t,tao2(2,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
hold on 
plot(t,tao(2,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
legend('ADRC controller','Robust ADP controller','Proposed method','Fontangle','italic','Fontweight','bold','Fontsize',16);
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
set (gcf,'position',[500,500,600,300] )
ylabel('\tauy [N/Nm]','Fontangle','italic');
xlabel('time [s]');
axis([-inf inf -5 5])

figure(9)
plot(t,tao3(3,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on 
plot(t,tao2(3,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
hold on 
plot(t,tao(3,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
legend('ADRC controller','Robust ADP controller','Proposed method','Fontangle','italic','Fontweight','bold','Fontsize',16);
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
set (gcf,'position',[500,500,600,300] )
ylabel('fz [N/Nm]','Fontangle','italic');
xlabel('time [s]');
axis([-inf inf -5 5])

figure(10)
plot(t,e3(1,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on 
plot(t,e2(1,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
hold on 
plot(t,e(1,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
legend('ADRC controller','Robust ADP controller','Proposed method','Fontangle','italic','Fontweight','bold','Fontsize',16);
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
ylabel('Errors-Δx [m]');
xlabel('time [s]');
set (gcf,'position',[500,500,600,300] )

figure(11)
plot(t,e3(2,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on 
plot(t,e2(2,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
hold on 
plot(t,e(2,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
legend('ADRC controller','Robust ADP controller','Proposed method','Fontangle','italic','Fontweight','bold','Fontsize',16);
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
ylabel('Errors-Δy [m]');
xlabel('time [s]');
set (gcf,'position',[500,500,600,300] )

figure(12)
plot(t,e3(3,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on 
plot(t,e2(3,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
hold on 
plot(t,e(3,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
legend('ADRC controller','Robust ADP controller','Proposed method','Fontangle','italic','Fontweight','bold','Fontsize',16);
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
ylabel('Errors-ΔL [m]');
xlabel('time [s]');
set (gcf,'position',[500,500,600,300] )






figure(13)
plot3(pos1(1,1),pos1(2,1),pos1(3,1),'k.','markersize',20)
hold on
plot3(pos1(1,N),pos1(2,N),pos1(3,N),'k*','markersize',5)
grid on
plot3(posd(1,:),posd(2,:),posd(3,:),'k--','LineWidth',1.5)
hold on
plot3(pos3(1,:),pos3(2,:),pos3(3,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on
plot3(pos2(1,:),pos2(2,:),pos2(3,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',0.8)
hold on
plot3(pos1(1,:),pos1(2,:),pos1(3,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)

hold on
time=0:0.01:1;
p1=pos1(1,1)*time;Y1=pos1(2,1)*sin(pi*time/2).*time.^2;Z1=pos1(3,1)*time;
pe=pos1(1,N)*time;Ye=pos1(2,N)*sin(pi*time/2).*time.^2;Ze=pos1(3,N)*time;
plot3(p1,Y1,Z1,'--','Color',[140,239,137]/255,'LineWidth',2)
hold on
plot3(pe,Ye,Ze,'--','Color',[0,174,158]/255,'LineWidth',2)
hold on
x11=-0.1:0.01:0.1;x22=zeros(1,21);
y11=-0.1:0.01:0.1;y22=zeros(1,21);
z11=zeros(1,21);
plot3(x22,y11,z11,'k',x11,y22,z11,'k','LineWidth',1.5)
legend('Start','End','Reference','ADRC controller','Robust ADP controller','Proposed method','Fontangle','italic','Fontweight','bold','Fontsize',12)
set(gca, 'Fontname', 'Times New Roman')

figure(14)
plot(t,w(1,:),t,w(2,:),t,w(3,:),t,w(4,:),t,w(5,:),t,w(6,:),'LineWidth',2)
legend('W1','W2','W3','W4','W5','W6','Fontangle','italic','Fontweight','bold','Fontsize',16);
xlabel('time [s]');
ylabel('W','Fontangle','italic');
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
set (gcf,'position',[500,500,600,300] )


figure(15)
plot(t,V,'LineWidth',2)
xlabel('time [s]');
ylabel('$\hat{J}$','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','Fontsize',12)
set (gcf,'position',[500,500,600,300] )







