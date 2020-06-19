Cclc
close all;
clear
%----------------------------System parameter------------------------------%
% In a certain period of time, five vehicle clusters are involved, and there
% are only one active CM transmitter to communicate with their CH receiver in each cluster. 
%----------------------------System parameter------------------------------%

V1=[36,30,38,40,30];% The speed vector of active CM transmitters in all clusters.
V2=[30,37,32,36,39];% The speed vector of CH receivers in all clusters.
v=zeros(5,5);% Initialize the relative velocity matrix between vehicles.
h=zeros(5,5);% Initialize the mobile links¡¯channel fast fading component in  previous time.
G=zeros(5,5);% Initialize average channel gain matrix between vehicles.
I=zeros(5,5);% Initialize interference channel gain matrix between vehicles. 
%(All diagonal elements are 0 elements, because diagonal elements are effective channel gains)
for i=1:5
    for j=1:5
        v(i,j)=abs(V1(i)-V2(j));
        h(i,j)=0.4+0.9*rand(1);
    end
end
% [6,1,4,0,3;0,7,2,6,9;8,1,6,2,1;10,3,8,4,1;0,7,2,6,9];Relative velocity matrix
h=[0.845418353066450,0.995224514302217,0.687960311535798,0.741632431164121,0.706496751675251;0.716599213586070,0.675355816604331,0.609143849450247,0.779761672779394,0.697456077357881;0.947490617895953,0.811444306977048,0.965654067256766,0.989572070823563,0.834170383464932;0.647590153069472,0.970613035572538,0.837424308992463,0.953446121969907,0.769790372435062;0.842902954353616,0.628305432990921,0.969908972473279,0.856831766456680,0.641799867159030];
Delta=1e-6;%Background noise
T=0.001;%CSI sampling period 1ms
c=3*1e+8;%
fc=5.9*1e+9;%Carrier frequency
j1=2*pi*fc*T/c*v;%The parameter of zero-order Bessel function
epsi=besselj(0,j1);%Calculate the value of Bessel function
a=(epsi.^2).*h;% CSI feedback part in the fast fading model
d=distance(5);%Relative distance matrix
Shadow=[2.16,0.1,0.1,1.28,1.76;0.1,2.16,0.2,1,2.;1,0.2,2.1,0.1,0.1;1.28,0.1,0.1,2.2,0.21;5.48,2.01,1.06,0.21,2.16];%
%Shadow fading matrix£» The value of element is set as 0.1-10.
%The parameter selection principle is to keep the large-scale fading portion 
%of the interference channel gain as consistent as possible, 
%so that it is easy to make balanced interference management 
%at the same interference threshold.
  theta=3;%Pathloss index
  l=Shadow.*d.^(-theta);%Large scale fading matrix
  G1=l.*a;
  G2=l.*(1-epsi.^2);%Average error channel gain matrix 
  G=G1+G2;
  MU=0.5;%The value of Mu in Bernstein approximation
  SIGMA=0.5;%The value of SIGMA in Bernstein approximation
  alfaj=l.*(1-epsi.^2);%The matrix of alfa
  betaj=l.*(1-epsi.^2);%The matrix of beta
 X=G1+MU*alfaj+betaj;%%The matrix of xi 
      E=0.1;          % Outage probability threshold;
  I=X+SIGMA*sqrt(-2*log(E))*alfaj; 
% --Constructing the interference coefficient matrix of Bernstein approximation--% 
    for i=1:5
      I(i,i)=0;
    end
    
global N
N =25;
arfa=zeros(N,5);%Successive convex approximation coefficient matrix  X
beta=zeros(N,5);%Successive convex approximation coefficient matrix  Y
SINR=zeros(N,5);
Pmax=log(0.5)*ones(N,5);
B=0.01;%Spectrum bandwidth(GHz)£¬also can be set as 1 or 0.001 which 
%affect the final energy efficiency convergence speed and precision. 
C1=linspace(1,1,5);
GF=linspace(1.5,1.5,5);%Power amplifier coefficient vector.
PC=0.05;%Fixed power consumption in the circuit
 S=zeros(N,5);
 s=zeros(N,5);
Temp=zeros(40,5);%Optimal power storage area for each iteration
 for F=1:10
 C=0.5*F;%Initial fixed price
 Ith=1e-6;%Interference threshold
 Z=1;
 P=zeros(N,5);
P(1,:)=[-8,-8,-8,-8,-8]; 
Lamda=zeros(N,5);
Lamda(1,:)=linspace(10,10,5);
for k=1:N
   for i=1:5
        I1(i)=exp(P(k,:))*G(:,i)+Delta-G(i,i)*exp(P(k,i));
        SINR(k,i)=G(i,i)*exp(P(k,i))/I1(i); 
        arfa(k,i)=SINR(k,i)/(1+SINR(k,i));
        Total_T(k)=B*dot(log(C1+SINR(k,:))/log(2),C1);
   end
  
       for i=1:5
        for j=1:5
            if j~=i
                S(k,j)=B/log(2)*arfa(k,j)*G(i,j)*SINR(k,j)/((exp(P(k,j))*G(j,j)))+Lamda(k,j)*I(i,j);%
                 else
             S(k,j)=0;
            end            
        end
        s(k,i)=sum(S(k,:))+C*GF(i);
  end
%  Total_T(k+1)=B*dot(log(C1+SINR(k,:))/log(2),C1); 
  P(k+1,:)=min(log(B.*arfa(k,:)/log(2))-log(s(k,:)),Pmax(k,:));
    if max(abs(P(k+1,:)-P(k,:)))<=1e-2
       P(k+1,:)=P(k,:);
    end
    for i=1:5
     Lamda(k+1,i)=max(0,Lamda(k,i)+100000/k*(exp(P(k,:))*I(:,i)-Ith));
     w(k,i)=exp(P(k,:))*I(:,i)-Ith;
    end 
    if  max(abs(w(k,:)))<=3e-7   
   Lamda(k+1,:)=Lamda(k,:);
    end
end
  SE(F)=Total_T(N);
  Temp(F,:)=P(N,:);
  EE(F)=SE(F)/(GF*exp(Temp(F,:)')+PC);
end
EE1=100*EE;
SE1=100*SE;

C=50+linspace(0,450,10);
figure
 [AX,H1,H2] = plotyy(C,EE1,C,SE1,'plot');
set(AX(1),'XColor','k','YColor','b','LineWidth',1.5);
set(AX(2),'XColor','k','YColor','r','LineWidth',1.5);
HH1=get(AX(1),'Ylabel');
set(HH1,'String','Average EE(bits/Hz/Joule)','fontsize',12);
set(HH1,'color','b');
HH2=get(AX(2),'Ylabel');
set(HH2,'String','Average SE(bits/s/Hz)','fontsize',12);
set(HH2,'color','r');
set(H1,'marker','s')
set(H1,'LineStyle','-','LineWidth',1.5);
set(H1,'color','b');
set(H2,'LineStyle','-','LineWidth',1.5);
set(H2,'color','r');
set(H2,'marker','^')
legend([H1,H2],{'EE';'SE'});
xlabel('Number of clusters','fontsize',12);
xlabel('Fixed price C','fontsize',12);
legend('EE','SE')
