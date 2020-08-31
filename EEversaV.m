clc
close all
clear
%----------------------------System parameter------------------------------%
%Fig. 6(b) in this paper.
% In a certain period of time, five vehicle clusters are involved, and there
% are only one active CM transmitter to communicate with their CH receiver in each cluster. 
%----------------------------System parameter------------------------------%
  L=zeros(5,5);
for F=1:6
V1=18*ones(1,5)+2*F*ones(1,5);%The speed vector of active CM transmitters in all clusters
V2=[20,20,20,20,20];%The speed vector of CH receivers in all clusters.
v=zeros(5,5);%Initialize the relative velocity matrix between vehicles.
h=zeros(5,5);%Initialize the mobile links¡¯channel fast fading component in  previous time.
G=zeros(5,5);%%Initialize average channel gain matrix between vehicles.
I=zeros(5,5);% Initialize interference channel gain matrix between vehicles. 
%(All diagonal elements are 0 elements, because diagonal elements are effective channel gains)
for i=1:5
    for j=1:5
         v(i,j)=abs(V1(i)-V2(j));
         h(i,j)=0.9+0.4*rand(1);
    end
end
h=[1.21026765186132,1.00186773622345,1.11230368038275,1.27044868578780,1.09550344842318;1.17912487696444,1.08933855474170,1.03750431438572,1.11650670059712,1.24530632319216;1.29435271127300,1.25859695058615,1.02514626462880,1.10131359710114,1.07083874991164;1.09534161637965,1.29533452164231,1.16447526737693,1.22477537780633,1.25255552041723;1.05006691209611,1.27092928278688,1.03153723834667,1.22352792893222,1.21881158135137];
Delta=1e-6;%%Background noise
T=0.001;%CSI sampling period 1ms
c=3*1e+8;%
fc=5.9*1e+9;%Carrier frequency
j1=2*pi*fc*T/c*v;%The parameter of zero-order Bessel function
epsi=besselj(0,j1);%Calculate the value of Bessel function
a=(epsi.^2).*h;% CSI feedback part in the fast fading model
d=distance(5);%Relative distance matrix
Shadow=[2.16,0.1,0.1,1.28,1.76;0.1,2.16,0.2,1,2;1,0.2,2.1,0.1,0.1;1.28,0.1,0.1,2.2,0.21;5.48,2.01,1.06,0.21,2.16];%
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
B=1;%Spectrum bandwidth£¬also can be set as 1 or 0.001(GHz) which 
%affect the final energy efficiency convergence speed and precision. 
C1=linspace(1,1,5);
GF=linspace(1.5,1.5,5);%Power amplifier coefficient vector.
PC=0.05;%%Fixed power consumption in the circuit
 S=zeros(N,5);
 s=zeros(N,5);
 eta=0.1;%%Initial value of energy efficiency
 Ith=1.3e-6;%Interference threshold
 Temp=zeros(40,5);%%Optimal power storage area for each iteration
  Z=1;
while 1
   P=zeros(N,5);
   P(1,:)=[-8,-8,-8,-8,-8]; 
   Lamda=zeros(N,5);
   Lamda(1,:)=linspace(20,20,5);
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
                S(k,j)=B/log(2)*arfa(k,j)*G(i,j)*SINR(k,j)/((exp(P(k,j))*G(j,j)))+Lamda(k,j)*I(i,j);
              else
             S(k,j)=0;
              end
          end
          s(k,i)=sum(S(k,:))+eta(Z)*GF(i);
  
  end
%  Total_T(k+1)=B*dot(log(C1+SINR(k,:))/log(2),C1); 
P(k+1,:)=min(log(B.*arfa(k,:)/log(2))-log(s(k,:)),Pmax(k,:));
   if max(abs(exp(P(k+1,:))-exp(P(k,:))))<=1e-5
  P(k+1,:)=P(k,:);
   end
    for i=1:5
   Lamda(k+1,i)=max(0,Lamda(k,i)+1000000/k*(exp(P(k,:))*I(:,i)-Ith));
  w(k,i)=exp(P(k,:))*I(:,i)-Ith;
    end 
    if  max(abs(w(k,:)))<=1.2e-7   
   Lamda(k+1,:)=Lamda(k,:);
    end

end
  D(Z)=Total_T(N);
  Temp(Z,:)=P(N,:);
  error =abs(D(Z)-eta(Z)*(GF*exp(Temp(Z,:)')+PC));
  eta(Z+1)=D(Z)/(GF*exp(Temp(Z,:)')+PC);
      if  error<=1e-3
    L(F,:)=exp(Temp(Z,:));
    Q(F)=D(Z);
    eta1(F)=eta(Z+1);
           break;
      end 
     Z=Z+1;
end
  end
  Total_T1=Q;
plot(exp(P(1:25,1)),'-+r');
hold on
plot(exp(P(1:25,2)),'-og');
plot(exp(P(1:25,3)),'-*b');
plot(exp(P(1:25,4)),'-sk');
plot(exp(P(1:25,5)),'-dc');
legend('P1','P2','P3','P4','P5')
xlabel('iteration');
ylabel('Power conversion value(W)');

figure
plot(Lamda(1:25,1),'-+r');
hold on
plot(Lamda(1:25,2),'-og');
plot(Lamda(1:25,3),'-*b');
plot(Lamda(1:25,4),'-sk');
plot(Lamda(1:25,5),'-dc');
legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5')
xlabel('iteration');
ylabel('Multiplier convergence');

figure
plot(eta,'-sb','linewidth',1.5) 
hold on
plot(D,'-^r','linewidth',1.5)
hold on
[AX,H1,H2] =plotyy(Z,eta,Z,D,@plot);% 
set(AX(1),'XColor','k','YColor','b','LineWidth',1.5);
set(AX(2),'XColor','k','YColor','r','LineWidth',1.5)
set(get(AX(1),'ylabel'),'color','blue','string', 'Average EE(bits/Joule/Hertz)','fontsize',12);
set(get(AX(2),'ylabel'),'color','red','string', 'Average SE(bits/sec/Hertz)','fontsize',12);
xlabel('iteration','fontsize',12);
legend('EE','SE')

 FF=[0,2,4,6,8,10];
 figure
 [AX,H1,H2] = plotyy(FF,eta1,FF,Total_T1,'plot');

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
xlabel('\Delta\upsilon (m/s)','fontsize',12);
legend('EE,C=C^*','SE,C=C^*') 
