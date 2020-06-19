%%%--------------Generalization of Topological Coordinates------------%%%%%
%%Suppose the width of the two lanes is 8m, where
%%CM_1(6,0),CM_2(6,8),CM_3(22,0),CM_4(22,8),CM_5(38,0);
%%CH_1(0,0),CH_2(12,8),CH_3(16,0),CH_4(28,8),CH_5(32,0);
function d=distance(M)

for n=1:M
    x(n)=8*n-6+2*(-1)^n;
    y(n)=4+4*(-1)^n;%CH receiver's X-Y coordinates
    a(n)=8*n-6+2*(-1)^n+6*(-1)^(n+1);
    b(n)=4+4*(-1)^n;%CM transmitter's X-Y coordinates
end
    for i=1:M
        for j=1:M
            d(i,j)=sqrt((a(i)-x(j))^2+(b(i)-y(j))^2);%Relative distance between two vehicles
        end
    end