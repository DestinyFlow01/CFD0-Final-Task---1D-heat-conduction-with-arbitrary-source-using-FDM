%material used in this task : Aluminum (Al)
clear all; clc
%Material Data
k = 237; rho = 2700; q = 3*10^5; cp = 897; L = 0.1; %in SI unit
alpha = sqrt(rho*cp/k);

%Case : bar clamped on both sides 

%Numerical analysis : 
dx = 0.01; dt = 0.5; C1 = dt/(alpha^2*dx^2);

%Explicit scheme : 
m = L/dx+1; %termasuk dua ujungny
T = zeros(1,m);
T(1) = 212; %Dirichlet BC on the left
T(m) = 152; %Dirichlet BC on the right
x = 0:dx:L;
Taft = T;

%heat generation matrix : 
Q = zeros(1,m);
gen = zeros(1,m);
n = 2; %orde osilasi profil heat generation

for i = 1:m 
   Q(i) = q*sin(n*pi*x(i)/L);
   gen(i) = Q(i)*dt/(k*alpha^2);
end

%plotting both analytical and numerics in some time simultaneeously for t =
%2, 4, 8, 16, 32, 64, 128, 256, steady state
t = 0;
T1 = T(1); T2 = T(m);
y = 10;
koef = zeros(1,y); %first y+1 terms of the transient

for i = 1:y
    first = -2/(pi^2*k*n*i)*(2*q*L^2*(-1)^(i+1)+pi*k*(T1+T2*(-1)^(i+1)));
    second = -q*L/(pi^2*k*i^2);
    koef(i) = first + second;
end

l = 0;
tfinal = 310;
for i = 1:tfinal/dt
    %numerical : 
    Tbef = Taft;
    for j = 2:m-1 
        Taft(j) = Tbef(j)+C1*(Tbef(j+1)-2*Tbef(j)+Tbef(j-1))+gen(j);
    end
    t = t+dt
    
    figure(10)
    plot(x,Taft)
    title('Numerik')
    hold on
    
    %Analitik
    figure(11)
    for e = 2:m-1
            transient = 0;
            steady = q*L/(k*pi^2*n^2)*(sin(n*pi*x(e)/L)-n*pi*x(e))+(T2-T1+q*L^2/(k*pi*n))*x(e)/L+T1;
            for f = 1:y
                transient = transient + koef(f)*exp(-(f*pi/(L*alpha))^2*t)*sin(f*pi*x(e)/L);
            end
            T(e) = transient + steady;
    end
    plot(x,T)
    title('Analitik')
    hold on
    
    if(t==2 | t==4 | t==8 | t == 16|t == 32|t==64 | t == 128 | t==256 |t==tfinal)
        l = l+1;
        figure(l)
        
        plot(x,Taft) %numerik
        hold on
        plot(x,T) %analitik
        legend('Numerik','Analitik')
        hold off
    end
end
hold off 
hold off