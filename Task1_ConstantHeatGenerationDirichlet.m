%material used in this task : Aluminum (Al)
clear all; clc
%Material Data
k = 237; rho = 2700; q = 3*10^5; cp = 897; L = 0.1; %in SI unit
alpha = sqrt(rho*cp/k);

%Case : bar clamped on both sides with 

%Numerical analysis : 
dx = 0.01; dt = 0.5; C1 = dt/(alpha^2*dx^2);
gen = q*dt/(k*alpha^2);

%Explicit scheme : 
n = L/dx+1; %termasuk dua ujungny
T = zeros(1,n);
T(1) = 212; %Dirichlet BC on the left
T(n) = 152; %Dirichlet BC on the right
x = 0:L/(n-1):L;
Taft = T;
figure(1)
t = 0;
for i = 1:160
    Tbef = Taft;
    for j = 2:n-1
        Taft(j) = Tbef(j)+C1*(Tbef(j+1)-2*Tbef(j)+Tbef(j-1))+gen;
    end
    t = t+dt;
    if(t==10 | t==40 | t==70) 
        plot(x,Taft)
        hold on
    end
    
end
hold off
title('Numerical')
legend('t = 10','t = 40','t = 70')

%plotting analytical and numerical result as comparison : 
T = -q.*x.^2/(2*k) + (T(n)-T(1)+q*L^2/(2*k)).*x/L+T(1);
figure(2)
plot(x,T);
hold on
plot(x,Taft);
hold off
legend('Analytical result','Numerical result')
title('Steady State condition comparison')