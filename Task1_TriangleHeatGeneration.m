%material used in this task : Aluminum (Al)
clear all; clc
%Material Data
k = 237; rho = 2700; q = 3*10^5; cp = 897; L = 0.1; %in SI unit
alpha = sqrt(rho*cp/k);

%Case : bar clamped on both sides with 

%Numerical analysis : 
dx = 0.01; dt = 0.5; C1 = dt/(alpha^2*dx^2);

%Explicit scheme : 
n = L/dx+1; %termasuk dua ujungny
T = zeros(1,n);
T(1) = 212; %Dirichlet BC on the left
T(n) = 152; %Dirichlet BC on the right
x = 0:L/(n-1):L;
Taft = T;

%heat generation matrix : 
Q = zeros(1,n);
gen = zeros(1,n);

for i = 1:n 
   if(i<=(n-1)/2) 
      Q(i) = 2*q*x(i)/L;
   else
      Q(i) = 2*q/L*(L-x(i));
   end
   gen(i) = Q(i)*dt/(k*alpha^2);
end

figure(1)
t = 0;
for i = 1:160
    Tbef = Taft;
    for j = 2:n-1
        Taft(j) = Tbef(j)+C1*(Tbef(j+1)-2*Tbef(j)+Tbef(j-1))+gen(j);
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
%temperature : 
a = q*L^2/k; 
A = [0 1 1;-0.5 0.5 1;1 -1 0];
T0 = T(1); T1 = T(n);
b = [T1+2*a/3; T0+a/6; 0.5*a];
koef = inv(A)*b;
C1 = koef(1)/L; C3 = koef(2)/L; C4 = koef(3);

for i = 1:n
   if(i <=floor(n/2))
       T(i) = -q*x(i)^3/(3*k*L)+T0 + C1*x(i);
   else
       T(i) = q*(x(i)^3/3-L*x(i)^2)/(k*L)+C3*x(i)+C4;
   end
end

%plotting :
figure(2)
plot(x,T);
hold on
plot(x,Taft);
hold off
legend('Analytical result','Numerical result')
title('Steady State condition comparison')