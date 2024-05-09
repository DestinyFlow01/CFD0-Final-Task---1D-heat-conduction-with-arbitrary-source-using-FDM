%material used in this task : Aluminum (Al)
clear all; clc
%Material Data
k = 237; rho = 2700; q = 3*10^3; cp = 897; L = 1000; %in SI unit
alpha = sqrt(rho*cp/k);

%Case : long bar with Gaussian heat transfer mode

%Numerical analysis : 
%dx = 2; dt = 0.5; C1 = dt/(alpha^2*dx^2);
dx = 0.1; C1 = 1.22321317973492e-05; dt = C1*(alpha*dx)^2;

%Explicit scheme : 
n = L/dx+1; %termasuk dua ujungny
T = zeros(1,n);
x = -L/2:dx:L/2;
Taft = T;

%heat generation matrix : 
gen = zeros(1,n);
a = q/k*sqrt(pi/2);
C2 = q*dt/(k*alpha^2);

for i = 1:n
   gen(i) = C2*x(i)*exp(-x(i)^2/2);
end

t = 0;
Taft(1)=a*(erf(x(1)/sqrt(2))+erf(x(1)/sqrt(2+4*t/alpha^2)));
Taft(n)=a*(erf(x(n)/sqrt(2))+erf(x(n)/sqrt(2+4*t/alpha^2)));
for i = 1:50000/dt
    
    Tbef = Taft;
    for j = 2:n-1
        Taft(j) = Tbef(j)+C1*(Tbef(j+1)-2*Tbef(j)+Tbef(j-1))+gen(j);
    end
    t = t+dt;
    Taft(1)=a*(erf(x(1)/sqrt(2))+erf(x(1)/sqrt(2+4*t/alpha^2)));
    Taft(n)=a*(erf(x(n)/sqrt(2))+erf(x(n)/sqrt(2+4*t/alpha^2)));
    if(t==1200)
        Taft1200 = Taft;
        
    elseif(t==2400)
        Taft2400 = Taft;
            
    elseif(t==3600)
        Taft3600 = Taft;
        
    elseif(t==4800)
        Taft4800 = Taft;
    end
end


%plotting analytical and numerical result as comparison : 
%temperature : 

T = zeros(1,n);

%steady state : 
for i = 1:n
   T(i)=a*erf(x(i)/sqrt(2)); 
end

%{
figure(1)
plot(x,T);
hold on
plot(x,Taft1200)
plot(x,Taft2400)
plot(x,Taft3600)
plot(x,Taft4800)
plot(x,Taft)
title('Steady State compared with numerical')
legend('Steady State','t = 20 min','t = 40 min','t = 60 min','t = 80 min','t = 5E4 s');
hold off

figure(2) 
plot(x,Taft1200)
hold on
for i = 1:n
   T(i)=a*(erf(x(i)/sqrt(2))-erf(x(i)/sqrt(2+4*1200/alpha^2))); 
end
plot(x,T);
title('t = 20 min')
legend('Numerical','Analytical')
hold off

figure(3) 
plot(x,Taft2400)
hold on
for i = 1:n
   T(i)=a*(erf(x(i)/sqrt(2))-erf(x(i)/sqrt(2+4*2400/alpha^2))); 
end
plot(x,T);
title('t = 40 min')
legend('Numerical','Analytical')
hold off

figure(4) 
plot(x,Taft3600)
hold on
for i = 1:n
   T(i)=a*(erf(x(i)/sqrt(2))-erf(x(i)/sqrt(2+4*3600/alpha^2))); 
end
plot(x,T);
title('t = 60 min')
legend('Numerical','Analytical')
hold off

figure(5) 
plot(x,Taft4800)
hold on
for i = 1:n
   T(i)=a*(erf(x(i)/sqrt(2))-erf(x(i)/sqrt(2+4*4800/alpha^2))); 
end
plot(x,T);
title('t = 80 min')
legend('Numerical','Analytical')
hold off
%}
figure(6) 
plot(x,Taft)
hold on
for i = 1:n
   T(i)=a*(erf(x(i)/sqrt(2))-erf(x(i)/sqrt(2+4*500000/alpha^2))); 
end
plot(x,T);
title('t = 5E4 s')
legend('Numerical','Analytical')
hold off