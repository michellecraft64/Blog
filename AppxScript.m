%Initial Value Problem (IVP)
%ODE y'=sin(t)^2*y solved by subfunction
h=0.01; %step-size (adjustable) AKA change in time or delta t (dt)
tEnd=5; %End of time interval (in s)
t=0:h:tEnd; %time vector
Lt=length(t); %length of time vector defined
y=2*exp(0.5*(t-sin(t).*cos(t))); %Exact solution
%Initializing vectors
y_FE=zeros(Lt,1); %Forward Euler y vector
y_FE(1)=y(1); %Initial value
y_RK4=zeros(Lt,1); %Runge-Kutta y vector
y_RK4(1)=y(1); %Initial value
%Recursive loop
for j=2:Lt %starting at second value
    y_FE(j)=y_FE(j-1)+h*fnc(t(j-1),y_FE(j-1)); %Forward Euler
    rhs1=fnc(t(j-1),y_RK4(j-1)); %First right hand side (RHS) of RK4 AKA k1
    rhs2=fnc(t(j-1)+0.5*h,y_RK4(j-1)+h*0.5*rhs1); %k2
    rhs3=fnc(t(j-1)+0.5*h,y_RK4(j-1)+h*0.5*rhs2); %k3
    rhs4=fnc(t(j-1)+h,y_RK4(j-1)+h*rhs3); %k4
    y_RK4(j)=y_RK4(j-1)+(1/6)*h*(rhs1+2*rhs2+2*rhs3+rhs4);
end
%Plot comparison of each solution
figure
hold on
plot(t,y_FE)
plot(t,y_RK4)
plot(t,y)
hold off
legend('Forward Euler','RK4','Exact')
xlabel('Time (s)')
ylabel('Y')
title('Approximation Comparison')
function x=fnc(t,y)
    x=sin(t)^2*y;
end