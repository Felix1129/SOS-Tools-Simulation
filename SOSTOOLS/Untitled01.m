clear
clc
[t,x] = ode45(@vdp1,[0 2],[0.15; 0.15]);
subplot(2,1,1),plot(t,x(:,1),'-'),xlabel('time'),ylabel('X1(t)'),grid on;
subplot(2,1,2),plot(t,x(:,2),'-'),xlabel('time'),ylabel('X2(t)'),grid on;
function dydt = vdp1(t,x)
dydt = [-(7/2+3/2*sin(x(1)))*x(1)-4*x(2);
        (19/2-21/2*sin(x(1)))*x(1)-2*x(2)];
end