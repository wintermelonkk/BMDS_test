clear
close all
clc

k1 = 100;% unit: /muM/min
k2 = 600; % unit: /min
k3 = 150; % unit:/min
x0 = [1;10;0;0]; % initial value for E=1uM S= 10uM ES=0uM P=0uM
h = 0.0001; % step size
tlist = 0:h:1; % total span for final demonstration is 1s
xsol(:,1) = x0;

%% Runge-Kutta function defintion
func = @(t,x) [-k1*x(1)*x(2)+k2*x(3)+k3*x(3);
            -k1*x(1)*x(2)+k2*x(3);
            -k2*x(3)+k1*x(1)*x(2)-k3*x(1)*x(4);
            k3*x(3)];

%% Runge-Kutta methods for question 2
ysol(:,1) = func(tlist(1),x0);
for i = 2:length(tlist)
    k1 = h * func(tlist(i),x0);
    k2 = h * func(tlist(i)+h/2,x0 + k1/2);
    k3 = h * func(tlist(i)+h/2,x0 +k2/2);
    k4 = h * func(tlist(i)+h,x0+k3);
    x1 = x0 + 1/6*(k1+2*k2+2*k3+k4);
    xsol(:,i) = x1;ysol(:,i) = func(tlist(i),x1);
    x0 = x1;
end


%% plot for question 2
figure
hold on
plot(tlist,xsol(1,:)');
plot(tlist,xsol(2,:)');
plot(tlist,xsol(3,:)');
plot(tlist,xsol(4,:)');
legend('E','S','ES','P');
ylabel('C_S, C_ES,C_S,C_P');
xlabel('time in sec');
%% plot for question 3
[~,ind] = max(ysol(4,:));
figure
plot(xsol(2,:),ysol(4,:),'r',xsol(2,ind),ysol(4,ind),'bo');
xlabel('C_S');
ylabel('V_P');

