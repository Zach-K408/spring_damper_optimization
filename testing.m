% Testing Parameters

% intial conditions
x0 = 1;
v0 = 0;

% time domain
t_start = 0;
t_end = 10;
t_step = 0.1;
t = t_start:t_step:t_end;

% forcing function
f0 = 0;
f = @(t) f0 + 0*t; 

% system properties
m = 100; % kg

% decision variables
c_min = 0;
c_max = 20;
c_step = 5;
c_vals = [0,5,15]; % kg/s

k = 15; % kg/s^2

% solver tolerances
alpha = odeset('RelTol',1e-6,'AbsTol',1e-6);

%function call

xc = zeros(length(c_vals),length(t));     % store position damping pair values
vc = 0*xc;                            % store velocity damping pair values

for i=1:length(c_vals)
    [~,xc(i,:),vc(i,:)] = function_eval(m,c_vals(i),k,f,t,x0,v0,alpha);
end

% plot forcing function
figure
sgtitle('Response to varying damping coefficents and constant forcing function','fontsize',14,'Interpreter','latex')

subplot(3,1,1)
fplot(f,[t_start t_end],'linewidth',2)
grid on; hold on
ylabel('Forcing Function (N)','fontsize',14,'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

% plot position
subplot(3,1,2)
plot(t,xc,'linewidth',2)
grid on; hold on
ylabel('Position (m)','fontsize',14,'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

% plot velocity
subplot(3,1,3)
plot(t,vc,'linewidth',2)
grid on; hold on
xlabel('Time (s)','fontsize',14,'Interpreter','latex')
ylabel('Velocity $\left( \frac{m}{s} \right)$','fontsize',14,'Interpreter','latex')

legend([repmat('$c= $ ',[length(c_vals) 1]) num2str(c_vals')],...
    'fontsize',14,'Interpreter','latex','location','southeast')

set(gca,'TickLabelInterpreter','latex')
set(gcf,'Position',[988 196 881 714])