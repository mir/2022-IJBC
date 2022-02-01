clear all;
close all;
clc;

% VCO input gain
K_vco = 250;

% Loop filter transfer function F(s) = (1+s tau_2)/(s tau_1), tau_1 > 0,
% tau_2 > 0
tau_2 = 0.0225;
tau_1 = 0.0633;

% The slope coefficient of piecewise-linear PD characteristic
% it becomes triangular with this coefficient
k=2/pi;

function y = draw_saddles_symmetric(ode, saddle, V, periods, options)
% draw_saddles_symmetric draws phase portrait of the system ode, 
% starting from saddle with eigenvectors V for given number of periods
% options are used to tune ODE solver

    % Integration time
    TIME_POSITIVE = 0:0.005:20;
    TIME_NEGATIVE = 0:-0.005:-20;
   
    v1 = 0.01 .* flip(V(:,1)');
    v2 = 0.01 .* flip(V(:,2)');
    
    % Calculate saddle separatrices and plot them
    [T_1, X_1] = ode45(ode, TIME_NEGATIVE, saddle+v1, options);
    [T_3, X_3] = ode45(ode, TIME_NEGATIVE, saddle-v1, options);
    for j=1:length(periods)
        plot(X_1(:,2)+periods(j),X_1(:,1));
        plot(X_3(:,2)+periods(j),X_3(:,1));
    end
    [T_2, X_2] = ode45(ode, TIME_POSITIVE, saddle+v2, options);
    [T_4, X_4] = ode45(ode, TIME_POSITIVE, saddle-v2, options);
    for j=1:length(periods)
       plot(X_2(:,2)+periods(j), X_2(:,1));
       plot(X_4(:,2)+periods(j), X_4(:,1));
    end
end

function y = sawtooth_diff(t)
% sawtooth_diff - the derivative of a triangular function 
% (sawtooth (T, WIDTH) from signal package)
  remain = abs(rem(t, 2*pi));
  if (pi/2 <= remain && remain <= 3*pi/2)
      y = -2/pi;
  else
      y = 2/pi;
  end
end

% PD characteristic
v_e = @(theta_e) (sawtooth(theta_e+pi/2,0.5));
% The derivative of PD characteristic
dv_e = @(theta_e) (sawtooth_diff(theta_e));
period = 2*pi;

% Parameters a, b, c from Theorem 1 and 2
a = sqrt(K_vco/tau_1)*tau_2;
b = sqrt(abs(a^2-4/k));
c = sqrt(a^2-4/k + 4*pi);
% Computing the lock-in range and the concervative lock-in range by 
% Theorem 1 and 2
syms x;
if a^2*k>4
    y_l = sqrt(pi)* ((c+b)/(c-b))^(a/(2*b));
    fcn = @(x) (x - (a-b)/2)^((b-a)/b) * (x - (a+b)/2)^((b+a)/b) -...
        pi*((c+b)/(c-b))^(a/b);
    init_param = [(a+b)/2, 10000000];
    % vpasolve numerically solves implicit equations with initial guess
    % init_param (it was proven that the equation has a unique solution
    % for x > (a+b)/2)
    d = vpasolve(fcn(x), x, init_param);
else 
    if a^2*k == 4
        y_l = sqrt(pi)*exp(a/2/sqrt(pi));
        d = a/2*(1 + 1/lambertw(a/2/sqrt(pi)*exp(-a/2/sqrt(pi))));
else
        y_l = sqrt(pi)*exp(a*atan(b/c)/b);
        fcn = @(x) ((x)^2 - a*x + 1/k) * exp( 2*a/b*atan((2*x-a)/b) - pi*a/b) -...
        pi*exp(2*a/b*atan(b/c));
        init_param = [a/2, 10000000];
        % vpasolve numerically solves implicit equations with initial guess
        % init_param (it was proven that the equation has a unique solution
        % for x > a/2)
        d = vpasolve(fcn(x), x, init_param);
    end
end
y_l_c = (d-0.5*(a-c))^((c-a)/c/2) * (d-0.5*(a+c))^((c+a)/c/2);


h = figure(1);
hold on;
grid on;

% Ploting two points which correspond the lock-in range and the conservative
% lock-in range
plot([0], [y_l], 'k.', 'MarkerSize', 20);
plot([-pi], [eval(y_l_c)], 'k.', 'MarkerSize', 20);

% System establishment for numerocal integration
% y = x(1), theta_e = x(2)
pll_s = @(t,x)([- a*dv_e(x(2))*x(1) - v_e(x(2));
                 x(1)]);

% One of the asymptotically stable equilibria
theta_eq = 0;
x_eq = 0;
    
% Draw phase portrait
saddles = [x_eq, -theta_eq+period/2-2*period;...
    x_eq, -theta_eq+period/2-period;...
    x_eq, -theta_eq+period/2;...
    x_eq, -theta_eq+period/2+period;...
    x_eq, -theta_eq+period/2+2*period];
focuses = [x_eq, theta_eq-period;...
    x_eq, theta_eq-period;...
    x_eq, theta_eq;...
    x_eq, theta_eq+period;...
    x_eq, theta_eq+2*period];
    
plot(-focuses(:,2), -focuses(:,1), 'r.', 'MarkerSize', 20);
plot(focuses(:,2), focuses(:,1), 'k.', 'MarkerSize', 20);
plot(-saddles(:,2), -saddles(:,1), 'r.', 'MarkerSize', 20);
plot(saddles(:,2), saddles(:,1), 'k.', 'MarkerSize', 20);

% Jacobian matrix of the system
A = [0 1; 
    2/pi 2*a/pi];

% Calculating saddle eigenvectors V
[V, D] = eig(A);

% Custom simulation options
options = odeset('MaxStep', 0.001, 'RelTol', 2e-7, 'AbsTol', 2e-7);
draw_saddles_symmetric(pll_s,...
    [x_eq,period/2-theta_eq],  ...
    V, ...
    [-2*period,-period,0,period,2*period],...
    options);

% Plot adjustments
axis([-2*pi 2*pi -5 5])
xticks([-4*pi -3*pi -2*pi -pi 0 pi 2*pi, 4*pi])
xticklabels({'-4\pi','-3\pi','-2\pi','-\pi','0','\pi','2\pi','4\pi'})

xlabel('\theta_e');
ylabel('y');




