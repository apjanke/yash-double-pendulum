function double_pendulum
% Make some plots of double pendulum behavior.
%
%

% TODO: What do all these short variable names mean?
l1 = 1;
l2 = 2;
m1 = 2;
m2 = 1;
g = 9.8;

% Initial condition

tspan = 50;
theta1 = 1.6;
theta1_prime = 0;
theta2 = 2.2;
theta2_prime = 0;

y0 = [theta1 theta1_prime theta2 theta2_prime];
[t, y] = ode45(@pend, [0 ,tspan], [2.5 0 1 0]);


% Position of mass 1 and mass 2

x1 = l1 * sin(y(:,1));
y1 = -l1 * cos(y(:,1));
x2 = l1 * sin(y(:,1)) + l2 * sin(y(:,3));
y2 = -l1 * cos(y(:,1)) - l2 * cos(y(:,3));

% Visualizing the result

fig1 = figure;
plot(x1, y1, 'linewidth', 2);
hold on
plot(x2, y2, 'r', 'linewidth', 2);
ax = gca;
set(ax, 'fontSize', 14);
xlabel('X', 'fontSize', 14);
ylabel('Y', 'fontSize', 14);
title('Chaotic Double Pendulum','fontsize', 14);
set(fig1, 'color', 'white');

fig2 = figure;
plot(y(:,1), 'linewidth', 2);
hold on
plot(y(:,3), 'r', 'linewidth', 2);
ax = gca;
set(ax, 'fontSize', 14);
legend('\theta_1', '\theta_2');
xlabel('time', 'fontSize', 14);
ylabel('theta', 'fontSize', 14);
title('\theta_1(t=0)=2.5 and \theta_2(t=0)=1.0','fontsize', 14);
set(fig2, 'color', 'white');


% Movie of double pendulum

fig3 = figure;
nFrames = 0;
fram = 0;
ax = gca;
set(fig3, 'color', 'white');

for i = 1:numel(y)
    nFrames = nFrames + 1;
    fram = fram + 1;
    plot(ax, 0, 0, '.', 'markersize', 20);
    hold on
    plot(ax, x1(i), y1(i), '.', 'markersize', 20);
    plot(ax, x2(i), y2(i), '.', 'markersize', 20);
    hold off
    line(ax, [0 x1(i)], [0 y1(i)], 'Linewidth', 2);
    axis(ax, [-(l1+l2) l1+l2 -(l1+l2) l1+l2]);
    line(ax, [x1(i) x2(i)], [y1(i) y2(i)], 'linewidth', 2);
    set(ax, 'fontSize', 12);
    xlabel(ax, 'X', 'fontSize', 12);
    ylabel(ax, 'Y', 'fontSize', 12);
    title(ax, 'Chaotic Motion', 'fontsize', 14);
    F = getframe;
end

movie(F, fram, 20);

end
