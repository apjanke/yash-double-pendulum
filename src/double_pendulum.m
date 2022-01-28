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
theta1_0 = 1.6;
theta1_prime_0 = 0;
theta2_0 = 2.2;
theta2_prime_0 = 0;

y_0 = [theta1_0 theta1_prime_0 theta2_0 theta2_prime_0];
[t, y] = ode45(@pend, [0, tspan], [2.5 0 1 0]);

% Here's the time series of each variable
[theta1, theta1_prime, theta2, theta2_prime] = deal(y(:,1), y(:,2), y(:,3), y(:,4));

% Position of mass 1 and mass 2

x1 = l1 * sin(y(:,1));
y1 = -l1 * cos(y(:,1));
x2 = l1 * sin(y(:,1)) + l2 * sin(y(:,3));
y2 = -l1 * cos(y(:,1)) - l2 * cos(y(:,3));

% Visualizing the result

fig1 = figure;
ax = gca;
set(fig1, 'color', 'white');
plot(x1, y1, 'linewidth', 2);
hold on
plot(x2, y2, 'r', 'linewidth', 2);
set(ax, 'fontSize', 14);
xlabel('X', 'fontSize', 14);
ylabel('Y', 'fontSize', 14);
title('Chaotic Double Pendulum','fontsize', 14);

fig2 = figure;
ax = gca;
set(fig2, 'color', 'white');
plot(y(:,1), 'linewidth', 2);
hold on
plot(y(:,3), 'r', 'linewidth', 2);
set(ax, 'fontSize', 14);
legend('\theta_1', '\theta_2');
xlabel('time', 'fontSize', 14);
ylabel('theta', 'fontSize', 14);
title('\theta_1(t=0)=2.5 and \theta_2(t=0)=1.0','fontsize', 14);

% Theta1 and 2 prime over time

fig4 = figure;
set(fig4, 'color', 'white');
ax = gca;
plot(t, [theta1_prime, theta2_prime]);
title(ax, 'theta1\_prime and theta2\_prime over time');
xlabel(ax, 'Time');
ylabel(ax, 'Theta 1 and 2 prime');


% Movie of double pendulum

% Actually we don't care about the movie; turn this off

if false
    
    fig3 = figure; %#ok<UNRCH>
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

end
