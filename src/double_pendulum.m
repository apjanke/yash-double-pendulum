function double_pendulum(l1, l2, m1, m2, tspan, theta1_0, theta1_prime_0, ...
    theta2_0, theta2_prime_0)
% Make some plots of double pendulum behavior.
%
% This function takes a bunch of inputs. See the my_double_pendulum
% function for an example of how to set them up.
%
% See also:
% MY_DOUBLE_PENDULUM

%#ok<*DEFNU>

% Gravitational constant.
g = 9.8;


% Here's the step function that calculates the pendulum behavior.

    function yprime = my_local_pend(t, y) %#ok<INUSL>
        % Pendulum action step function
        
        % NOTE: I do not know why we are ignoring the "t" input argument!
        yprime = zeros(4, 1);
        
        a = (m1 + m2) * l1;
        b = m2 * l2 * cos(y(1) - y(3));
        c = m2 * l1 * cos(y(1) - y(3));
        d = m2 * l2;
        e = -m2 * l2 * y(4) * y(4) * sin(y(1) - y(3)) - g * (m1 + m2) * sin(y(1));
        f = m2 * l1 * y(2) * y(2) * sin(y(1) - y(3)) - m2 * g * sin(y(3));
        
        yprime(1) = y(2);
        yprime(3) = y(4);
        yprime(2) = (e*d-b*f) / (a*d-c*b);
        yprime(4) = (a*f-c*e) / (a*d-c*b);
        
    end

y_0 = [theta1_0 theta1_prime_0 theta2_0 theta2_prime_0];
[t, y] = ode45(@my_local_pend, [0, tspan], y_0);

% Here's the time series of each variable
[theta1, theta1_prime, theta2, theta2_prime] = deal(y(:,1), y(:,2), y(:,3), y(:,4));

% Position of mass 1 and mass 2

x1 = l1 * sin(y(:,1));
y1 = -l1 * cos(y(:,1));
x2 = l1 * sin(y(:,1)) + l2 * sin(y(:,3));
y2 = -l1 * cos(y(:,1)) - l2 * cos(y(:,3));

% Graphics and plotting

fig = figure;
set(fig, 'color', 'white');

ax = subplot(2, 2, 1);
plot_xy_procession_path(ax, x1, y1, x2, y2);
ax = subplot(2, 2, 2);
plot_thetas_vs_time(ax, t, theta1, theta2, theta1_0, theta2_0);
ax = subplot(2, 2, 3);
plot_theta_primes_vs_time(ax, t, theta1_prime, theta2_prime);

% Movie of double pendulum
%
% Actually we don't care about the movie; turn this off.
% plot_movie(x1, y1, x2, y2, y, l1, l2);

end

function plot_xy_procession_path(ax, x1, y1, x2, y2)
plot(ax, x1, y1, 'linewidth', 2);
hold on
plot(ax, x2, y2, 'r', 'linewidth', 2);
hold off
xlabel(ax, 'X', 'fontSize', 14);
ylabel(ax, 'Y', 'fontSize', 14);
title(ax, 'Chaotic Double Pendulum', 'fontsize', 14);
end

function plot_thetas_vs_time(ax, t, theta1, theta2, theta1_0, theta2_0)
plot(ax, t, theta1, 'linewidth', 2);
hold on
plot(ax, t, theta2, 'r', 'linewidth', 2);
hold off
legend(ax, '\theta_1', '\theta_2');
xlabel(ax, 'time', 'fontSize', 14);
ylabel(ax, 'theta', 'fontSize', 14);
title(ax, ...
  sprintf('\\theta_1(t=0)=%.1f and \\theta_2(t=0)=%.1f', theta1_0, theta2_0), ...
  'fontsize', 14);
end

function plot_theta_primes_vs_time(ax, t, theta1_prime, theta2_prime)
plot(ax, t, [theta1_prime, theta2_prime]);
title(ax, '\theta1\_prime and \theta2\_prime over time');
legend(ax, 'theta1 prime', 'theta2 prime');
xlabel(ax, 'Time');
ylabel(ax, 'theta primes');
end

function plot_movie(x1, y1, x2, y2, y, l1, l2)
    fig = figure;
    nFrames = 0;
    fram = 0;
    ax = gca;
    apply_plot_cosmetics(fig, ax);
    
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
        xlabel(ax, 'X', 'fontSize', 12);
        ylabel(ax, 'Y', 'fontSize', 12);
        title(ax, 'Chaotic Motion', 'fontsize', 14);
        F = getframe;
    end
    
    movie(F, fram, 20);

end


