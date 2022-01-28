function my_double_pendulum
% Calls double_pendulum with particular variables.
%
% This is the main entry function to call. This is a source file where you
% can twiddle the variables for initial conditions and whatnot, without
% modifying the core logic.

% TODO: What do all these short variable names mean?
l1 = 1;
l2 = 2;
m1 = 2;
m2 = 1;

% Initial condition

theta1_0 = 1.6;
theta1_prime_0 = 0;
theta2_0 = 2.2;
theta2_prime_0 = 0;

% How much time, I think?

tspan = 50;

double_pendulum(l1, l2, m1, m2, tspan, theta1_0, theta1_prime_0, ...
    theta2_0, theta2_prime_0);

end