clear;clc;close all;

% Laplace variable
s = tf('s');


%% EXAMPLE 1
% 
H = (s^2+0.1*s+7.5)/(s^4+0.12*s^3+9*s^2);
figure;bode(H);
figure;bode2(H);



%% EXAMPLE 2

H = 1/(s*(s+0.1));
figure;bode(H);
figure;bode2(H);



%% EXAMPLE 3

H = tf([-0.1,-2.4,-181,-1950],[1,3.3,990,2600]);
figure;bode(H);
figure;bode2(H);