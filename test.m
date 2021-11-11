clear;clc;close all;


% Laplace variable
s = tf('s');

% system transfer function
G = 1/(s*(s+0.1));

opts = bodeoptions('cstprefs');
opts.Magunits = 'abs';
opts.Magscale = 'log';

figure;
margin(G,opts);
grid on;


figure;
bode(G,opts);
grid on;

figure;
[bode_plot,mag_plot,phase_plot] = bode2(G);
mag_plot
phase_plot
bode_plot