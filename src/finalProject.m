close all; clear all; clc;

%%Creation of The sinusoidal with random noise added
f_sample = 2.5e6;
f_signal = 125e3;
samples = 1024;
w0 = 2 * pi * f_signal;
A = 1;
phi = pi/4;
k = 0.5;
t = (0:samples-1) / f_sample;
signal = A * sin(w0 * t + phi) + k;
signal = signal.';

nois_mean = 0.0;
noise_std = 0.01

noise = 0.5*(rand(samples,1)-0.5);

signal = signal + noise;
figure(); hold on;
subplot(2,1,1);
plot(t, signal,'r.');

subplot(2,1,2);hold on;
plot(t, signal,'r.');
xlim( [0 (samples/f_sample)/16] );

%%signal reconstruction
E = [cos(w0 * t).', sin(w0 * t).', ones(samples,1)];

p = inv(E.'*E) * E.' * signal;

alpha = p(1);
beta = p(2);
k_prime = p(3);

A_prime = sqrt(alpha^2 + beta^2);
phi_prime = atan(alpha/beta);

reconstructed_signal = A_prime * sin(w0*t + phi_prime) + k_prime;

subplot(2,1,1);hold on;
plot(t, reconstructed_signal,'b.');
title_formated = sprintf("A_real %0.3e A_estimated %0.3e Error %0.3e",A,A_prime, A-A_prime);
title(title_formated);
subplot(2,1,2);hold on;
plot(t, reconstructed_signal,'b');
xlim( [0 (samples/f_sample)/16] )


