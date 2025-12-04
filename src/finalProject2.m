close all; clear all; clc;

%%Creation of The sinusoidal with random noise added
f_sample = 2.5e6;
f_signal = 125e3;
n_samples = 1024;
w0 = 2 * pi * f_signal;
##A = 1;
##phi = pi/4;
##k = 0.5;
tn = (0:n_samples-1) / f_sample;


n_simulations = 1e3;

% noise
%noise_mu = 0.0 ;
%noise_std = 0.1 ;

A_signal_mu = 1.0;   A_signal_std = 0.1;
phi_mu = 0.0;        phi_std = pi/4;
c_mu = 0.0;          c_std = 0.1;
fsignal_mu = 125e3;  fsignal_std = 0.005 * 125e3; % 125 Hz
noise_mu = 0.0;      noise_std = 0.1;

A_ref = zeros(n_simulations,1);
fsignal_ref =  zeros(n_simulations,1);
phi_ref = zeros(n_simulations,1);
c_ref = zeros(n_simulations,1);

A_estimation = zeros(n_simulations,1);
phi_estimation = zeros(n_simulations,1);
c_estimation = zeros(n_simulations,1);
fsignal_estimation =  zeros(n_simulations,1);

% distribution frequency
%fsignal_mu = 125e3;
%fsignal_std = 0.005 * 125e3;
E = [ (sin(w0*tn)).'  (cos(w0*tn)).'  ones(n_samples,1) ] ;
for ii=1:n_simulations

    A_ref(ii,1) = A_signal_mu + A_signal_std * randn(1,1);
    fsignal_ref(ii,1) = fsignal_mu + fsignal_std * randn(1,1);
    phi_ref(ii,1) = phi_mu + phi_std * randn(1,1);
    c_ref(ii,1) = c_mu + c_std * randn(1,1);

    fsignal_estimation(ii,1) = f_signal;
    w0_ref = 2*pi*fsignal_ref(ii,1) ;
    signal = A_ref(ii,1)*sin(w0_ref*tn + phi_ref(ii,1)) + c_ref(ii,1) ;

    % add noise
    signal = signal + ( noise_mu +  noise_std * randn(1,n_samples) );

    % estimation
    p = ( inv(E'*E)*E' ) * signal.' ;
    A_estimation(ii,1) = sqrt( p(1)^2 + p(2)^2  );
    phi_estimation(ii,1) = atan2( p(2) , p(1) );
    c_estimation(ii,1) = p(3);


end


figure;
A_error = A_estimation-A_ref;
plot(A_error)

%A_error_mean = mean(A_error);
%A_error_std = std(A_error);
%str_aux = sprintf('mean %e and std %e \n',A_error_mean,A_error_std);
%title(str_aux)

x = [A_ref, phi_ref, c_ref, fsignal_ref];
y = [A_estimation, phi_estimation, c_estimation, fsignal_estimation];

z = [x'; y'];
muX = mean(x);
muY = mean(y);
gam = cov(z');

n_params = 4;
Sigma_xx = gam(1:n_params, 1:n_params);
Sigma_xy = gam(1:n_params, n_params+1:end);
Sigma_yx = gam(n_params+1:end, 1:n_params);
Sigma_yy = gam(n_params+1:end, n_params+1:end);
Sigma_yy_reg = Sigma_yy + 1e-9 * eye(size(Sigma_yy));

%A = 1;
%phi = pi/4;
%c = 0.5;
%f = 125e3;
%w = 2 * pi * f;
%

A = A_signal_mu + A_signal_std * randn(1,1);
f = fsignal_mu + fsignal_std * randn(1,1);
phi = phi_mu + phi_std * randn(1,1);
c = c_mu + c_std * randn(1,1);

w = 2*pi*f ;
signal = A*sin(w*tn + phi) + c ;

    % add noise
signal = signal + ( noise_mu +  noise_std * randn(1,n_samples) );

p = ( inv(E'*E)*E' ) * signal.' ;
A_2 = sqrt( p(1)^2 + p(2)^2  );
phi_2 = atan2( p(2) , p(1) );
c_2 = p(3);
y1 = [A_2,phi_2,c_2,125e3]
rec = muX' + Sigma_xy * inv(Sigma_yy_reg) * (y1' - muY')

sig_one = A_2 * sin(2*pi*125e3*tn + phi_2) + c_2;
reconstructed_signal = rec(1)*sin(2*pi*rec(4)*tn + rec(2)) + rec(3);

figure;
hold on;
plot(tn,signal,'r.');
plot(tn, reconstructed_signal,'b');
plot(tn,sig_one,'g')
xlim( [0 (n_samples/f_sample)/16] );
title("Signal Estimation 1%");
xlabel("Time");
ylabel("Amplitude");
hold off;
legend("Noisy Signal","Schur Complement Estimation", "Least Squares Estimation");

