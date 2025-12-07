
close all; clear; clc;

%% ============================================================
% 1) Monte-Carlo training
%% ============================================================

% Model sampling (used for training signals)
f_sample_train = 2.5e6;    % training sample rate (Hz)
f_signal_nom   = 125e3;    % nominal frequency used for LS basis in training
n_samples      = 1024;
tn_train       = (0:n_samples-1) / f_sample_train;

n_simulations = 5e4; % increase for stable covariances

% Parameter distributions (recommended)
A_signal_mu = 1.0;    A_signal_std = 0.1;
phi_mu      = 0.0;    phi_std      = pi/12;
c_mu        = 0.0;    c_std        = 0.05;
fsignal_mu  = 125e3;  fsignal_std  = 0.05 * 125e3; % tune to real-world drift
noise_mu    = 0.0;    noise_std    = 0.05;

% Preallocate storage
A_ref = zeros(n_simulations,1);
phi_ref = zeros(n_simulations,1);
c_ref = zeros(n_simulations,1);
fsignal_ref = zeros(n_simulations,1);

A_est = zeros(n_simulations,1);
phi_est = zeros(n_simulations,1);
c_est = zeros(n_simulations,1);
f_est = zeros(n_simulations,1);

% Precompute LS operator for the assumed LS basis (fixed nominal f_signal_nom)
E_train = [sin(2*pi*f_signal_nom*tn_train(:)), cos(2*pi*f_signal_nom*tn_train(:)), ones(n_samples,1)];
LS_train = (E_train' * E_train) \ E_train';

% Monte-Carlo loop
for ii = 1:n_simulations
    % draw parameters
    A_ref(ii) = A_signal_mu + A_signal_std * randn;
    phi_ref(ii) = phi_mu + phi_std * randn;
    c_ref(ii) = c_mu + c_std * randn;
    fsignal_ref(ii) = fsignal_mu + fsignal_std * randn;

    % generate noisy signal
    w = 2*pi*fsignal_ref(ii);
    sig = A_ref(ii) * sin(w*tn_train + phi_ref(ii)) + c_ref(ii);
    sig = sig + noise_mu + noise_std * randn(1, n_samples);

    % LS estimate using assumed basis (same approach as online LS)
    p = LS_train * sig.';      % 3x1
    A_est(ii) = hypot(p(1), p(2));
    phi_est(ii) = atan2(p(2), p(1));
    c_est(ii) = p(3);

    % Simulate the online FFT-based frequency estimator on the training signal
    % use zero-padding to improve frequency resolution (same method online)
    Nfft = 4 * n_samples;
    win = hann(n_samples);
    Y = fft((sig(:) .* win), Nfft);
    mag = abs(Y(1:floor(Nfft/2)));
    [~, idx] = max(mag);
    % quadratic interpolation
    if idx>1 && idx < floor(Nfft/2)
        alpha = mag(idx-1); beta = mag(idx); gamma = mag(idx+1);
        p_quadr = 0.5*(alpha - gamma) / (alpha - 2*beta + gamma);
    else
        p_quadr = 0;
    end
    f_est(ii) = ((idx-1) + p_quadr) * (f_sample_train / Nfft);
end

% Build training x (true) and y (measured) using sin/cos for phase and include f_est
% x = [A_true, sin(phi_true), cos(phi_true), c_true, f_true]
% y = [A_est, sin(phi_est), cos(phi_est), c_est, f_est]
x = [A_ref, sin(phi_ref), cos(phi_ref), c_ref, fsignal_ref];
y = [A_est, sin(phi_est), cos(phi_est), c_est, f_est];

% Means and covariance (stacked)
muX = mean(x, 1);   % 1 x nx
muY = mean(y, 1);   % 1 x ny

z = [x'; y'];
gam = cov(z');      % covariance of stacked vector

nx = size(x,2);     % 5
ny = size(y,2);     % 5

Sigma_xx = gam(1:nx, 1:nx);
Sigma_xy = gam(1:nx, nx+1 : nx+ny);
Sigma_yy = gam(nx+1 : nx+ny, nx+1 : nx+ny);

% Regularize Sigma_yy and check conditioning
reg = 1e-6;
Sigma_yy_reg = Sigma_yy + reg * eye(ny);
if cond(Sigma_yy_reg) > 1e8
    reg = 1e-3;
    Sigma_yy_reg = Sigma_yy + reg * eye(ny);
    warning('Increased Sigma_yy regularization to %g due to poor conditioning.', reg);
end

%% ============================================================
% 2) Real-time loop (reads files, computes f_sample per file)
%% ============================================================

directory_samples = '/home/emilie/WaveformEstimationUsingSchurComplement/samples';
opts = [];

fprintf("Ready. Watching folder...\n");

% plotting setup
figure('Name','Live Waveform Processing','NumberTitle','off','Color','w');
ax = axes; hold on;
h_raw = plot(nan,nan,'r.','MarkerSize',6);
h_ls  = plot(nan,nan,'g-','LineWidth',1.2);
h_sc  = plot(nan,nan,'b-','LineWidth',1.2);
ylim([-1.5 1.5]);
xlabel("Time [s]"); ylabel("Voltage [V]");
title("Raw (red), LS (green), Schur (blue)");
legend({'raw','LS','Schur'});
h_freq_text = text(0.02, 0.95, '', 'Units','normalized','FontSize',12,'Color','k','FontWeight','bold');

% keep small history of frequency estimates
max_hist = 200;
freq_hist = nan(1, max_hist);
hist_idx = 0;

while true
    listing = dir(directory_samples);
    if numel(listing) < 4
        pause(0.25);
        fprintf("Waiting for new data...\n");
        continue;
    end

    fprintf("Collecting new file...\n");
    name_file = listing(end-1).name;
    fullpath = fullfile(directory_samples, name_file);

    if isempty(opts)
        opts = detectImportOptions(fullpath);
    end

    M = readmatrix(fullpath, opts);

    % clean directory: preserve ".", "..", last file
    for k = 3:(numel(listing)-1)
        delete(fullfile(directory_samples, listing(k).name));
    end

    % extract data
    time = M(1:end-1,1);
    voltage = M(1:end-1,2);

    % compute actual sampling frequency from timestamps (critical)
    dt_vec = diff(time);
    dt = mean(dt_vec);
    if any(dt_vec <= 0)
        error('Non-monotonic or zero timestamp differences in file.');
    end
    f_sample_live = 1 / dt;

    N = length(time);

    % LS estimate using the SAME assumed LS basis used in training (f_signal_nom)
    % (this keeps A_est & phi_est consistent with training)
    E_live = [sin(2*pi*f_signal_nom*time), cos(2*pi*f_signal_nom*time), ones(N,1)];
    LS_live = (E_live' * E_live) \ E_live';
    p = LS_live * voltage;

    A_LS   = hypot(p(1), p(2));
    phi_LS = atan2(p(2), p(1));
    c_LS   = p(3);

    % Online FFT frequency estimate (use f_sample_live and zero-padding)
    Nfft = 4 * N;
    win = hann(N);
    Y = fft((voltage(:) .* win), Nfft);
    mag = abs(Y(1:floor(Nfft/2)));
    [~, idx] = max(mag);
    if idx>1 && idx < floor(Nfft/2)
        alpha = mag(idx-1); beta = mag(idx); gamma = mag(idx+1);
        p_quadr = 0.5*(alpha - gamma) / (alpha - 2*beta + gamma);
    else
        p_quadr = 0;
    end
    f_est_fft = ((idx-1) + p_quadr) * (f_sample_live / Nfft);

    % Build measured vector y_meas EXACTLY as training y:
    % y = [A_est, sin(phi_est), cos(phi_est), c_est, f_est]
    y_meas = [A_LS, sin(phi_LS), cos(phi_LS), c_LS, f_est_fft];

    % Schur reconstruction (affine map)
    rec = muX' + Sigma_xy * (Sigma_yy_reg \ (y_meas' - muY'));

    % Recover parameters from rec (note ordering follows x definition)
    A_SC = rec(1);
    sinphi_SC = rec(2);
    cosphi_SC = rec(3);
    c_SC = rec(4);
    f_SC = rec(5);

    % Recover phase
    phi_SC = atan2(sinphi_SC, cosphi_SC);

    % Build signals for plotting
    sig_LS = A_LS * sin(2*pi*f_signal_nom*time + phi_LS) + c_LS;   % LS uses nominal basis (training-consistent)
    sig_SC = A_SC * sin(2*pi*f_SC*time + phi_SC) + c_SC;

    % Update plots
    set(h_raw, 'XData', time, 'YData', voltage);
    set(h_ls,  'XData', time, 'YData', sig_LS);
    set(h_sc,  'XData', time, 'YData', sig_SC);

    % Update frequency text and history
    hist_idx = hist_idx + 1;
    if hist_idx > max_hist, hist_idx = 1; end
    freq_hist(hist_idx) = f_SC / 1e3; % store in kHz

    set(h_freq_text, 'String', sprintf('f_{Schur} = %.3f kHz | f_{FFT} = %.3f kHz | f_s = %.3f MHz', ...
        f_SC/1e3, f_est_fft/1e3, f_sample_live/1e6));

    drawnow;

    % optional small pause
    pause(0.05);
end
