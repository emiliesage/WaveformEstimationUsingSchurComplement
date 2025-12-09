close all; clear; clc;
%Chat Gpt was used to find an optimization method to reduce error and get rid of negative coupling

%% ============================================================
% 1) Monte-Carlo training
%% ============================================================

% Model sampling (used for training signals)
f_sample_train = 2.5e6;    % training sample rate (Hz)
f_signal_nom   = 125e3;    % nominal frequency used for LS basis in training
n_samples      = 1024;
tn_train       = (0:n_samples-1) / f_sample_train;

n_simulations = 5e3; % increase for stable covariances

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
fprintf('Running Monte-Carlo training (%d simulations)...\n', n_simulations);
for ii = 1:n_simulations
    if mod(ii, 5000) == 0
        fprintf('  Progress: %d/%d\n', ii, n_simulations);
    end
    
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

    % Estimate frequency from phase progression across the window
    % Use instantaneous frequency method: measure phase change over time
    window_size = 256;
    hop_size = 128;
    num_windows = floor((n_samples - window_size) / hop_size) + 1;
    
    if num_windows >= 2
        phases = zeros(num_windows, 1);
        times = zeros(num_windows, 1);
        
        for ww = 1:num_windows
            start_idx = (ww-1)*hop_size + 1;
            end_idx = start_idx + window_size - 1;
            
            seg = sig(start_idx:end_idx);
            t_seg = tn_train(start_idx:end_idx);
            
            % LS on this segment using nominal frequency
            E_seg = [sin(2*pi*f_signal_nom*t_seg(:)), cos(2*pi*f_signal_nom*t_seg(:)), ones(window_size,1)];
            LS_seg = (E_seg' * E_seg) \ E_seg';
            p_seg = LS_seg * seg(:);
            
            phases(ww) = atan2(p_seg(2), p_seg(1));
            times(ww) = mean(t_seg);
        end
        
        % Unwrap phases and estimate frequency from phase slope
        phases_unwrap = unwrap(phases);
        
        % Linear fit to get phase rate
        P = polyfit(times, phases_unwrap, 1);
        phase_rate = P(1); % rad/s
        
        % Frequency estimate: f = (f_nominal + phase_rate/(2*pi))
        f_est(ii) = f_signal_nom + phase_rate / (2*pi);
    else
        % Fallback if not enough windows
        f_est(ii) = f_signal_nom;
    end
end

fprintf('Monte-Carlo training complete.\n');

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

fprintf('Training statistics computed. Schur complement ready.\n\n');

%% ============================================================
% 2) Real-time loop (reads files, computes f_sample per file)
%% ============================================================

directory_samples = '/home/emilie/WaveformEstimationUsingSchurComplement/samples';
opts = [];

fprintf("Ready. Watching folder: %s\n", directory_samples);

% plotting setup
figure('Name','Live Waveform Processing','NumberTitle','off','Color','w');
ax = axes; hold on;
h_raw = plot(nan,nan,'r.','MarkerSize',6);
h_ls  = plot(nan,nan,'g-','LineWidth',1.2);
h_sc  = plot(nan,nan,'b-','LineWidth',1.2);
ylim([-1.5 1.5]);
xlabel("Time [s]"); ylabel("Voltage [V]");
title("Raw (red), LS (green), Schur (blue)");
legend('raw','LS','Schur','Location', 'southeast');
h_freq_text = text(0.02, 0.95, '', 'Units','normalized','FontSize',10,'Color','k','FontWeight','bold');

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

    % Estimate frequency from phase progression (same method as training)
    window_size = min(256, floor(N/4));
    hop_size = max(1, floor(window_size/2));
    num_windows = floor((N - window_size) / hop_size) + 1;
    
    if num_windows >= 2
        phases = zeros(num_windows, 1);
        times = zeros(num_windows, 1);
        
        for ww = 1:num_windows
            start_idx = (ww-1)*hop_size + 1;
            end_idx = start_idx + window_size - 1;
            
            seg = voltage(start_idx:end_idx);
            t_seg = time(start_idx:end_idx);
            
            % LS on this segment using nominal frequency
            E_seg = [sin(2*pi*f_signal_nom*t_seg), cos(2*pi*f_signal_nom*t_seg), ones(window_size,1)];
            LS_seg = (E_seg' * E_seg) \ E_seg';
            p_seg = LS_seg * seg;
            
            phases(ww) = atan2(p_seg(2), p_seg(1));
            times(ww) = mean(t_seg);
        end
        
        % Unwrap phases and estimate frequency from phase slope
        phases_unwrap = unwrap(phases);
        
        % Linear fit to get phase rate
        P = polyfit(times, phases_unwrap, 1);
        phase_rate = P(1); % rad/s
        
        % Frequency estimate
        f_est_phase = f_signal_nom + phase_rate / (2*pi);
    else
        % Fallback if not enough windows
        f_est_phase = f_signal_nom;
    end

    % Build measured vector y_meas EXACTLY as training y:
    % y = [A_est, sin(phi_est), cos(phi_est), c_est, f_est]
    y_meas = [A_LS, sin(phi_LS), cos(phi_LS), c_LS, f_est_phase];

    % Schur reconstruction (affine map)
    rec = muX' + Sigma_xy * (Sigma_yy_reg \ (y_meas' - muY'));

    % Recover parameters from rec (note ordering follows x definition)
    A_SC = rec(1);
    sinphi_SC = rec(2);
    cosphi_SC = rec(3);
    c_SC = rec(4);
    f_SC_initial = rec(5);

    % Recover phase
    phi_SC = atan2(sinphi_SC, cosphi_SC);

    %% ========================================
    %  OPTIMIZE Schur parameters using grid search
    %% ========================================
    
    % Search around the Schur-predicted frequency
    f_search_range = 2e3;  % +/- 2 kHz around Schur estimate
    f_search_step = 100;   % 100 Hz steps
    
    f_grid = max(f_SC_initial - f_search_range, 50e3) : f_search_step : (f_SC_initial + f_search_range);
    mse_grid = zeros(size(f_grid));
    
    for jj = 1:length(f_grid)
        f_test = f_grid(jj);
        sig_test = A_SC * sin(2*pi*f_test*time + phi_SC) + c_SC;
        mse_grid(jj) = mean((voltage - sig_test).^2);
    end
    
    [~, idx_min] = min(mse_grid);
    f_SC = f_grid(idx_min);
    
    % Refine with local quadratic fit if not at boundary
    if idx_min > 1 && idx_min < length(f_grid)
        f1 = f_grid(idx_min-1); mse1 = mse_grid(idx_min-1);
        f2 = f_grid(idx_min);   mse2 = mse_grid(idx_min);
        f3 = f_grid(idx_min+1); mse3 = mse_grid(idx_min+1);
        
        % Quadratic interpolation for sub-grid accuracy
        denom = (f1-f2)*(f1-f3)*(f2-f3);
        if abs(denom) > 1e-10
            A_q = (f3*(mse2-mse1) + f2*(mse1-mse3) + f1*(mse3-mse2)) / denom;
            B_q = (f3^2*(mse1-mse2) + f2^2*(mse3-mse1) + f1^2*(mse2-mse3)) / denom;
            if abs(A_q) > 1e-10
                f_refined = -B_q / (2*A_q);
                if f_refined >= f1 && f_refined <= f3
                    f_SC = f_refined;
                end
            end
        end
    end

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
    
    mse_Schur = mean((voltage - sig_SC).^2);
    mse_LS = mean((voltage - sig_LS).^2);
    
    set(h_freq_text, 'String', sprintf('f_{Schur} = %.3f kHz | f_{LS} = %.3f kHz | MSE_{Schur} = %.5f | MSE_{LS} = %.5f', ...
        f_SC/1e3, f_signal_nom/1e3, mse_Schur, mse_LS));

    xlim([min(time), max(time)]);
    drawnow;

    % optional small pause
    pause(0.05);
end
