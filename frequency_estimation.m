function frequency_estimation()
% FREQUENCY_ESTIMATION  Estimate speech formants using AR spectral analysis.
%
%   Supports file input or live recording. Applies a pre-emphasis filter,
%   estimates AR models using three methods (Yule-Walker, Levinson-Durbin,
%   and Burg), plots the resulting spectra, and marks up to seven formants
%   on each plot.
%
%   Requirements: MATLAB Signal Processing Toolbox
%
%   Usage:
%       frequency_estimation()
%
%   The function will prompt you to choose between loading an audio file
%   and making a live recording.  It then produces a three-panel figure
%   showing the AR spectrum obtained by each method with formant
%   frequencies annotated.

    fprintf('Frequency Estimation in Audio Signal using MATLAB\n');
    fprintf('--------------------------------------------------\n');
    choice = input('Select input source (1 = audio file, 2 = live recording): ');

    if choice == 1
        %% File input
        [fname, fpath] = uigetfile( ...
            {'*.wav;*.mp3;*.flac;*.ogg', 'Audio Files (*.wav, *.mp3, *.flac, *.ogg)'; ...
             '*.*',                       'All Files (*.*)'}, ...
            'Select an audio file');

        if isequal(fname, 0)
            error('frequency_estimation:noFile', 'No file selected.');
        end

        [y, fs] = audioread(fullfile(fpath, fname));
        fprintf('Loaded "%s"  (%.2f s, %d Hz, %d ch)\n', ...
                fname, size(y,1)/fs, fs, size(y,2));

    elseif choice == 2
        %% Live recording
        fs       = 44100;
        nBits    = 16;
        nChans   = 1;
        duration = input('Recording duration in seconds: ');

        if duration <= 0
            error('frequency_estimation:badDuration', ...
                  'Duration must be a positive number.');
        end

        fprintf('Recording for %.1f second(s)... ', duration);
        recObj = audiorecorder(fs, nBits, nChans);
        recordblocking(recObj, duration);
        fprintf('Done.\n');
        y = getaudiodata(recObj);

    else
        error('frequency_estimation:badChoice', ...
              'Invalid choice. Enter 1 (file) or 2 (recording).');
    end

    %% Pre-processing -------------------------------------------------------

    % Mix down to mono
    if size(y, 2) > 1
        y = mean(y, 2);
    end
    y = y(:);               % ensure column vector

    % Normalize to [-1, 1]
    peak = max(abs(y));
    if peak > 0
        y = y / peak;
    end

    % Pre-emphasis filter: H(z) = 1 - 0.97 z^{-1}
    % Boosts high-frequency content and compensates for spectral tilt.
    preEmph = 0.97;
    y_pe = filter([1, -preEmph], 1, y);

    %% AR model parameters --------------------------------------------------
    modelOrder  = 16;    % AR model order (standard for speech analysis)
    nFft        = 2048;  % FFT points for spectral display
    maxFormants = 7;     % Maximum number of formants to detect

    %% Estimate AR coefficients with three methods --------------------------

    % 1. Yule-Walker (uses the biased autocorrelation estimator)
    [a_yw, ~] = aryule(y_pe, modelOrder);

    % 2. Levinson-Durbin (solved from biased autocorrelation directly)
    r          = xcorr(y_pe, modelOrder, 'biased');   % lags -P … P
    r_onesided = r(modelOrder + 1 : end);             % lags  0 … P
    [a_ld, ~]  = levinson(r_onesided, modelOrder);

    % 3. Burg (minimises forward and backward prediction errors)
    [a_bg, ~] = arburg(y_pe, modelOrder);

    methodNames  = {'Yule-Walker', 'Levinson-Durbin', 'Burg'};
    coefficients = {a_yw, a_ld, a_bg};

    %% Plot -----------------------------------------------------------------
    figure('Name',        'Frequency Estimation using AR Methods', ...
           'NumberTitle', 'off', ...
           'Position',    [100, 100, 960, 720]);

    colors = lines(3);   % distinct colours for each subplot

    for i = 1:3
        a = coefficients{i};

        % Frequency response of the all-pole AR model:  H(z) = 1 / A(z)
        [H, f] = freqz(1, a, nFft, fs);
        mag_dB = 20 * log10(abs(H) + eps);
        mag_dB = mag_dB - max(mag_dB);   % normalise peak to 0 dB

        % Detect formants
        [fFreqs, ~] = findFormants(abs(H), f, maxFormants);

        % Draw spectrum
        subplot(3, 1, i);
        plot(f, mag_dB, 'Color', colors(i, :), 'LineWidth', 1.4);
        hold on;

        % Annotate each detected formant
        for k = 1 : length(fFreqs)
            [~, idx]  = min(abs(f - fFreqs(k)));
            amp_at_f  = mag_dB(idx);
            marker_c  = [0.85, 0.1, 0.1];   % dark red

            plot(fFreqs(k), amp_at_f, 'v', ...
                 'MarkerSize',      8, ...
                 'MarkerFaceColor', marker_c, ...
                 'MarkerEdgeColor', marker_c * 0.6);

            text(fFreqs(k), amp_at_f + 2.5, ...
                 sprintf('F%d\n%.0f Hz', k, fFreqs(k)), ...
                 'HorizontalAlignment', 'center', ...
                 'FontSize',            7, ...
                 'Color',               marker_c);
        end

        hold off;
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        title(sprintf('%s Method', methodNames{i}));
        xlim([0, min(fs/2, 8000)]);
        ylim([-60, 5]);
        grid on;
    end

    sgtitle('Formant Detection using AR Spectral Methods', 'FontSize', 13);
end


% =========================================================================
function [formantFreqs, formantAmps] = findFormants(Hmag, freq, maxN)
% FINDFORMANTS  Locate formant peaks in an AR magnitude spectrum.
%
%   [formantFreqs, formantAmps] = findFormants(Hmag, freq, maxN)
%
%   Inputs:
%       Hmag  - Linear magnitude of the AR frequency response (column vector)
%       freq  - Corresponding frequency axis in Hz (column vector)
%       maxN  - Maximum number of formants to return
%
%   Outputs:
%       formantFreqs - Formant frequencies in Hz, sorted ascending
%       formantAmps  - Corresponding linear magnitude values
%
%   Only peaks in the typical speech formant range [80, 4500] Hz are
%   considered.

    fMin = 80;     % Hz  – lower bound for formant search
    fMax = 4500;   % Hz  – upper bound  (F1-F4 nearly always below 4 kHz)

    % Zero out components outside the region of interest so findpeaks
    % does not select them.
    H_roi          = Hmag;
    H_roi(freq < fMin | freq > fMax) = 0;

    minProminence = 0.05 * max(H_roi);   % at least 5 % of the ROI peak

    [~, locs] = findpeaks(H_roi, ...
                          'MinPeakProminence', minProminence, ...
                          'SortStr',           'descend', ...
                          'NPeaks',            maxN);

    % Sort selected peaks by ascending frequency
    locs         = sort(locs);
    formantFreqs = freq(locs);
    formantAmps  = Hmag(locs);
end
