% plot_ber_curves.m
% Script to plot BER curves from simulation results

clear;
close all;

results_dir = '../res';
files = dir(fullfile(results_dir, 'results_*.txt'));

figure;
hold on;
grid on;

markers = {'-o', '-s', '-^', '-x', '-d'};
colors = lines(numel(files));

for k = 1:length(files)
    filename = files(k).name;
    filepath = fullfile(results_dir, filename);

    % Load data
    data = load(filepath);
    snr = data(:,1);
    ber = data(:,2);

    % Determine coding type from filename
    if contains(filename, 'ONLY_CC')
        label = 'Convolutional Code';
    elseif contains(filename, 'ONLY_RS')
        label = 'Reed-Solomon Code';
    elseif contains(filename, 'RS_AND_CC')
        label = 'RS + Convolutional Code';
    else
        label = 'Unknown';
    end

    % Extract code rate from filename, e.g., '_2_3' => 1/2
    rate_match = regexp(filename, '_([0-9]+)_([0-9]+)', 'tokens');
    if ~isempty(rate_match)
        rate_num = rate_match{1}{2}; % denominator in the filename
        rate_den = rate_match{1}{1}; % numerator in the filename
        code_rate = [rate_den '/' rate_num];
    else
        code_rate = 'Unknown';
    end

    % Extract additional parameters for more detailed labeling
    intlv_match = regexp(filename, 'intlv(\d+)', 'tokens');
    if ~isempty(intlv_match)
        intlv = intlv_match{1}{1};
    else
        intlv = 'N/A';
    end
    basis = 'Dual Basis';
    if contains(filename, 'noDualBasis')
        basis = 'No Dual Basis';
    end

    % Plot
    semilogy(snr, ber, markers{mod(k-1,length(markers))+1}, 'LineWidth', 1.5, 'Color', colors(k,:),...
        'DisplayName', [label ', Rate=' code_rate ', intlv=' intlv ', ' basis]);
end

xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs Eb/N0');
legend('Location', 'southwest');
set(gca, 'YScale', 'log');
xlim([min(snr) max(snr)]);
ylim([1e-6 1]);

hold off;