%% real data
N_neur = 338;
bin_size = 0.04; % sec
real = struct;
predict = struct;
rate = struct;
corr = struct;
norm = struct;
diff = struct;

real.ideal = readmatrix('~/Downloads/personalRepo/Cov/Data3/idealTraces.csv');
real.denoiseIdeal = readmatrix('~/Downloads/personalRepo/Cov/Data3/denoiseIdeal.csv');
real.soma = readmatrix('~/Downloads/personalRepo/Cov/Data3/soma.csv');
real.evt1 = readmatrix('~/Downloads/personalRepo/Cov/Data3/evt.csv');
real.locs = readmatrix('~/Downloads/personalRepo/Cov/Data3/locs.csv');
%% Real spike rates
frame_size = floor(600 / (20/bin_size)); % number of frames per bin_size
real.rate_evt = zeros([N_neur (600 / frame_size)]); % N_neur by number of bins matrix
for kk = 1:N_neur
    for jj = 1:frame_size:(600 - frame_size)
        real.rate_evt(kk, floor(jj / frame_size)) = sum((real.evt1(kk, :) > jj/30) & (real.evt1(kk, :) < (jj + frame_size)/30))/ bin_size;
    end
end
%% Load predictions
predict.a_soma = readmatrix('~/Downloads/personalRepo/Cov/Data3/mls_0.03_soma.csv');
predict.a_ideal = readmatrix('~/Downloads/personalRepo/Cov/Data3/mls_0.03_ideal.csv');
predict.a_ideal_denoise = readmatrix('~/Downloads/personalRepo/Cov/Data3/mls_0.03_idealdenoise.csv');

predict.suite_soma = readmatrix('~/Downloads/personalRepo/Cov/Data3/suite2P_soma.csv');
predict.suite_ideal = readmatrix('~/Downloads/personalRepo/Cov/Data3/suite2P_ideal.csv');
predict.suite_ideal_denoise = readmatrix('~/Downloads/personalRepo/Cov/Data3/suite2P_denoise.csv');
     predict.suite_soma(predict.suite_soma==0) = NaN;
     predict.suite_ideal(predict.suite_ideal==0) = NaN;
     predict.suite_ideal_denoise(predict.suite_ideal==0) = NaN;
     predict.suite_soma = predict.suite_soma.';
     predict.suite_ideal = predict.suite_ideal.';
     predict.suite_ideal_denoise = predict.suite_ideal_denoise.';
predict.five_soma = readmatrix('~/Downloads/personalRepo/Cov/Data3/5_soma_predict.csv');
predict.five_ideal = readmatrix('~/Downloads/personalRepo/Cov/Data3/5_ideal_predict.csv');
predict.five_ideal_denoise = readmatrix('~/Downloads/personalRepo/Cov/Data3/5_denoiseIdeal_predict.csv');
    predict.five_soma = predict.five_soma.';
    predict.five_ideal = predict.five_ideal.';
    predict.five_ideal_denoise = predict.five_ideal.';
    predict.five_soma = predict.five_soma(2:size(predict.five_soma, 1), 2:size(predict.five_soma, 2));
    predict.five_ideal = predict.five_ideal(2:size(predict.five_ideal, 1), 2:size(predict.five_ideal, 2));
    predict.five_ideal_denoise = predict.five_ideal_denoise(2:size(predict.five_ideal_denoise, 1), 2:size(predict.five_ideal_denoise, 2));
    predict.five_soma = [predict.five_soma zeros([size(predict.five_soma, 1) 1])];
    predict.five_ideal = [predict.five_ideal zeros([size(predict.five_ideal, 1) 1])];
    predict.five_ideal_denoise = [predict.five_ideal_denoise zeros([size(predict.five_ideal_denoise, 1) 1])];
predict.cascade_soma = load('~/Downloads/personalRepo/Cov/Data3/predictions_soma.mat');
predict.cascade_ideal = load('~/Downloads/personalRepo/Cov/Data3/predictions_idealTraces.mat');
predict.cascade_ideal_denoise = load('~/Downloads/personalRepo/Cov/Data3/predictions_denoiseIdeal.mat');
predict.cascade_soma = predict.cascade_soma.spike_rates;
predict.cascade_ideal = predict.cascade_ideal.spike_rates;
predict.cascade_ideal_denoise = predict.cascade_ideal_denoise.spike_rates;
%% mlspike spike rates
[predict.a_s_zero, predict.a_i_zero] = non_zero_trace(real.soma, real.ideal, N_neur);
[rate.rate_mls_a_s, rate.rate_mls_a_i] = mls_to_spike_rate(predict.a_soma,predict.a_ideal, N_neur, predict.a_s_zero, predict.a_i_zero, frame_size, bin_size);
[~, predict.a_denoisei_zero] = non_zero_trace(real.soma, real.denoiseIdeal, N_neur);
[~, rate.rate_mls_a_denoisei] = mls_to_spike_rate(predict.a_soma,predict.a_ideal_denoise, N_neur, predict.a_s_zero, predict.a_denoisei_zero, frame_size, bin_size);

%% suite2P spike rates
[rate.rate_suite_s,rate.rate_suite_i] =  deconv_to_spike_rate(predict.suite_soma, predict.suite_ideal, N_neur, bin_size, frame_size);
[~,rate.rate_suite_denoisei] =  deconv_to_spike_rate(predict.suite_soma, predict.suite_ideal_denoise, N_neur, bin_size, frame_size);

%[rate.rate_suite_s,rate.rate_suite_i] = to_spike_rate(predict.suite_soma,predict.suite_ideal, N_neur, bin_size, frame_size);
%[~,~] = prob_to_spike_rate(predict.suite_soma, predict.suite_ideal, N_neur, bin_size, frame_size);

%% cascade rates
[rate.rate_cas_s, rate.rate_cas_i] = prob_to_spike_rate(predict.cascade_soma, predict.cascade_ideal, N_neur, bin_size, frame_size);
[~, rate.rate_cas_denoisei] = prob_to_spike_rate(predict.cascade_soma, predict.cascade_ideal_denoise, N_neur, bin_size, frame_size);
%% 5 rates
[rate.rate_5_s, rate.rate_5_i] = prob_to_spike_rate(predict.five_soma,predict.five_ideal, N_neur, bin_size, frame_size);
[~, rate.rate_5_denoisei] = prob_to_spike_rate(predict.five_soma,predict.five_ideal_denoise, N_neur, bin_size, frame_size);
%% cov
% real correlations between spike rates
[corr.corr_real] = real_cov(real.rate_evt);
[corr.corr_a_soma] = real_cov(rate.rate_mls_a_s);
[corr.corr_a_ideal] = real_cov(rate.rate_mls_a_i);
[corr.corr_a_denoiseideal] = real_cov(rate.rate_mls_a_denoisei);
% cascade correlations between spike rates
rate.rate_cas_s(isnan(rate.rate_cas_s)) = 0;
rate.rate_cas_i(isnan(rate.rate_cas_i)) = 0;
rate.rate_cas_denoisei(isnan(rate.rate_cas_denoisei)) = 0;
corr.corr_cas_soma = real_cov(rate.rate_cas_s);
corr.corr_cas_ideal = real_cov(rate.rate_cas_i);
corr.corr_cas_denoiseideal = real_cov(rate.rate_cas_denoisei);
%suite2p correlations between spike rates
corr.corr_suite_soma = real_cov(rate.rate_suite_s);
corr.corr_suite_ideal = real_cov(rate.rate_suite_i);
corr.corr_suite_denoiseideal = real_cov(rate.rate_suite_denoisei);
% five correlations between spike rates
corr.corr_5_s = real_cov(rate.rate_5_s);
corr.corr_5_i = real_cov(rate.rate_5_i);
corr.corr_5_denoisei = real_cov(rate.rate_5_denoisei);
% comparing
[norm.norm_mls_a_s, norm.norm_mls_a_i] = cov_compare(corr.corr_real, corr.corr_a_soma, corr.corr_a_ideal, N_neur, real.rate_evt, rate.rate_mls_a_s, rate.rate_mls_a_i);
[diff.diff_mls_a_s, diff.diff_mls_a_i] = cov_diff(corr.corr_real, corr.corr_a_soma, corr.corr_a_ideal, N_neur);

[norm.norm_cas_s, norm.norm_cas_i] = cov_compare(corr.corr_real, corr.corr_cas_soma, corr.corr_cas_ideal, N_neur, real.rate_evt, rate.rate_cas_s, rate.rate_cas_i);
[diff.diff_cas_s, diff.diff_cas_i] = cov_diff(corr.corr_real, corr.corr_cas_soma, corr.corr_cas_ideal, N_neur);

[norm.norm_suite_s, norm.norm_suite_i] = cov_compare(corr.corr_real, corr.corr_suite_soma, corr.corr_suite_ideal, N_neur, real.rate_evt, rate.rate_suite_s, rate.rate_suite_denoisei);
[diff.diff_suite_s, diff.diff_suite_i] = cov_diff(corr.corr_real, corr.corr_suite_soma, corr.corr_suite_ideal, N_neur);

[norm.norm_5_s, norm.norm_5_i] = cov_compare(corr.corr_real, corr.corr_5_s, corr.corr_5_i, N_neur, real.rate_evt, rate.rate_5_s, rate.rate_5_i);
[diff.diff_5_s, diff.diff_5_i] = cov_diff(corr.corr_real, corr.corr_5_s, corr.corr_5_i, N_neur);
% comparing w/ deenoise
[~, norm.norm_mls_a_denoisei] = cov_compare(corr.corr_real, corr.corr_a_soma, corr.corr_a_denoiseideal, N_neur, real.rate_evt, rate.rate_mls_a_s, rate.rate_mls_a_denoisei);
[~, diff.diff_mls_a_denoisei] = cov_diff(corr.corr_real, corr.corr_a_soma, corr.corr_a_denoiseideal, N_neur);

[~, norm.norm_cas_denoisei] = cov_compare(corr.corr_real, corr.corr_cas_soma, corr.corr_cas_denoiseideal, N_neur, real.rate_evt, rate.rate_cas_s, rate.rate_cas_denoisei);
[~, diff.diff_cas_denoisei] = cov_diff(corr.corr_real, corr.corr_cas_soma, corr.corr_cas_denoiseideal, N_neur);

[~, norm.norm_suite_denoisei] = cov_compare(corr.corr_real, corr.corr_suite_soma, corr.corr_suite_denoiseideal, N_neur, real.rate_evt, rate.rate_suite_s, rate.rate_suite_denoisei);
[~, diff.diff_suite_denoisei] = cov_diff(corr.corr_real, corr.corr_suite_soma, corr.corr_suite_denoiseideal, N_neur);

[~, norm.norm_5_denoisei] = cov_compare(corr.corr_real, corr.corr_5_s, corr.corr_5_denoisei, N_neur, real.rate_evt, rate.rate_5_s, rate.rate_5_denoisei);
[~, diff.diff_5_denoisei] = cov_diff(corr.corr_real, corr.corr_5_s, corr.corr_5_denoisei, N_neur);
%% roc id ideal
h = figure;
id = isolateVisibleSomaMod(real.locs, 0);
val = struct();
fields00 = real_roc_id(corr.corr_real, corr.corr_a_denoiseideal, 0.0001, id);
fields01 = real_roc_id(corr.corr_real, corr.corr_a_denoiseideal, 0.001, id);
fields0 = real_roc_id(corr.corr_real, corr.corr_a_denoiseideal, 0.01, id);
fields1 = real_roc_id(corr.corr_real, corr.corr_a_denoiseideal, 0.1, id);
fields2 = real_roc_id(corr.corr_real, corr.corr_a_denoiseideal, 0.2, id);
fields3 = real_roc_id(corr.corr_real, corr.corr_a_denoiseideal, 0.3, id);
fields4 = real_roc_id(corr.corr_real, corr.corr_a_denoiseideal, 0.4, id);
fields5 = real_roc_id(corr.corr_real, corr.corr_a_denoiseideal, 0.5, id);
fields6 = real_roc_id(corr.corr_real, corr.corr_a_denoiseideal, 0.6, id);
fields7 = real_roc_id(corr.corr_real, corr.corr_a_denoiseideal, 0.7, id);
fields8 = real_roc_id(corr.corr_real, corr.corr_a_denoiseideal, 0.8, id);
y = [fields00.true fields01.true fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields00.zero fields01.zero fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0001'; '.001'; '0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];

val.mls_i= trapz(flip(x), flip(y));
scatter(x, y, [], 'red', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);
hold on
xlim([0, 1]);
ylim([0, 1]);
plot(x, y, 'Color', 'red');
ylabel('True Positive');
xlabel('False Positive');
hold off
% savefig('Data3/ideal/roc_real_mls');
fields0 = real_roc_id(corr.corr_real, corr.corr_cas_denoiseideal, 0.01, id);
fields1 = real_roc_id(corr.corr_real, corr.corr_cas_denoiseideal, 0.1, id);
fields2 = real_roc_id(corr.corr_real, corr.corr_cas_denoiseideal, 0.2, id);
fields3 = real_roc_id(corr.corr_real, corr.corr_cas_denoiseideal, 0.3, id);
fields4 = real_roc_id(corr.corr_real, corr.corr_cas_denoiseideal, 0.4, id);
fields5 = real_roc_id(corr.corr_real, corr.corr_cas_denoiseideal, 0.5, id);
fields6 = real_roc_id(corr.corr_real, corr.corr_cas_denoiseideal, 0.6, id);
fields7 = real_roc_id(corr.corr_real, corr.corr_cas_denoiseideal, 0.7, id);
fields8 = real_roc_id(corr.corr_real, corr.corr_cas_denoiseideal, 0.8, id);
y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];

hold on
val.cas_i= trapz(flip(x), flip(y));
scatter(x, y, [], 'green', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);
% hold on
xlim([0, 1]);
ylim([0, 1]);
plot(x, y, 'Color', 'green');
ylabel('True Positive');
xlabel('False Positive');
hold off
% savefig('Data3/ideal/roc_real_cas');
fields00 = real_roc_id(corr.corr_real, corr.corr_suite_denoiseideal, 0.0001, id);
fields01 = real_roc_id(corr.corr_real, corr.corr_suite_denoiseideal, 0.001, id);
fields0 = real_roc_id(corr.corr_real, corr.corr_suite_denoiseideal, 0.01, id);
fields1 = real_roc_id(corr.corr_real, corr.corr_suite_denoiseideal, 0.1, id);
fields2 = real_roc_id(corr.corr_real, corr.corr_suite_denoiseideal, 0.2, id);
fields3 = real_roc_id(corr.corr_real, corr.corr_suite_denoiseideal, 0.3, id);
fields4 = real_roc_id(corr.corr_real, corr.corr_suite_denoiseideal, 0.4, id);
fields5 = real_roc_id(corr.corr_real, corr.corr_suite_denoiseideal, 0.5, id);
fields6 = real_roc_id(corr.corr_real, corr.corr_suite_denoiseideal, 0.6, id);
fields7 = real_roc_id(corr.corr_real, corr.corr_suite_denoiseideal, 0.7, id);
fields8 = real_roc_id(corr.corr_real, corr.corr_suite_denoiseideal, 0.8, id);
y = [fields00.true fields01.true fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields00.zero fields01.zero fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0001'; '.001'; '0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];
hold on
val.suite_i= trapz(flip(x), flip(y));
scatter(x, y, [], 'blue', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);
% hold on
xlim([0, 1]);
ylim([0, 1]);
plot(x, y), 'Color', 'blue';
ylabel('True Positive');
xlabel('False Positive');
hold off
% savefig('Data3/ideal/roc_real_suite');
fields0 = real_roc_id(corr.corr_real, corr.corr_5_denoisei, 0.01, id);
fields1 = real_roc_id(corr.corr_real, corr.corr_5_denoisei, 0.1, id);
fields2 = real_roc_id(corr.corr_real, corr.corr_5_denoisei, 0.2, id);
fields3 = real_roc_id(corr.corr_real, corr.corr_5_denoisei, 0.3, id);
fields4 = real_roc_id(corr.corr_real, corr.corr_5_denoisei, 0.4, id);
fields5 = real_roc_id(corr.corr_real, corr.corr_5_denoisei, 0.5, id);
fields6 = real_roc_id(corr.corr_real, corr.corr_5_denoisei, 0.6, id);
fields7 = real_roc_id(corr.corr_real, corr.corr_5_denoisei, 0.7, id);
fields8 = real_roc_id(corr.corr_real, corr.corr_5_denoisei, 0.8, id);
y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];
hold on
val.five_i= trapz(flip(x), flip(y));
scatter(x, y, [], 'magenta', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);
hold on
xlim([0, 1]);
ylim([0, 1]);
plot(x, y, 'Color', 'magenta');
ylabel('True Positive');
xlabel('False Positive');
legend( 'MLSpike',  'Cascade', 'Suite2P', 'Five', 'Location', 'southeast');
title('Set: Data3 - denoised ideal w/ visible somas');
hold off
% savefig('Data3/ideal/roc_real_5');
savefig('~/Downloads/personalRepo/Cov/Data3/ideal/cumulative_roc_ideal_denoise');
saveas(h, '~/Downloads/Data3_idealdenoise_visible_roc.png', 'png');
%% roc id ideal
h = figure;
id = isolateVisibleSomaMod(real.locs, 0);
val = struct();
fields00 = real_roc_id(corr.corr_real, corr.corr_a_ideal, 0.0001, id);
fields01 = real_roc_id(corr.corr_real, corr.corr_a_ideal, 0.001, id);
fields0 = real_roc_id(corr.corr_real, corr.corr_a_ideal, 0.01, id);
fields1 = real_roc_id(corr.corr_real, corr.corr_a_ideal, 0.1, id);
fields2 = real_roc_id(corr.corr_real, corr.corr_a_ideal, 0.2, id);
fields3 = real_roc_id(corr.corr_real, corr.corr_a_ideal, 0.3, id);
fields4 = real_roc_id(corr.corr_real, corr.corr_a_ideal, 0.4, id);
fields5 = real_roc_id(corr.corr_real, corr.corr_a_ideal, 0.5, id);
fields6 = real_roc_id(corr.corr_real, corr.corr_a_ideal, 0.6, id);
fields7 = real_roc_id(corr.corr_real, corr.corr_a_ideal, 0.7, id);
fields8 = real_roc_id(corr.corr_real, corr.corr_a_ideal, 0.8, id);
y = [fields00.true fields01.true fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields00.zero fields01.zero fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0001'; '.001'; '0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];

val.mls_i= trapz(flip(x), flip(y));
scatter(x, y, [], 'red', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);
hold on
xlim([0, 1]);
ylim([0, 1]);
plot(x, y, 'Color', 'red');
ylabel('True Positive');
xlabel('False Positive');
hold off
% savefig('Data3/ideal/roc_real_mls');
fields0 = real_roc_id(corr.corr_real, corr.corr_cas_ideal, 0.01, id);
fields1 = real_roc_id(corr.corr_real, corr.corr_cas_ideal, 0.1, id);
fields2 = real_roc_id(corr.corr_real, corr.corr_cas_ideal, 0.2, id);
fields3 = real_roc_id(corr.corr_real, corr.corr_cas_ideal, 0.3, id);
fields4 = real_roc_id(corr.corr_real, corr.corr_cas_ideal, 0.4, id);
fields5 = real_roc_id(corr.corr_real, corr.corr_cas_ideal, 0.5, id);
fields6 = real_roc_id(corr.corr_real, corr.corr_cas_ideal, 0.6, id);
fields7 = real_roc_id(corr.corr_real, corr.corr_cas_ideal, 0.7, id);
fields8 = real_roc_id(corr.corr_real, corr.corr_cas_ideal, 0.8, id);
y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];

hold on
val.cas_i= trapz(flip(x), flip(y));
scatter(x, y, [], 'green', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);
% hold on
xlim([0, 1]);
ylim([0, 1]);
plot(x, y, 'Color', 'green');
ylabel('True Positive');
xlabel('False Positive');
hold off
% savefig('Data3/ideal/roc_real_cas');
fields00 = real_roc_id(corr.corr_real, corr.corr_suite_ideal, 0.0001, id);
fields01 = real_roc_id(corr.corr_real, corr.corr_suite_ideal, 0.001, id);
fields0 = real_roc_id(corr.corr_real, corr.corr_suite_ideal, 0.01, id);
fields1 = real_roc_id(corr.corr_real, corr.corr_suite_ideal, 0.1, id);
fields2 = real_roc_id(corr.corr_real, corr.corr_suite_ideal, 0.2, id);
fields3 = real_roc_id(corr.corr_real, corr.corr_suite_ideal, 0.3, id);
fields4 = real_roc_id(corr.corr_real, corr.corr_suite_ideal, 0.4, id);
fields5 = real_roc_id(corr.corr_real, corr.corr_suite_ideal, 0.5, id);
fields6 = real_roc_id(corr.corr_real, corr.corr_suite_ideal, 0.6, id);
fields7 = real_roc_id(corr.corr_real, corr.corr_suite_ideal, 0.7, id);
fields8 = real_roc_id(corr.corr_real, corr.corr_suite_ideal, 0.8, id);
y = [fields00.true fields01.true fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields00.zero fields01.zero fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0001'; '.001'; '0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];
hold on
val.suite_i= trapz(flip(x), flip(y));
scatter(x, y, [], 'blue', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);
% hold on
xlim([0, 1]);
ylim([0, 1]);
plot(x, y), 'Color', 'blue';
ylabel('True Positive');
xlabel('False Positive');
hold off
% savefig('Data3/ideal/roc_real_suite');
fields0 = real_roc_id(corr.corr_real, corr.corr_5_i, 0.01, id);
fields1 = real_roc_id(corr.corr_real, corr.corr_5_i, 0.1, id);
fields2 = real_roc_id(corr.corr_real, corr.corr_5_i, 0.2, id);
fields3 = real_roc_id(corr.corr_real, corr.corr_5_i, 0.3, id);
fields4 = real_roc_id(corr.corr_real, corr.corr_5_i, 0.4, id);
fields5 = real_roc_id(corr.corr_real, corr.corr_5_i, 0.5, id);
fields6 = real_roc_id(corr.corr_real, corr.corr_5_i, 0.6, id);
fields7 = real_roc_id(corr.corr_real, corr.corr_5_i, 0.7, id);
fields8 = real_roc_id(corr.corr_real, corr.corr_5_i, 0.8, id);
y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];
hold on
val.five_i= trapz(flip(x), flip(y));
scatter(x, y, [], 'magenta', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);
hold on
xlim([0, 1]);
ylim([0, 1]);
plot(x, y, 'Color', 'magenta');
ylabel('True Positive');
xlabel('False Positive');
legend( 'MLSpike',  'Cascade', 'Suite2P', 'Five', 'Location', 'southeast');
title('Set: Data3 - ideal w/ visible somas');
hold off
% savefig('Data3/ideal/roc_real_5');
savefig('~/Downloads/personalRepo/Cov/Data3/ideal/cumulative_roc_ideal_visible');
saveas(h, '~/Downloads/Data3_ideal_visible_roc.png', 'png');

%% roc mod
h = figure;
val = struct();
fields00 = real_roc(corr.corr_real, corr.corr_a_soma, 0.0001);
fields01 = real_roc(corr.corr_real, corr.corr_a_soma, 0.001);
fields0 = real_roc(corr.corr_real, corr.corr_a_soma, 0.01);
fields1 = real_roc(corr.corr_real, corr.corr_a_soma, 0.1);
fields2 = real_roc(corr.corr_real, corr.corr_a_soma, 0.2);
fields3 = real_roc(corr.corr_real, corr.corr_a_soma, 0.3);
fields4 = real_roc(corr.corr_real, corr.corr_a_soma, 0.4);
fields5 = real_roc(corr.corr_real, corr.corr_a_soma, 0.5);
fields6 = real_roc(corr.corr_real, corr.corr_a_soma, 0.6);
fields7 = real_roc(corr.corr_real, corr.corr_a_soma, 0.7);
fields8 = real_roc(corr.corr_real, corr.corr_a_soma, 0.8);
y = [fields00.true fields01.true fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields00.zero fields01.zero fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0001'; '.001'; '0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];
val.mls_s= trapz(flip(x), flip(y));
scatter(x, y, [], 'red', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);

hold on

xlim([0, 1]);
ylim([0, 1]);
plot(x, y, 'Color', 'red');
ylabel('True Positive');
xlabel('False Positive');
hold off
%savefig('Data3/soma/roc_real_mls');

fields0 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.01);
fields1 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.1);
fields2 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.2);
fields3 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.3);
fields4 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.4);
fields5 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.5);
fields6 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.6);
fields7 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.7);
fields8 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.8);
fields9 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.9);
fields10 = real_roc(corr.corr_real, corr.corr_suite_soma, 1.0);

y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true fields9.true fields10.true];
x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero fields9.zero fields10.zero];
c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'; '0.90'; '1.00'];
hold on
val.cas_s= trapz(flip(x), flip(y));
scatter(x, y, [], 'green','filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);

xlim([0, 1]);
ylim([0, 1]);
plot(x, y, 'Color', 'green');
ylabel('True Positive');
xlabel('False Positive');
hold off
%savefig('Data3/soma/roc_real_cas');
fields0 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.01);
fields1 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.1);
fields2 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.2);
fields3 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.3);
fields4 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.4);
fields5 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.5);
fields6 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.6);
fields7 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.7);
fields8 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.8);
fields9 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.9);
fields10 = real_roc(corr.corr_real, corr.corr_suite_soma, 1.0);

y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];
hold on
val.suite_s= trapz(flip(x), flip(y));
scatter(x, y,[], 'blue', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);

xlim([0, 1]);
ylim([0, 1]);
plot(x, y, 'Color', 'blue');
ylabel('True Positive');
xlabel('False Positive');
hold off
%savefig('Data3/soma/roc_real_suite');
fields0 = real_roc(corr.corr_real, corr.corr_5_s, 0.01);
fields1 = real_roc(corr.corr_real, corr.corr_5_s, 0.1);
fields2 = real_roc(corr.corr_real, corr.corr_5_s, 0.2);
fields3 = real_roc(corr.corr_real, corr.corr_5_s, 0.3);
fields4 = real_roc(corr.corr_real, corr.corr_5_s, 0.4);
fields5 = real_roc(corr.corr_real, corr.corr_5_s, 0.5);
fields6 = real_roc(corr.corr_real, corr.corr_5_s, 0.6);
fields7 = real_roc(corr.corr_real, corr.corr_5_s, 0.7);
fields8 = real_roc(corr.corr_real, corr.corr_5_s, 0.8);
y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];

hold on
val.five_s= trapz(flip(x), flip(y));
scatter(x, y, [], 'magenta', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);
% hold on
xlim([0, 1]);
ylim([0, 1]);
plot(x, y, 'Color', 'magenta');
ylabel('True Positive');
xlabel('False Positive');
legend( 'MLSpike',  'Cascade', 'Suite2P', 'Five', 'Location', 'southeast');
title('Set: Data3 - soma');
hold off
% savefig('Data3/soma/roc_real_5');
% savefig('~/Downloads/personalRepo/Cov/Data3/soma/cumulative_roc');

%% roc mod ideal
h = figure;
fields00 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.0001);
fields01 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.001);
fields0 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.01);
fields1 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.1);
fields2 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.2);
fields3 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.3);
fields4 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.4);
fields5 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.5);
fields6 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.6);
fields7 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.7);
fields8 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.8);
y = [fields00.true fields01.true fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields00.zero fields01.zero fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0001'; '.001'; '0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];

val.mls_i= trapz(flip(x), flip(y));
scatter(x, y, [], 'red', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);
hold on
xlim([0, 1]);
ylim([0, 1]);
plot(x, y, 'Color', 'red');
ylabel('True Positive');
xlabel('False Positive');
hold off
% savefig('Data3/ideal/roc_real_mls');
fields0 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.01);
fields1 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.1);
fields2 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.2);
fields3 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.3);
fields4 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.4);
fields5 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.5);
fields6 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.6);
fields7 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.7);
fields8 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.8);
y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];

hold on
val.cas_i= trapz(flip(x), flip(y));
scatter(x, y, [], 'green', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);
% hold on
xlim([0, 1]);
ylim([0, 1]);
plot(x, y, 'Color', 'green');
ylabel('True Positive');
xlabel('False Positive');
hold off
% savefig('Data3/ideal/roc_real_cas');
fields00 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.0001);
fields01 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.001);
fields0 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.01);
fields1 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.1);
fields2 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.2);
fields3 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.3);
fields4 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.4);
fields5 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.5);
fields6 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.6);
fields7 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.7);
fields8 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.8);
y = [fields00.true fields01.true fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields00.zero fields01.zero fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0001'; '.001'; '0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];
hold on
val.suite_i= trapz(flip(x), flip(y));
scatter(x, y, [], 'blue', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);
% hold on
xlim([0, 1]);
ylim([0, 1]);
plot(x, y), 'Color', 'blue';
ylabel('True Positive');
xlabel('False Positive');
hold off
% savefig('Data3/ideal/roc_real_suite');
fields0 = real_roc(corr.corr_real, corr.corr_5_i, 0.01);
fields1 = real_roc(corr.corr_real, corr.corr_5_i, 0.1);
fields2 = real_roc(corr.corr_real, corr.corr_5_i, 0.2);
fields3 = real_roc(corr.corr_real, corr.corr_5_i, 0.3);
fields4 = real_roc(corr.corr_real, corr.corr_5_i, 0.4);
fields5 = real_roc(corr.corr_real, corr.corr_5_i, 0.5);
fields6 = real_roc(corr.corr_real, corr.corr_5_i, 0.6);
fields7 = real_roc(corr.corr_real, corr.corr_5_i, 0.7);
fields8 = real_roc(corr.corr_real, corr.corr_5_i, 0.8);
y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true];
x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero];
c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'];
hold on
val.five_i= trapz(flip(x), flip(y));
scatter(x, y, [], 'magenta', 'filled', 'HandleVisibility', 'off');
text(x+0.001, y + 0.02, c);
hold on
xlim([0, 1]);
ylim([0, 1]);
plot(x, y, 'Color', 'magenta');
ylabel('True Positive');
xlabel('False Positive');
legend( 'MLSpike',  'Cascade', 'Suite2P', 'Five', 'Location', 'southeast');
title('Set: Data3 - denoised ideal');
hold off
% savefig('Data3/ideal/roc_real_5');
savefig('~/Downloads/personalRepo/Cov/Data3/ideal/cumulative_roc_ideal');
saveas(h, '~/Downloads/Data3_idealdenoise_roc.png', 'png');
