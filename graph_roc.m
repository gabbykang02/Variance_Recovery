%% real data
N_neur = 338;
bin_size = 0.04; % sec
real = struct;
predict = struct;
rate = struct;
corr = struct;

real.ideal = readmatrix('~/Downloads/personalRepo/Cov/Data3/idealTraces.csv');
real.soma = readmatrix('~/Downloads/personalRepo/Cov/Data3/soma.csv');
real.evt1 = readmatrix('~/Downloads/personalRepo/Cov/Data3/evt.csv');
%% Real spike rates
frame_size = floor(600 / (20/bin_size)); % number of frames per bin_size
real.rate_evt = zeros([N_neur (600 / frame_size)]); % N_neur by number of bins matrix
for curr_neur = 1:N_neur
    for curr_frame = 1:frame_size:(600 - frame_size)
        real.rate_evt(curr_neur, floor(curr_frame / frame_size)) = sum((real.evt1(curr_neur, :) > curr_frame/30) & (real.evt1(curr_neur, :) < (curr_frame + frame_size)/30))/ bin_size;
    end
end
%% Load predictions
predict.a_soma = readmatrix('~/Downloads/personalRepo/Cov/Data3/mls_0.03_soma.csv');
predict.a_ideal = readmatrix('~/Downloads/personalRepo/Cov/Data3/mls_0.03_ideal.csv');

predict.suite_soma = readmatrix('~/Downloads/personalRepo/Cov/Data3/suite2P_soma.csv');
predict.suite_ideal = readmatrix('~/Downloads/personalRepo/Cov/Data3/suite2P_ideal.csv');
     predict.suite_soma(predict.suite_soma==0) = NaN;
     predict.suite_ideal(predict.suite_ideal==0) = NaN;
     predict.suite_soma = predict.suite_soma.';
     predict.suite_ideal = predict.suite_ideal.';
predict.five_soma = readmatrix('~/Downloads/personalRepo/Cov/Data3/5_soma_predict.csv');
predict.five_ideal = readmatrix('~/Downloads/personalRepo/Cov/Data3/5_ideal_predict.csv');
    predict.five_soma = predict.five_soma.';
    predict.five_ideal = predict.five_ideal.';
    predict.five_soma = predict.five_soma(2:size(predict.five_soma, 1), 2:size(predict.five_soma, 2));
    predict.five_ideal = predict.five_ideal(2:size(predict.five_ideal, 1), 2:size(predict.five_ideal, 2));
    predict.five_soma = [predict.five_soma zeros([size(predict.five_soma, 1) 1])];
    predict.five_ideal = [predict.five_ideal zeros([size(predict.five_ideal, 1) 1])];
predict.cascade_soma = load('~/Downloads/personalRepo/Cov/Data3/predictions_soma.mat');
predict.cascade_ideal = load('~/Downloads/personalRepo/Cov/Data3/predictions_idealTraces.mat');
predict.cascade_soma = predict.cascade_soma.spike_rates;
predict.cascade_ideal = predict.cascade_ideal.spike_rates;
%% mlspike spike rates
[predict.a_s_zero, predict.a_i_zero] = non_zero_trace(real.soma, real.ideal, N_neur);
[rate.rate_mls_a_s, rate.rate_mls_a_i] = mls_to_spike_rate(predict.a_soma,predict.a_ideal, N_neur, predict.a_s_zero, predict.a_i_zero, frame_size, bin_size);
%% suite2P spike rates
[rate.rate_suite_s,rate.rate_suite_i] =  deconv_to_spike_rate(predict.suite_soma, predict.suite_ideal, N_neur, bin_size, frame_size);
%[rate.rate_suite_s,rate.rate_suite_i] = to_spike_rate(predict.suite_soma,predict.suite_ideal, N_neur, bin_size, frame_size);
%[~,~] = prob_to_spike_rate(predict.suite_soma, predict.suite_ideal, N_neur, bin_size, frame_size);
%% cascade rates
[rate.rate_cas_s, rate.rate_cas_i] = prob_to_spike_rate(predict.cascade_soma, predict.cascade_ideal, N_neur, bin_size, frame_size);
%% 5 rates
[rate.rate_5_s, rate.rate_5_i] = prob_to_spike_rate(predict.five_soma,predict.five_ideal, N_neur, bin_size, frame_size);
%% calculate correlations between neurons
% real correlations between spike rates
[corr.corr_real] = real_cov(real.rate_evt);
% mlspike correlations between spike rates
[corr.corr_a_soma] = real_cov(rate.rate_mls_a_s);
[corr.corr_a_ideal] = real_cov(rate.rate_mls_a_i);
% cascade correlations between spike rates
% Requires conversion of NaN to 0
rate.rate_cas_s(isnan(rate.rate_cas_s)) = 0;
rate.rate_cas_i(isnan(rate.rate_cas_i)) = 0;
corr.corr_cas_soma = real_cov(rate.rate_cas_s);
corr.corr_cas_ideal = real_cov(rate.rate_cas_i);
%suite2p correlations between spike rates
corr.corr_suite_soma = real_cov(rate.rate_suite_s);
corr.corr_suite_ideal = real_cov(rate.rate_suite_i);
% five correlations between spike rates
corr.corr_5_s = real_cov(rate.rate_5_s);
corr.corr_5_i = real_cov(rate.rate_5_i);
%% cumulative_roc graph - soma
val = struct();
h = figure;
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
fields9 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.9);
fields10 = real_roc(corr.corr_real, corr.corr_cas_soma, 1.0);

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
fields00000 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.000001);
fields0000 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.00001);
fields000 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.0001);
fields00 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.001);
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

y = [fields00000.true fields0000.true fields000.true fields00.true fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true fields9.true fields10.true];
x = [fields00000.zero fields0000.zero fields000.zero fields00.zero fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero fields9.zero fields10.zero];
c = ['    '; '    '; '    '; '    ';'0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'; '0.90'; '1.00'];
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
% saveas(h, '~/Downloads/Data3_soma_roc.png', 'png');
%% cumulative_roc graph - idealTraces
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
fields00000 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.000001);
fields0000 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.00001);
fields000 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.0001);
fields00 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.001);
fields0 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.01);
fields1 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.1);
fields2 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.2);
fields3 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.3);
fields4 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.4);
fields5 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.5);
fields6 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.6);
fields7 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.7);
fields8 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.8);
fields9 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.9);
fields10 = real_roc(corr.corr_real, corr.corr_suite_ideal, 1.0);

y = [fields00000.true fields0000.true fields000.true fields00.true fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true fields8.true fields9.true fields10.true];
x = [fields00000.zero fields0000.zero fields000.zero fields00.zero fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero fields8.zero fields9.zero fields10.zero];

c = ['    '; '    '; '    '; '    ';'0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'; '0.80'; '0.90'; '1.00'];
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
title('Set: Data3 - ideal');
hold off
% savefig('Data3/ideal/roc_real_5');
% savefig('~/Downloads/personalRepo/Cov/Data3/ideal/cumulative_roc');
% saveas(h, '~/Downloads/Data3_ideal_roc.png', 'png');
