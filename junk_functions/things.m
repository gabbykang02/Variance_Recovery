%% real data
N_neur = 338;
bin_size = 0.04; % sec
real = struct;
predict = struct;
rate = struct;
corr = struct;
norm = struct;
diff = struct;

real.ideal = readmatrix('~/Downloads/personalRepo/Cov/~/Downloads/personalRepo/Cov/Data3/idealTraces.csv');
real.soma = readmatrix('~/Downloads/personalRepo/Cov/~/Downloads/personalRepo/Cov/Data3/soma.csv');
real.evt1 = readmatrix('~/Downloads/personalRepo/Cov/~/Downloads/personalRepo/Cov/Data3/evt.csv');
%% Real spike rates
frame_size = floor(600 / (20/bin_size)); % number of frames per bin_size
real.rate_evt = zeros([N_neur (600 / frame_size)]); % N_neur by number of bins matrix
for kk = 1:N_neur
    for jj = 1:frame_size:(600 - frame_size)
        real.rate_evt(kk, floor(jj / frame_size)) = sum((real.evt1(kk, :) > jj/30) & (real.evt1(kk, :) < (jj + frame_size)/30))/ bin_size;
    end
end
%% Load predictions
predict.a_soma = readmatrix('~/Downloads/personalRepo/Cov/~/Downloads/personalRepo/Cov/Data3/mls_0.03_soma.csv');
predict.a_ideal = readmatrix('~/Downloads/personalRepo/Cov/~/Downloads/personalRepo/Cov/Data3/mls_0.03_ideal.csv');

predict.suite_soma = readmatrix('~/Downloads/personalRepo/Cov/~/Downloads/personalRepo/Cov/Data3/suite2P_soma.csv');
predict.suite_ideal = readmatrix('~/Downloads/personalRepo/Cov/~/Downloads/personalRepo/Cov/Data3/suite2P_ideal.csv');
     predict.suite_soma(predict.suite_soma==0) = NaN;
     predict.suite_ideal(predict.suite_ideal==0) = NaN;
     predict.suite_soma = predict.suite_soma.';
     predict.suite_ideal = predict.suite_ideal.';
predict.five_soma = readmatrix('~/Downloads/personalRepo/Cov/~/Downloads/personalRepo/Cov/Data3/5_soma_predict.csv');
predict.five_ideal = readmatrix('~/Downloads/personalRepo/Cov/~/Downloads/personalRepo/Cov/Data3/5_ideal_predict.csv');
    predict.five_soma = predict.five_soma.';
    predict.five_ideal = predict.five_ideal.';
    predict.five_soma = predict.five_soma(2:size(predict.five_soma, 1), 2:size(predict.five_soma, 2));
    predict.five_ideal = predict.five_ideal(2:size(predict.five_ideal, 1), 2:size(predict.five_ideal, 2));
    predict.five_soma = [predict.five_soma zeros([size(predict.five_soma, 1) 1])];
    predict.five_ideal = [predict.five_ideal zeros([size(predict.five_ideal, 1) 1])];
predict.cascade_soma = load('~/Downloads/personalRepo/Cov/~/Downloads/personalRepo/Cov/Data3/predictions_soma.mat');
predict.cascade_ideal = load('~/Downloads/personalRepo/Cov/~/Downloads/personalRepo/Cov/Data3/predictions_idealTraces.mat');
predict.cascade_soma = predict.cascade_soma.spike_rates;
predict.cascade_ideal = predict.cascade_ideal.spike_rates;
%% mlspike spike rates
[predict.a_s_zero, predict.a_i_zero] = non_zero_trace(real.soma, real.ideal, N_neur);
[rate.rate_mls_a_s, rate.rate_mls_a_i] = mls_to_spike_rate(predict.a_soma,predict.a_ideal, N_neur, predict.a_s_zero, predict.a_i_zero, frame_size, bin_size);
%% suite2P spike rates
[rate.rate_suite_s,rate.rate_suite_i] = to_spike_rate(predict.suite_soma,predict.suite_ideal, N_neur, bin_size, frame_size);
%[rate.rate_suite_s,rate.rate_suite_i] = prob_to_spike_rate(predict.suite_soma, predict.suite_ideal, N_neur, bin_size, frame_size);
%% cascade rates
[rate.rate_cas_s, rate.rate_cas_i] = prob_to_spike_rate(predict.cascade_soma, predict.cascade_ideal, N_neur, bin_size, frame_size);
%% 5 rates
[rate.rate_5_s, rate.rate_5_i] = prob_to_spike_rate(predict.five_soma,predict.five_ideal, N_neur, bin_size, frame_size);
%% cov
% real correlations between spike rates
[corr.corr_real] = real_cov(real.rate_evt);
[corr.corr_a_soma] = real_cov(rate.rate_mls_a_s);
[corr.corr_a_ideal] = real_cov(rate.rate_mls_a_i);
% cascade correlations between spike rates
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
% comparing
[norm.norm_mls_a_s, norm.norm_mls_a_i] = cov_compare(corr.corr_real, corr.corr_a_soma, corr.corr_a_ideal, N_neur, real.rate_evt, rate.rate_mls_a_s, rate.rate_mls_a_i);
[diff.diff_mls_a_s, diff.diff_mls_a_i] = cov_diff(corr.corr_real, corr.corr_a_soma, corr.corr_a_ideal, N_neur);

[norm.norm_cas_s, norm.norm_cas_i] = cov_compare(corr.corr_real, corr.corr_cas_soma, corr.corr_cas_ideal, N_neur, real.rate_evt, rate.rate_cas_s, rate.rate_cas_i);
[diff.diff_cas_s, diff.diff_cas_i] = cov_diff(corr.corr_real, corr.corr_cas_soma, corr.corr_cas_ideal, N_neur);

[norm.norm_suite_s, norm.norm_suite_i] = cov_compare(corr.corr_real, corr.corr_suite_soma, corr.corr_suite_ideal, N_neur, real.rate_evt, rate.rate_suite_s, rate.rate_suite_i);
[diff.diff_suite_s, diff.diff_suite_i] = cov_diff(corr.corr_real, corr.corr_suite_soma, corr.corr_suite_ideal, N_neur);

[norm.norm_5_s, norm.norm_5_i] = cov_compare(corr.corr_real, corr.corr_5_s, corr.corr_5_i, N_neur, real.rate_evt, rate.rate_5_s, rate.rate_5_i);
[diff.diff_5_s, diff.diff_5_i] = cov_diff(corr.corr_real, corr.corr_5_s, corr.corr_5_i, N_neur);
%% roc mod
% val = struct();
% fields0 = real_roc(corr.corr_real, corr.corr_a_soma, 0.01);
% fields1 = real_roc(corr.corr_real, corr.corr_a_soma, 0.1);
% fields2 = real_roc(corr.corr_real, corr.corr_a_soma, 0.2);
% fields3 = real_roc(corr.corr_real, corr.corr_a_soma, 0.3);
% fields4 = real_roc(corr.corr_real, corr.corr_a_soma, 0.4);
% fields5 = real_roc(corr.corr_real, corr.corr_a_soma, 0.5);
% fields6 = real_roc(corr.corr_real, corr.corr_a_soma, 0.6);
% fields7 = real_roc(corr.corr_real, corr.corr_a_soma, 0.7);
% y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true];
% x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero];
% c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'];
% 
% scatter(x, y, 'filled');
% text(x+0.001, y + 0.02, c);
% 
% hold on
% 
% xlim([0, 1]);
% ylim([0, 1]);
% plot(x, y, 'DisplayName', 'MLSpike');
% ylabel('True Positive');
% xlabel('False Positive');
% hold off
% %savefig('~/Downloads/personalRepo/Cov/Data3/soma/roc_real_mls');
% fields0 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.01);
% fields1 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.1);
% fields2 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.2);
% fields3 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.3);
% fields4 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.4);
% fields5 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.5);
% fields6 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.6);
% fields7 = real_roc(corr.corr_real, corr.corr_cas_soma, 0.7);
% y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true];
% x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero];
% c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'];
% hold on
% scatter(x, y, 'filled');
% text(x+0.001, y + 0.02, c);
% 
% xlim([0, 1]);
% ylim([0, 1]);
% plot(x, y, 'DisplayName', 'Cascade');
% ylabel('True Positive');
% xlabel('False Positive');
% hold off
% %savefig('~/Downloads/personalRepo/Cov/Data3/soma/roc_real_cas');
% fields0 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.01);
% fields1 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.1);
% fields2 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.2);
% fields3 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.3);
% fields4 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.4);
% fields5 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.5);
% fields6 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.6);
% fields7 = real_roc(corr.corr_real, corr.corr_suite_soma, 0.7);
% y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true];
% x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero];
% c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'];
% hold on
% scatter(x, y, 'filled');
% text(x+0.001, y + 0.02, c);
% 
% xlim([0, 1]);
% ylim([0, 1]);
% plot(x, y, 'DisplayName', 'Suite2P');
% ylabel('True Positive');
% xlabel('False Positive');
% hold off
% %savefig('~/Downloads/personalRepo/Cov/Data3/soma/roc_real_suite');
% fields0 = real_roc(corr.corr_real, corr.corr_5_s, 0.01);
% fields1 = real_roc(corr.corr_real, corr.corr_5_s, 0.1);
% fields2 = real_roc(corr.corr_real, corr.corr_5_s, 0.2);
% fields3 = real_roc(corr.corr_real, corr.corr_5_s, 0.3);
% fields4 = real_roc(corr.corr_real, corr.corr_5_s, 0.4);
% fields5 = real_roc(corr.corr_real, corr.corr_5_s, 0.5);
% fields6 = real_roc(corr.corr_real, corr.corr_5_s, 0.6);
% fields7 = real_roc(corr.corr_real, corr.corr_5_s, 0.7);
% y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true];
% x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero];
% c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'];
% 
% hold on
% scatter(x, y, 'filled');
% text(x+0.001, y + 0.02, c);
% % hold on
% xlim([0, 1]);
% ylim([0, 1]);
% plot(x, y, 'DisplayName', 'Five');
% ylabel('True Positive');
% xlabel('False Positive');
% legend('Location', 'southeast');
% hold off
% % savefig('~/Downloads/personalRepo/Cov/Data3/soma/roc_real_5');
% savefig('~/Downloads/personalRepo/Cov/Data3/soma/cumulative_roc');
% %% roc mod ideal
% fields0 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.01);
% fields1 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.1);
% fields2 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.2);
% fields3 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.3);
% fields4 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.4);
% fields5 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.5);
% fields6 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.6);
% fields7 = real_roc(corr.corr_real, corr.corr_a_ideal, 0.7);
% y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true];
% x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero];
% c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'];
% 
% scatter(x, y, 'filled');
% legend('Location', 'southeast');
% text(x+0.001, y + 0.02, c);
% hold on
% xlim([0, 1]);
% ylim([0, 1]);
% plot(x, y, 'DisplayName', 'MLSpike');
% ylabel('True Positive');
% xlabel('False Positive');
% hold off
% % savefig('~/Downloads/personalRepo/Cov/Data3/ideal/roc_real_mls');
% fields0 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.01);
% fields1 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.1);
% fields2 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.2);
% fields3 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.3);
% fields4 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.4);
% fields5 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.5);
% fields6 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.6);
% fields7 = real_roc(corr.corr_real, corr.corr_cas_ideal, 0.7);
% y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true];
% x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero];
% c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'];
% 
% hold on
% scatter(x, y, 'filled');
% text(x+0.001, y + 0.02, c);
% % hold on
% xlim([0, 1]);
% ylim([0, 1]);
% plot(x, y, 'DisplayName', 'Cascade');
% ylabel('True Positive');
% xlabel('False Positive');
% hold off
% % savefig('~/Downloads/personalRepo/Cov/Data3/ideal/roc_real_cas');
% fields0 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.01);
% fields1 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.1);
% fields2 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.2);
% fields3 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.3);
% fields4 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.4);
% fields5 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.5);
% fields6 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.6);
% fields7 = real_roc(corr.corr_real, corr.corr_suite_ideal, 0.7);
% y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true];
% x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero];
% c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'];
% 
% hold on
% scatter(x, y, 'filled');
% text(x+0.001, y + 0.02, c);
% % hold on
% xlim([0, 1]);
% ylim([0, 1]);
% plot(x, y, 'DisplayName', 'Suite2P');
% ylabel('True Positive');
% xlabel('False Positive');
% hold off
% % savefig('~/Downloads/personalRepo/Cov/Data3/ideal/roc_real_suite');
% fields0 = real_roc(corr.corr_real, corr.corr_5_i, 0.01);
% fields1 = real_roc(corr.corr_real, corr.corr_5_i, 0.1);
% fields2 = real_roc(corr.corr_real, corr.corr_5_i, 0.2);
% fields3 = real_roc(corr.corr_real, corr.corr_5_i, 0.3);
% fields4 = real_roc(corr.corr_real, corr.corr_5_i, 0.4);
% fields5 = real_roc(corr.corr_real, corr.corr_5_i, 0.5);
% fields6 = real_roc(corr.corr_real, corr.corr_5_i, 0.6);
% fields7 = real_roc(corr.corr_real, corr.corr_5_i, 0.7);
% y = [fields0.true fields1.true fields2.true fields3.true fields4.true fields5.true fields6.true fields7.true];
% x = [fields0.zero fields1.zero fields2.zero fields3.zero fields4.zero fields5.zero fields6.zero fields7.zero];
% c = ['0.01'; '0.10'; '0.20'; '0.30'; '0.40'; '0.50'; '0.60'; '0.70'];
% hold on
% scatter(x, y, 'filled');
% text(x+0.001, y + 0.02, c);
% hold on
% xlim([0, 1]);
% ylim([0, 1]);
% plot(x, y, 'DisplayName', 'Five');
% ylabel('True Positive');
% xlabel('False Positive');
% hold off
% % savefig('~/Downloads/personalRepo/Cov/Data3/ideal/roc_real_5');
% savefig('~/Downloads/personalRepo/Cov/Data3/ideal/cumulative_roc');
% %% histogram
% adjust = ones(N_neur * N_neur, 1);
% adjusty = ones(N_neur * N_neur, 1);
% for kk = 1:N_neur*N_neur
%     adjust(kk) = mean(real.rate_evt((ceil(kk / N_neur)), :));
%     thing = mod(kk, N_neur);
%     if (thing == 0)
%         thing = N_neur;
%     end
%     adjusty(kk) = mean(real.rate_evt(thing, :));
% end
%%
corr.corr_a_soma(isnan(corr.corr_a_soma)) = 2;
corr.corr_suite_soma(isnan(corr.corr_suite_soma)) = 2;
corr.corr_5_s(isnan(corr.corr_5_s)) = 2;
corr.corr_cas_soma(isnan(corr.corr_cas_soma)) = 2;

corr.corr_a_ideal(isnan(corr.corr_a_ideal)) = 2;
corr.corr_suite_ideal(isnan(corr.corr_suite_ideal)) = 2;
corr.corr_5_i(isnan(corr.corr_5_i)) = 2;
corr.corr_cas_ideal(isnan(corr.corr_cas_ideal)) = 2;
k_bins = 10;
[image, cluster_idx, image_order] = cluster(corr.corr_real, corr.corr_a_soma, N_neur, k_bins);
imagesc(image);
savefig('~/Downloads/personalRepo/Cov/Data3/soma/mls_corr_map');
[image, cluster_idx, image_order] = cluster(corr.corr_real, corr.corr_suite_soma, N_neur, k_bins);
imagesc(image);
savefig('~/Downloads/personalRepo/Cov/Data3/soma/suite_corr_map');
[image, cluster_idx, image_order] = cluster(corr.corr_real, corr.corr_5_s, N_neur, k_bins);
imagesc(image);
savefig('~/Downloads/personalRepo/Cov/Data3/soma/5_corr_map');
[image, cluster_idx, image_order] = cluster(corr.corr_real, corr.corr_cas_soma, N_neur, k_bins);
imagesc(image);
savefig('~/Downloads/personalRepo/Cov/Data3/soma/cas_corr_map');

[image, cluster_idx, image_order] = cluster(corr.corr_real, corr.corr_a_ideal, N_neur, k_bins);
imagesc(image);
savefig('~/Downloads/personalRepo/Cov/Data3/ideal/mls_corr_map');
[image, cluster_idx, image_order] = cluster(corr.corr_real, corr.corr_suite_ideal, N_neur, k_bins);
imagesc(image);
savefig('~/Downloads/personalRepo/Cov/Data3/ideal/suite_corr_map');
[image, cluster_idx, image_order] = cluster(corr.corr_real, corr.corr_5_i, N_neur, k_bins);
imagesc(image);
savefig('~/Downloads/personalRepo/Cov/Data3/ideal/5_corr_map');
[image, cluster_idx, image_order] = cluster(corr.corr_real, corr.corr_cas_ideal, N_neur, k_bins);
imagesc(image);
savefig('~/Downloads/personalRepo/Cov/Data3/ideal/cas_corr_map');
[rate, num, stats] = compare_img(real.rate_evt, cluster_idx, image_order);

%imagesc(image < 0.1 & image > -0.1);



% writematrix(image, '~/Downloads/personalRepo/Cov/Data3/ideal/mls_image.csv');
% writematrix(cluster_idx, '~/Downloads/personalRepo/Cov/Data3/ideal/mls_cluster_idx.csv');
% writematrix(image_order, '~/Downloads/personalRepo/Cov/Data3/ideal/mls_clusterorder.csv');
% writecell(rate, '~/Downloads/personalRepo/Cov/Data3/ideal/mls_seperate_spike.csv');
% writematrix(stats, '~/Downloads/personalRepo/Cov/Data3/ideal/mls_mean_std.csv');
% savefig(h, '~/Downloads/personalRepo/Cov/Data3/ideal/mls_fig.fig');

% writematrix(image, '~/Downloads/personalRepo/Cov/Data3/ideal/suite_image.csv');
% writematrix(cluster_idx, '~/Downloads/personalRepo/Cov/Data3/ideal/suite_cluster_idx.csv');
% writematrix(image_order, '~/Downloads/personalRepo/Cov/Data3/ideal/suite_clusterorder.csv');
% writecell(arate, '~/Downloads/personalRepo/Cov/Data3/ideal/suite_seperate_spike.csv');
% writematrix(stats, '~/Downloads/personalRepo/Cov/Data3/ideal/suite_mean_std.csv');
% savefig(h, '~/Downloads/personalRepo/Cov/Data3/ideal/suite_fig.fig');

% writematrix(image, '~/Downloads/personalRepo/Cov/Data3/ideal/five_image.csv');
% writematrix(cluster_idx, '~/Downloads/personalRepo/Cov/Data3/ideal/five_cluster_idx.csv');
% writematrix(image_order, '~/Downloads/personalRepo/Cov/Data3/ideal/five_clusterorder.csv');
% writecell(arate, '~/Downloads/personalRepo/Cov/Data3/ideal/five_seperate_spike.csv');
% writematrix(stats, '~/Downloads/personalRepo/Cov/Data3/ideal/five_mean_std.csv');
% savefig(h, '~/Downloads/personalRepo/Cov/Data3/ideal/five_fig.fig');

% writematrix(image, '~/Downloads/personalRepo/Cov/Data3/ideal/cas_image.csv');
% writematrix(cluster_idx, '~/Downloads/personalRepo/Cov/Data3/ideal/cas_cluster_idx.csv');
% writematrix(image_order, '~/Downloads/personalRepo/Cov/Data3/ideal/cas_clusterorder.csv');
% writecell(rate, '~/Downloads/personalRepo/Cov/Data3/ideal/cas_seperate_spike.csv');
% writematrix(stats, '~/Downloads/personalRepo/Cov/Data3/ideal/cas_mean_std.csv');
% savefig(h, '~/Downloads/personalRepo/Cov/Data3/ideal/cas_fig.fig');

%%
figure;

adjustz= reshape(diff_mls_a_s, N_neur * N_neur, 1);
scatter3(0, 0, 0);
for kk = 1:(N_neur * N_neur)
    [x y z] = sphere();
    hold on
    x = x/100 +adjust(kk);
    y = y/100 +adjusty(kk);
    z = z/100 + adjustz(kk);
    surf(x, y, z);
    hold off
end
%scatter3(adjust, adjusty, reshape(diff_mls_a_s, N_neur * N_neur, 1));
%%
figure;
%scatter(adjust, reshape(diff_mls_a_s, N_neur * N_neur, 1));
%scatter(adjust, reshape(diff_cas_s, N_neur * N_neur, 1));
%scatter(adjust, reshape(diff_suite_s, N_neur * N_neur, 1));
%scatter(adjust, reshape(diff_5_s, N_neur * N_neur, 1));
%scatter(adjust, reshape(diff_mls_a_i, N_neur * N_neur, 1));
%scatter(adjust, reshape(diff_cas_i, N_neur * N_neur, 1));
%scatter(adjust, reshape(diff_suite_i, N_neur * N_neur, 1));
%scatter(adjust, reshape(diff_5_i, N_neur * N_neur, 1));
%%
figure;
%histogram(reshape(norm_mls_a_s, N_neur * N_neur, 1), 'BinWidth', 0.005);
%histogram(reshape(norm_cas_s, N_neur * N_neur, 1), 'BinWidth', 0.005);
%histogram(reshape(norm_suite_s, N_neur * N_neur, 1), 'BinWidth', 0.005);
%histogram(reshape(norm_5_s, N_neur * N_neur, 1), 'BinWidth', 0.005);
%histogram(reshape(norm_mls_a_i, N_neur * N_neur, 1), 'BinWidth', 0.005);
%histogram(reshape(norm_cas_i, N_neur * N_neur, 1), 'BinWidth', 0.005);
%histogram(reshape(norm_suite_i, N_neur * N_neur, 1), 'BinWidth', 0.005);
%histogram(reshape(norm_5_i, N_neur * N_neur, 1), 'BinWidth', 0.005);

