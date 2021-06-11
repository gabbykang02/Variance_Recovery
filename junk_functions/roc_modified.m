function fields = roc_modified(corr_real, corr_method, threshold)
%ROC_MODIFIED Summary of this function goes here
%   Detailed explanation goes here
fields = struct();
image = corr_real - corr_method;
image = image < threshold & image > -threshold;

temp = corr_real > 0.1 | corr_real < -0.1;
correlated_in_method = corr_method > 0.1 | corr_method < -0.1;
% how many were detected and important, true positive and relevant
fields.true_det = sum(sum(temp & image));
fields.total_true = sum(temp(:));
fields.true = fields.true_det/fields.total_true;

% how many were detected and not importnat
fields.zero_corr_det = sum(sum(~temp & image));
fields.total_zero_corr = sum(~temp(:));
% Neurons uncorrelated, but detected as having correlation
fields.total_zero_det = sum(sum(~temp));
fields.zero_corr_det = sum(image & (correlated_in_method));
% correlated and not detected
% fields.zero_corr_det = sum(sum(temp & ~image));
% fields.total_zero_corr = sum(temp(:));
fields.zero = fields.zero_corr_det/fields.total_zero_corr;

end

