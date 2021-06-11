function fields = real_roc_id(corr_real, corr_method, threshold, id)
fields = struct();
corr_real = abs(corr_real);
for kk = 1:338
    for jj = 1:338
        if (jj <= kk) 
            corr_real(kk, jj) = NaN;
        end
    end
end
temp_real1= zeros(338);
temp_method1 = zeros(338);
idx = 1;
for kk = 1:size(id, 1)
    temp_real1(idx, :) = corr_real(id(idx), :);
    temp_method1(idx, :) = corr_method(id(idx), :);
    idx = idx + 1;
end
temp_real1 = temp_real1(1:size(id, 1), :);
temp_method1 = temp_method1(1:size(id, 1), :);
temp_real= zeros(size(id, 1), 338);
temp_method = zeros(size(id, 1), 338);
idx = 1;
for kk = 1:size(id, 1)
    temp_real(:, idx) = temp_real1(:, id(idx));
    temp_method(:, idx) = temp_method1(:, id(idx));
    idx = idx + 1;
end
temp_real = temp_real(:, 1:size(id, 1));
temp_method = temp_method(:, 1:size(id, 1));
corr_method = abs(corr_method);
% true positive
fields.num_correlated = sum(sum(temp_real > 0.2));
fields.true_positive = sum(sum(temp_method > threshold & temp_real > 0.2));
fields.true = fields.true_positive/fields.num_correlated;

% false positive
fields.non_correlated = sum(sum(temp_real <= 0.2));
fields.false_positive = sum(sum(temp_method > threshold & temp_real <= 0.2));
fields.zero = fields.false_positive/fields.non_correlated;

end

