function fields = real_roc(corr_real, corr_method, threshold)
fields = struct();
corr_real = abs(corr_real);
for kk = 1:338
    for jj = 1:338
        if (jj <= kk) 
            corr_real(kk, jj) = NaN;
        end
    end
end    
corr_method = abs(corr_method);
% true positive
fields.num_correlated = sum(sum(corr_real > 0.3));
fields.true_positive = sum(sum(corr_method > threshold & corr_real > 0.3));
fields.true = fields.true_positive/fields.num_correlated;

% false positive
fields.non_correlated = sum(sum(corr_real <= 0.3));
fields.false_positive = sum(sum(corr_method > threshold & corr_real <= 0.3));
fields.zero = fields.false_positive/fields.non_correlated;

end

