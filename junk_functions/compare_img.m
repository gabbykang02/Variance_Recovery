function [rate,num_per_cluster, mean_std_stats] = compare_img(real, cluster_idx, order)
num_per_cluster = zeros(size(order, 1), 1);
for kk = 1:size(num_per_cluster, 1)
    num_per_cluster(kk) = sum(cluster_idx == order(kk));
end
rate = cell(size(order, 1), 1);

%for 1:num bins
for kk = 1:size(order, 1)
    rate{kk, 1} = zeros(num_per_cluster(kk), 1);
    idx = 1;
    %for 1:N_neur
    for jj = 1:size(cluster_idx, 1)
        if (cluster_idx(jj) == order(kk))
            rate{kk, 1}(idx) = mean(real(jj, :));
            idx = idx + 1;
        end
    end
end

exp1 = zeros(kk, 1);
exp2 = zeros(kk, 1);
exp3 = zeros(kk, 1);
exp4 = zeros(kk, 1);
for kk = 1:size(order, 1)
    exp1(kk) = mean(rate{kk, 1}(:));
    exp2(kk) = std(rate{kk, 1}(:));
    exp3(kk) = max(rate{kk, 1}(:));
    exp4 (kk)= min(rate{kk, 1}(:));
end
mean_std_stats = [exp1 exp2 exp3 exp4 num_per_cluster];
end

