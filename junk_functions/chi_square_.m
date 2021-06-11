spike_rate_thresh = 0.5;
total = 383;
counts = zeros([4 10]);
r_counts = zeros(10, 1);
c_counts = zeros(4, 1);
expected = zeros(4, 10);

num_set = "2";
type = "ideal";
method = "cas";
if num_set == "1"
    rate = readmatrix('Data/' + type + '/' + method + '_seperate_spike.csv');
else
    rate = readmatrix('Data' + num_set + '/' + type + '/' + method + '_seperate_spike.csv');
end

df = 9*3;
for k = 1:10
    r_counts(k) = sum(~isnan(rate(k, :)));
    for j = 1:4
        counts(j, k) = sum(rate(k, :) < j * 0.5 & rate(k, :) > ((j-1)* 0.2));
        c_counts(j) = c_counts(j) + counts(j, k);
    end
end
for k = 1:10
    for j = 1:4
            expected(j, k) = r_counts(k) * c_counts(j) / total;
    end
end

%chi swuare analysis
chi = 0;
for k = 1:10
    for j = 1:4
        if (expected(j, k) ~= 0)
            chi = chi + ((counts(j, k) - expected(j, k))^2 /  expected(j, k));
        end
    end
end
p = chi2cdf(chi, df);
