N_neur = 383;

%% MLSpike
N_neur = 383;
soma = readmatrix('~/Downloads/personalRepo/Cov/Data/soma.csv');
ideal = readmatrix('~/Downloads/personalRepo/Cov/Data/denoiseIdeal.csv');
[s_zero,i_zero] = non_zero_trace(soma, ideal, N_neur);
non_soma = zeros(sum(s_zero ~= 0), 600);
non_ideal = zeros(sum(i_zero ~= 0), 600);
for kk = 1:N_neur
    if (s_zero(kk) ~= 0)
        non_soma(s_zero(kk), :) = soma(kk, :);
    end
    if (i_zero(kk) ~= 0)
        non_ideal((i_zero(kk)), :) = ideal(kk, :);
    end
end

run_all_mls(non_soma,non_ideal);




