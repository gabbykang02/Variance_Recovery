function [soma_comp, ideal_comp] = cov_diff(real_corr, soma_corr, ideal_corr, N_neur)
soma_comp = zeros(N_neur);
ideal_comp = zeros(N_neur);
for kk = 1:N_neur
    for jj = 1: N_neur
        soma_comp(kk, jj) = real_corr(kk, jj) - soma_corr(kk, jj);
        ideal_comp(kk, jj) = real_corr(kk, jj) - ideal_corr(kk, jj);
        %ideal_comp(kk, jj) = norm([real_corr(kk, jj); ideal_comp(kk, jj)]);
    end
end
end

