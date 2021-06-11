function [soma_comp, ideal_comp] = cov_compare(real_corr, soma_corr, ideal_corr, N_neur, real_rate, soma_rate, ideal_rate)
soma_comp = zeros(N_neur);
ideal_comp = zeros(N_neur);
for kk = 1:N_neur
    for jj = 1: N_neur
        %(x, y) coordinates of kk's real rate vs correlation between kk and jj
        a = [mean(real_rate(kk, :)) real_corr(kk, jj)];
        b = [mean(soma_rate(kk, :)) soma_corr(kk, jj)];
        c = [mean(ideal_rate(kk, :)) ideal_corr(kk, jj)];
        soma_comp(kk, jj) = norm( b - a );
        ideal_comp(kk, jj) = norm(c - a);
        %soma_comp(kk, jj) = norm([real_corr(kk, jj); soma_corr(kk, jj)]);
        %ideal_comp(kk, jj) = norm([real_corr(kk, jj); ideal_corr(kk, jj)]);
    end
end
end

