function [s,i] = to_spike_rate(soma,ideal, N_neur, bin_size, frame_size)
%%
% ~40ms, smaller later
%%
s = zeros([N_neur floor(600 / (floor(frame_size)))]);
i = zeros([N_neur floor(600 / (floor(frame_size)))]);

for kk = 1:N_neur
    for jj = 1:(floor(frame_size)):(size(soma, 2) - (floor(frame_size)))
        s(kk, floor(jj / (floor(frame_size)))) = sum(~isnan(soma(kk, jj:(jj+(floor(frame_size)) - 1)))) / bin_size;
        i(kk, floor(jj / (floor(frame_size)))) = sum(~isnan(ideal(kk, jj:(jj+(floor(frame_size)) - 1)))) / bin_size;
    end
end
end

