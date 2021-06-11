function [soma_rate, ideal_rate] = deconv_to_spike_rate(soma_all,ideal_all, N_neur, bin_size, frame_size)
soma_rate = zeros([N_neur floor(600 / frame_size)]);
ideal_rate = zeros([N_neur floor(600 / frame_size)]);
%scatter(1/30*[1:600], 0.15*soma_max(1, :));
for kk = 1:N_neur
    soma_all(soma_all == 0) = NaN;
    ideal_all(ideal_all == 0) = NaN;
    min_s = min(soma_all(kk, :));
    min_i = min(ideal_all(kk, :));
    
%     % remove baselinee by subtracting minimum
%     soma_all(kk, :) = soma_all(kk, :) - min_s;
%     % recalculate new minimum
%     temp = soma_all(kk, :);
%     temp(temp == 0) = NaN;
%     soma_all(kk, :) = temp;
%     min_s = min(soma_all(kk, :));
%     % repeat for ideal
%     % remove baselinee by subtracting minimum
%     ideal_all(kk, :) = ideal_all(kk, :) - min_i;
%     % recalculate new minimum
%     temp = ideal_all(kk, :);
%     temp(temp == 0) = NaN;
%     ideal_all(kk, :) = temp;
%     min_i = min(ideal_all(kk, :));
    for jj = 1:frame_size:(600-frame_size)
        
        soma_rate(kk, floor(jj / frame_size)) = 1/bin_size * soma_all(kk, jj) / min_s;
        ideal_rate(kk, floor(jj / frame_size)) = 1/bin_size * ideal_all(kk, jj) / min_i;
        %ideal_rate(kk, floor(jj / frame_size)) = 1/bin_size *trapz(1/30*([jj:(jj+frame_size)]), ideal_all(kk, jj:(jj+frame_size)));
    end
end
soma_rate(isnan(soma_rate)) = 0;
ideal_rate(isnan(ideal_rate)) = 0;
end

