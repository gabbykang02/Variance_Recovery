function [soma_rate, ideal_rate] = prob_to_spike_rate(soma_all,ideal_all, N_neur, bin_size, frame_size)
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts probability predictions into spike rates
% 
% INPUTS:
% soma_all = predictions based on soma data, N_neur x nt
% ideal_all = predictions based on idealTraces, N_neur x nt
% N_neur = total number of neurons
% frame_size = # of frames in a bin
% bin_size = # of seconds in a bin
% 
% OUTPUTS:
% s = soma prediction spike rates
% i = idealTraces prediction spike rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
soma_rate = zeros([N_neur floor(600 / frame_size)]);
ideal_rate = zeros([N_neur floor(600 / frame_size)]);

for curr_neur = 1:N_neur
    for curr_frame = 1:frame_size:(600-frame_size)
        soma_rate(curr_neur, floor(curr_frame / frame_size)) = 1/bin_size * trapz(1/30*([curr_frame:(curr_frame+frame_size)]), soma_all(curr_neur, curr_frame:(curr_frame+frame_size)));
        ideal_rate(curr_neur, floor(curr_frame / frame_size)) = 1/bin_size *trapz(1/30*([curr_frame:(curr_frame+frame_size)]), ideal_all(curr_neur, curr_frame:(curr_frame+frame_size)));
    end
end

end

