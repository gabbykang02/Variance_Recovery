function [s,i] = mls_to_spike_rate(soma_all,ideal_all, N_neur, s_count, i_count, frame_size, bin_size)
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts MLSpike predictions into spike rates
% 
% INPUTS:
% soma_all = MLspike predictions based on soma data, non_zero_soma_activity neurons by nt
% ideal_all = MLspike predictions based on idealTraces, non_zero_idealTrace_activity neurons by nt
% N_neur = total number of neurons
% s_count = N_neur x 1 matrix with indices of all non_zero_soma_activity
% i_count = N_neur x 1 matrix with indices of all non_zero_idealTrace_activity
% frame_size = # of frames in a bin
% bin_size = # of seconds in a bin
% 
% OUTPUTS:
% s = MLSpike soma prediction spike rates
% i = MLSpike idealTraces prediction spike rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
s = zeros([N_neur floor(600 / frame_size)]);
i = zeros([N_neur floor(600 / frame_size)]);

for curr_neur = 1:N_neur
    for curr_frame = 1:frame_size:(600 - frame_size)
        if (s_count(curr_neur) ~= 0)
            % temp is the curr_neur row from the MLSpike prediction matrix
            temp = soma_all(s_count(curr_neur), :);
            % Counts the number of events in temp that occur within the bin
            s(curr_neur, floor(curr_frame / frame_size)) = sum((temp > ((curr_frame - frame_size) / 30)) & (temp < (curr_frame /30))) / bin_size;
        end
        if (i_count(curr_neur) ~= 0)
            temp = ideal_all(i_count(curr_neur), :);
            i(curr_neur, floor(curr_frame / frame_size)) = sum((temp > ((curr_frame - frame_size) / 30)) & (temp < (curr_frame /30))) / bin_size;
           
        end
    end
end
end

