function [s,i] = non_zero_trace(soma, ideal, N_neur)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In MLSpike, Neurons with 0 detected activity had to be removed. The
% function, non_zero_trace, returns a N_neur size matrix indicating what
% index in the MLSpike predictions corresponds to the soma neurons.
% 
%   EXAMPLE: if s(10) = 8, the predictions at ml_spike_soma_predictions(8)
%            corresponds with the real_soma_activity(10)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
s = zeros(N_neur, 1);
i = zeros(N_neur, 1);

ideal = ideal(1:N_neur,:);
ideal = ideal.';
ideal = mat2cell(ideal, [600], ones(N_neur, 1));

soma = soma(1:N_neur, :);

i_idx = 1;
s_idx = 1;
for kk = 1:N_neur
    if (sum(ideal{1, kk}) ~= 0)
        i(kk) = i_idx;
        i_idx = i_idx + 1;
    end
    if (sum(soma(kk, :)) ~= 0)
        s(kk) = s_idx;
        s_idx = s_idx + 1;
    end
end

end

