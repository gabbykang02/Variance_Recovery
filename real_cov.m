function [co] = real_cov(or)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the correlation between neuron spike rates. 
% 
% Takes in a N_neur x nt matrix, outputs a N_neur x N_neur matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
or = or.';
co = corrcoef(or);

end

