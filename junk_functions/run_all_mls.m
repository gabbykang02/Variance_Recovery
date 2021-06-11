function run_all_mls(soma,ideal)
dt = 1/30;
%soma = soma.';
ideal = ideal.';
par = tps_mlspikes('par');
% (indicate the frame duration of the data)
par.dt = dt;
% (set physiological parameters)
par.a = 0.07; % DF/F for one spike
par.tau = 1; % decay time constant (second)
par.saturation = 0.1; % OGB dye saturation
% (set noise parameters)
par.finetune.sigma = .02; % a priori level of noise (if par.finetune.sigma
                          % is left empty, MLspike has a low-level routine 
                          % to try estimating it from the data
par.drift.parameter = .01; % if par.drift parameter is not set, the 
                           % algorithm assumes that the baseline remains
                           % flat; it is also possible to tell the
                           % algorithm the value of the baseline by setting
                           % par.F0
% (do not display graph summary)
par.dographsummary = false;

par.a = 0.03;
%[spikest fit drift] = spk_est(soma,par);
[idealEst fit drift] = spk_est(ideal, par);
temp1 = idealEst.';
%temp2 = spikest.';
writecell(temp1, '~/Downloads/personalRepo/Cov/Data/mls_0.03_idealdenoise.csv');
%writecell(temp2, '~/Downloads/personalRepo/Cov/Data/mls_0.03_soma.csv');

alpha_avg = robustSTD(soma);
alpha_avg = mean(alpha_avg);
par.a = alpha_avg;
%[spikest fit drift] = spk_est(soma,par);
alpha_avg = robustSTD(ideal);
alpha_avg = mean(alpha_avg);
par.a = alpha_avg;
[idealEst fit drift] = spk_est(ideal, par);
temp1 = idealEst.';
%temp2 = spikest.';

%writematrix(alpha_avg, '~/Downloads/personalRepo/Cov/Data/mls_avg.csv');
writecell(temp1, '~/Downloads/personalRepo/Cov/Data/mls_a_idealdenoise.csv');
%writecell(temp2, '~/Downloads/personalRepo/Cov/Data/mls_a_soma.csv');

par.a = 0.07;
%[spikest fit drift] = spk_est(soma,par);
[idealEst fit drift] = spk_est(ideal, par);
temp1 = idealEst.';
%temp2 = spikest.';
writecell(temp1, '~/Downloads/personalRepo/Cov/Data/mls_0.07_idealdenoise.csv');
%writecell(temp2, '~/Downloads/personalRepo/Cov/Data/mls_0.07_soma.csv');
end

