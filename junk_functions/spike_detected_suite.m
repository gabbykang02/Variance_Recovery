function [soma,ideal] = spike_detected_suite(real_evt, soma_evt, ideal_evt)
soma = struct();
ideal = struct();
real_evt = real_evt(1:338, :);
temp = real_evt.';
soma.total = sum(~isnan(temp)).';
ideal.total = soma.total;

soma.found = zeros(size(real_evt, 1), 1);
ideal.found = zeros(size(real_evt, 1), 1);

soma.acc_rate = zeros(size(real_evt, 1), 1);
ideal.acc_rate = zeros(size(real_evt, 1), 1);
for n_neur = 1:size(real_evt)
    evt_idx = 1;
    for kk = 1:600   
        real = real_evt(n_neur, evt_idx);
        if (~(isnan(real)) && ~isnan(soma_evt(n_neur, kk)) &&  kk/30 > real - 0.3 && kk/30 < real + 0.3)
            evt_idx = evt_idx + 1;
            soma.found(n_neur) = soma.found(n_neur) + 1;
        end
        if (~isnan(real) && kk/30 > real) 
            while(kk / 30 > real)
                evt_idx = evt_idx + 1;
                real = real_evt(n_neur, evt_idx);
            end
        end
    end
    soma.acc_rate(n_neur) = soma.found(n_neur)/soma.total(n_neur);
    
    evt_idx = 1;
    for kk = 1:600   
        real = real_evt(n_neur, evt_idx);
        if (~(isnan(real)) && ~isnan(ideal_evt(n_neur, kk)) &&  kk/30 > real - 0.3 && kk/30 < real + 0.3)
            evt_idx = evt_idx + 1;
            ideal.found(n_neur) = ideal.found(n_neur) + 1;
        end
        if (~isnan(real) && kk/30 > real) 
            while(kk / 30 > real)
                evt_idx = evt_idx + 1;
                real = real_evt(n_neur, evt_idx);
            end
        end
    end
    ideal.acc_rate(n_neur) = ideal.found(n_neur)/ideal.total(n_neur);
end



