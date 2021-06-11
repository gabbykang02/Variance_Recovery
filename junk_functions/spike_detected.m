function [soma,ideal] = spike_detected(real_evt, soma_evt, ideal_evt, soma_all, ideal_all)
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
num_offset_s = 0;
num_offset_i = 0;
for n_neur = 1:338%size(real_evt)
    evt_idx = 1;

    for kk = 1:size(soma_evt, 2)
        while (sum(sum(soma_all(n_neur, :))) == 0) 
            n_neur = n_neur + 1;
            num_offset_s = num_offset_s + 1;
        end
        real = real_evt(n_neur, evt_idx);
        if ( ~(isnan(real)) && soma_evt(n_neur - num_offset_s, kk) > real - 0.3 && soma_evt(n_neur - num_offset_s, kk) < real + 0.3)
            evt_idx = evt_idx + 1;
            soma.found(n_neur) = soma.found(n_neur) + 1;
        end
        if (~isnan(real) && soma_evt(n_neur - num_offset_s, kk) > real) 
            while(soma_evt(n_neur - num_offset_s, kk) > real)
                evt_idx = evt_idx + 1;
                real = real_evt(n_neur, evt_idx);
            end
        end
    end

    soma.acc_rate(n_neur) = soma.found(n_neur)/soma.total(n_neur);
end

for n_neur = 1:338%size(real_evt)
    evt_idx = 1;
    for kk = 1:size(ideal_evt, 2)
        while (sum(sum(ideal_all(n_neur, :))) == 0) 
            n_neur = n_neur + 1;
            num_offset_i = num_offset_i + 1;
        end
        real = real_evt(n_neur, evt_idx);
        if (~(isnan(real)) && ideal_evt(n_neur - num_offset_i, kk) > real - 0.3 && ideal_evt(n_neur - num_offset_i, kk) < real + 0.3)
            evt_idx = evt_idx + 1;
            ideal.found(n_neur) = ideal.found(n_neur) + 1;
        end
        if (~isnan(real) && ideal_evt(n_neur - num_offset_i, kk) > real) 
            while(ideal_evt(n_neur - num_offset_i, kk) > real)
                evt_idx = evt_idx + 1;
                real = real_evt(n_neur, evt_idx);
            end
        end
    end
    ideal.acc_rate(n_neur) = ideal.found(n_neur)/ideal.total(n_neur);
end



