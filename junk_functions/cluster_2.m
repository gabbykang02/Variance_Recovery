function [image, idx, A] = cluster_2(real, predict, N_neur, k_bins)
%CLUSTER Summary of this function goes here
%   Detailed explanation goes here
temp = real - predict;
[idx, c] = kmeans(temp, k_bins);

%sum up number of elements not retained
image_graph_sutff = temp < 0.3 & temp > -0.3;

true_neg = zeros(k_bins, 1);
false_neg = zeros(k_bins, 1);
true_pos = zeros(k_bins, 1);
false_pos = zeros(k_bins, 1);
num_in = zeros(k_bins, 1);
for k = 1:k_bins
    for kk = 1:N_neur
        if (idx(kk) == k)
            for second_kk = 1:N_neur
                true_neg(k) = true_neg(k) + (predict(kk, second_kk) > 0.3 | predict(kk, second_kk) < -0.3) & (real(kk, second_kk) <= 0.3 & real(kk, second_kk) >= -0.3);
            end
            for second_kk = 1:N_neur
                false_neg(k) = false_neg(k) - (predict(kk, second_kk) <= 0.3 & predict(kk, second_kk) >= -0.3) & (real(kk, second_kk) > 0.3 | real(kk, second_kk) < -0.3);
            end
            for second_kk = 1:N_neur
                true_pos(k) = true_pos(k) +  (predict(kk, second_kk) > 0.3 | predict(kk, second_kk) < -0.3) & (real(kk, second_kk) > 0.3 | real(kk, second_kk) < -0.3);
            end
            for second_kk = 1:N_neur
                false_pos(k) = false_pos(k) - (predict(kk, second_kk) <= 0.3 & predict(kk, second_kk) >= -0.3) & (real(kk, second_kk) <= 0.3 & real(kk, second_kk) >= -0.3);
            end
            num_in(k) = num_in(k) + 1;
        end
    end
%     for kk = 1:N_neur
%         num_in = 0;
%         if (idx(kk) == k && sum(image_graph_sutff(kk, :)) <= sum(idx == k))
%             not_retain(k) = not_retain(k) + 1;
%             num_in = num_in + 1;
%          end
%     end

%     num_in = 0;
%     for kk = 1:N_neur
%         if (idx(kk) == k)
%             num_in = num_in + 1;
%             
%             not_retain(k) = not_retain(k) - image_graph_sutff(kk, kk) + sum(image_graph_sutff(kk, :)) + sum(image_graph_sutff(:, kk));
%         end
%     end
%     not_retain(k) = ((N_neur*N_neur - num_in*num_in) - not_retain(k)) / num_in;
%not_retain(k) = not_retain(k) / num_in;
end
%[B A] = sort(not_retain, 'descend');
comp = (true_neg );%/num_in + true_pos/num_in;
[B A] = sort(comp);





count = 1;
% Change temp colors
% 400 = lime green
% 100 = light blue
% 0 = dark blue
% 600 = yellow

% nonzero correlation not retained, false neg
temp((predict <= 0.3 & predict >= -0.3) & (real > 0.3 | real < -0.3)) = 400;
% zero correlation not retained, true neg
temp((predict > 0.3 | predict < -0.3) & (real <= 0.3 & real >= -0.3)) = 0;
% zero correlation retained, false pos
temp((predict <= 0.3 & predict >= -0.3) & (real <= 0.3 & real >= -0.3)) = 100;
% nonzero correlation retained, true pos
temp((predict > 0.3 | predict < -0.3) & (real > 0.3 | real < -0.3)) = 600;

temp_adjust = zeros(N_neur);
% reorder to copy 1 first by columns
figure
for k = 1:k_bins
    for kk = 1:N_neur
        if (idx(kk) == A(k))
            temp_adjust(count, :) = temp(kk, 1:N_neur);
            count = count + 1;
        end
    end
end
count = 1;
image = temp_adjust;
%copy by rows
for k = 1:k_bins
    for kk = 1:N_neur
        if (idx(kk) == A(k))
            image(:, count) = temp_adjust(1:N_neur, kk);
            count = count + 1;
        end
    end
end
end

