function [image, idx, A] = cluster(real, predict, N_neur, k_bins)
%CLUSTER Summary of this function goes here
%   Detailed explanation goes here
temp = real - predict;
[idx, c] = kmeans(temp, k_bins);

%sum up number of elements not retained
image_graph_sutff = temp < 0.1 & temp > -0.1;

not_retain = zeros(k_bins, 1);
for k = 1:k_bins
    for kk = 1:N_neur
        num_in = 0;
        if (idx(kk) == k)
            for second_kk = 1:N_neur
                not_retain(k) = not_retain(k) + image_graph_sutff(kk, second_kk);
            end
            num_in = num_in + 1;
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
[B A] = sort(not_retain);





count = 1;
% Change temp colors
% 400 = lime green
% 100 = light blue
% 0 = dark blue
% 600 = yellow

% nonzero correlation not retained
temp((predict <= 0.1 & predict >= -0.1) & (real > 0.1 | real < -0.1)) = 400;
% zero correlation not retained
temp((predict > 0.1 | predict < -0.1) & (real <= 0.1 & real >= -0.1)) = 0;
% zero correlation retained
temp((predict <= 0.1 & predict >= -0.1) & (real <= 0.1 & real >= -0.1)) = 100;
% nonzero correlation retained
temp((predict > 0.1 | predict < -0.1) & (real > 0.1 | real < -0.1)) = 600;

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

