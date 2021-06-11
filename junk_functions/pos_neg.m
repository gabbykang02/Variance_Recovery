function [roc_count, roc_countideal] = pos_neg(soma, ideal, evt, threshold)
roc_count = zeros(size(soma, 1), 3);
roc_countideal = zeros(size(soma, 1), 3);
for kk = 1:size(soma, 1)
    for jj = 1:sum(~isnan(evt{1, kk}))
        % second column = total number of real events
        roc_count(kk, 2) = sum(~isnan(evt{1, kk}));
        roc_countideal(kk, 2) = sum(~isnan(evt{1, kk}));
        % third column = totla number of predicted events
        roc_count(kk, 3) = sum(~isnan(soma(kk, :)));
        roc_countideal(kk, 3) = sum(~isnan(ideal(kk, :)));
        % first column = total number of correectly identified events
        temp = find(soma(kk, :) < evt{1, kk}(1, jj) + threshold ...
            & soma(kk, :) > evt{1, kk}(1, jj) - threshold, 1);
        temp2 = find(ideal(kk, :) < evt{1, kk}(1, jj) + threshold ...
            & ideal(kk, :) > evt{1, kk}(1, jj) - threshold, 1);
        if ~isempty(temp)
            roc_count(kk, 1) = roc_count(kk, 1) + 1;
        end
        if ~isempty(temp2)
            roc_countideal(kk, 1) = roc_countideal(kk, 1) + 1;
        end
    end
end
end

