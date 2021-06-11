function new_arr = mk_cov(arr)
%MK_COV Summary of this function goes here
%   Detailed explanation goes here
new_arr = zeros(size(arr));
new_arr(:, 1) = arr(:, 1);
arr = arr.';
for kk = 2:size(arr,2)
    A = [ zeros(kk - 1, 1).' arr(1, 1:(size(arr, 2) - kk + 1))];

    for jj = 1:(kk-1) %(size(arr, 2)):-1:(kk + 1)
        A(jj) = arr(size(arr, 2) - kk+1+jj);
        %A(jj) = arr(size(arr, 2) - jj+1);
    end
    new_arr(:, kk) = A.'; 
end
end

