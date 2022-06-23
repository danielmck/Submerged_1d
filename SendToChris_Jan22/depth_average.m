% Function to depth average a variable X, with N points at num_times
% different time values. Returns a 1D matrix of length num_times. N points
% can either be a single value or a matrix of num_times different values
function depth_ave = depth_average(X,N,num_times)
    depth_ave= zeros(num_times,1);
    if (length(N) == num_times)  
        for j=1:num_times
            depth_ave(j) = (sum(X(:,j))-(X(1,j)+X(N(j),j)/2))/N;
        end
    elseif (length(N) == 1)
        for j=1:num_times
            depth_ave(j) = (sum(X(:,j))-(X(1,j)+X(N,j)/2))/N;
        end
    else
        error('N must have length 1 or length num_times');
    end
end
