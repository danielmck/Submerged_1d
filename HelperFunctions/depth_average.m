function depth_ave = depth_average(X,N,num_times)
    dz = 1/(N-0.5);
    depth_ave= zeros(num_times,1);
    for j=1:num_times
        depth_ave(j) = (sum(X(:,j))-(X(1,j)+X(N,j)/2))*dz/(1-dz);
    end
end
