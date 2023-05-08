function erosion_point = get_erosion_time(u_p,tol)
% Gets the index of the time which corresponds to the point of deposition
% across the domain.
    n_points = size(u_p,2);
    erosion_point = sum(abs(u_p)<tol,2);
    erosion_point = erosion_point/n_points;
    end
