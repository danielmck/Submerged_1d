function crit_press = get_critical_value(u_p,p_e,erosion,tol)
% Finds the critical value of a field at the point at which the flow
% deposits
    n_points = size(u_p,2);
    n_times = size(u_p,1);
    flow_points = sum(abs(u_p)<tol,1);
    crit_press = zeros(n_points,1);
    if (~erosion)
        flow_points = n_times-flow_points+1;
    end
    flow_points(flow_points>n_times)=n_times;
    if (erosion)
        for i=2:n_points
            if ((abs(u_p(flow_points(i),i))<tol)&&(abs(u_p(flow_points(i)-1,i))>tol))
                break
            end
            while (abs(u_p(flow_points(i),i))<tol && flow_points(i)>1)
                flow_points(i) = flow_points(i)-1;
            end
            while (abs(u_p(flow_points(i)-1,i))>tol && flow_points(i)<n_times)
                flow_points(i) = flow_points(i)+1;
            end
        end
    else
        for i=2:n_points
            if ((abs(u_p(flow_points(i),i))>tol)&&(abs(u_p(flow_points(i)-1,i))<tol))
                break
            end
            while (abs(u_p(flow_points(i),i))>tol && flow_points(i)<n_times)
                flow_points(i) = flow_points(i)+1;
            end
            while (abs(u_p(flow_points(i)-1,i))<tol && flow_points(i)>1)
                flow_points(i) = flow_points(i)-1;
            end
        end
    end
    for i=1:n_points
        crit_press(i,1) = p_e(flow_points(i),i);
    end
end