function multi_wave_create
    lambda = 50;
%     tau0 = 0;
    theta = 10;
    npts = 50;
    h_crit_mat = zeros(npts,npts);
    theta_start = 0;
    theta_stop = 100;
    theta_list = linspace(theta_start,theta_stop,npts);
    Fr_start = 0.6;
    Fr_stop = 5;
    Fr_list = linspace(Fr_start,Fr_stop,npts);
    i_fail = [];
    j_fail = [];
    Fr_fail = [];
    theta_fail = [];
    for i=1:npts
        Fr = Fr_list(i);
        for j=1:npts
            tau0 = theta_list(j);
            params = [Fr,theta,tau0,lambda,false,1.0,true];
            if (j == 1)
                if (i==1)
                    [col_xi,col_y] = construct_depo_wave(true,params,false);
                else
                    try
                        [col_xi,col_y] = construct_depo_wave(true,params,true,col_xi,col_y,col_params);
                        catch ME
                            switch ME.identifier
                                case 'Daniel:MaxIt'
                                    try
                                        [col_xi,col_y] = construct_depo_wave(true,params,false);
                                        catch ME
                                            switch ME.identifier
                                                case 'Daniel:MaxIt'
                                                    i_fail(end+1) = i;
                                                    j_fail(end+1) = j;
                                                    Fr_fail(end+1) = Fr;
                                                    theta_fail(end+1) = theta;
                                                otherwise
                                                    rethrow(ME)
                                            end
                                    end
                                otherwise
                                    rethrow(ME)
                            end
                    end
                end
                col_params = params;
                curr_params = params;
                curr_xi = col_xi;
                curr_y = col_y;
            else
            try
                [curr_xi,curr_y] = construct_depo_wave(true,params,true,curr_xi,curr_y,curr_params);
                catch ME
                    switch ME.identifier
                        case 'Daniel:MaxIt'
                            try
                                [curr_xi,curr_y] = construct_depo_wave(true,params,false);
                                catch ME
                                    switch ME.identifier
                                        case 'Daniel:MaxIt'
                                            i_fail(end+1) = i;
                                            j_fail(end+1) = j;
                                            Fr_fail(end+1) = Fr;
                                            theta_fail(end+1) = theta;
                                        otherwise
                                            rethrow(ME)
                                    end
                            end
                        otherwise
                            rethrow(ME)
                    end
            end
                curr_params = params;
            end
%             [curr_xi,curr_y] = construct_depo_wave(true,params,false);
            crit_h_rel = (curr_y(2,1)*Fr)^(2/3);
            h_crit_mat(i,j) = crit_h_rel;
        end
    end
    save("Results/multi_comp_tau0_Fr.txt","h_crit_mat","-ascii")
    fail_out = horzcat(i_fail,j_fail,Fr_fail,theta_fail);
    save("Results/fail_cond.txt","fail_out","-ascii")
end