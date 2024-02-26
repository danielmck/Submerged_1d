function multi_lambda_create
%     lambda = 50;
    tau0 = 0;
    theta = 10;
    npts = 50;
    lambda_list = linspace(10,500,npts); %[10,20,50,75,100,120,150,200,250,300,400,500];
    Fr = 1;
    uw_list = zeros(size(lambda_list));

    h_fig = figure();
    SetPaperSize(10,10);
    hold on
%     u_fig = figure();
%     SetPaperSize(15,7);
    hold on
    
    for j=1:size(lambda_list,2)
        lambda = lambda_list(j);
        curr_params = [Fr,theta,tau0,lambda,false,1.0,true];
        if (j == 1)
            [curr_xi,curr_y] = construct_depo_wave(true,curr_params,false);
        else
            try
                [curr_xi,curr_y] = construct_depo_wave(true,curr_params,true,last_xi,last_y,last_params);
                catch ME
                    switch ME.identifier
                        case 'Daniel:MaxIt'
                            try
                                [curr_xi,curr_y] = construct_depo_wave(true,curr_params,false);
                                    catch ME
                                        switch ME.identifier
                                            case 'Daniel:MaxIt'
                                                error("Non convergence for $\lambda=$"+num2str(lambda))
                                            otherwise
                                                rethrow(ME)
                                        end
                             end
                                otherwise
                                    rethrow(ME)
                end
            end
        end
        last_params = curr_params;
        last_xi = curr_xi;
        last_y = curr_y;
        uw_list(j) = curr_y(1,1);
%         figure(h_fig)
%         plot(lambda*curr_xi,curr_y(4,:), 'color', 'k')
%         figure(u_fig)
%         plot(lambda*curr_xi,-curr_y(2,:)./curr_y(4,:)+curr_y(1,:), 'color', 'k')
    end
    figure(h_fig)
    plot(lambda_list,uw_list, 'color', 'k')
    xlabel("$\lambda$")
    ylabel("$u_w$")
%     xticks([0,50,100,150,200,250,300,350,400,450,500])
%     yticks([0.5,1,1.5,2.0,2.5,3])
    xlim([10,500])
    xticks([10,100,200,300,400,500])
    title("$\theta="+num2str(theta)+"^{\circ}$, $Fr="+num2str(Fr)+"$, $\tau_0 ="+num2str(tau0)+"$Pa")
    exp_graph(h_fig,"Two_eqn_lambda_comp_uw.pdf")
%     figure(u_fig)
%     xlabel("$\xi$")
%     ylabel("$u$")
%     xticks([0,50,100,150,200,250,300,350,400,450,500])
%     yticks([0,0.5,1,1.5,2.0,2.5])
%     exp_graph(u_fig,"Two_eqn_lambda_comp_u.pdf")
end