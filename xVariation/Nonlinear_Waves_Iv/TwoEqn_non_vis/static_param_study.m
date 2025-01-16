function static_param_study
%     n_step = 50;
    min_lambda = load('Fr_tau0_min_lambda.txt'); %zeros(n_step);
    n_step = size(min_lambda,1);
    i_init = max(sum(min_lambda(:,1)~=0),1);
    j_init = 1; %sum(min_lambda(i_init,:)~=0)+1;
    if j_init > n_step
        i_init=i_init+1;
        j_init=1;
    end
    Fr_vals = linspace(0.6,5,n_step);
    tau0_vals = linspace(0,100,n_step);
    lambda_vals=linspace(50,1000,96);
    theta = 12;
    params_master = [Fr_vals(1),theta,tau0_vals(1),lambda_vals(1),false,1,true];
    [xi_master,y_master]=construct_depo_wave(true,params_master,false);
    fail_Fr=[];
    fail_tau0=[];
    fail_theta=[];
    for i=i_init:n_step
        tau0=tau0_vals(i);
        params_in = [Fr_vals(1),theta,tau0,lambda_vals(1),false,1,true];
        [xi_master,y_master]=construct_depo_wave(true,params_in,true,xi_master,y_master,params_master);
        params_master = params_in;
        if i == i_init
            j_min = j_init;
        else
            j_min = 1;
        end
        for j = j_min:n_step
            Fr = Fr_vals(j);
            if (j==j_min)
                params_loop = [Fr,theta,tau0,lambda_vals(1),false,1,true];
                [xi_loop,y_loop]=construct_depo_wave(true,params_loop,true,xi_master,y_master,params_master);
            else
                params_loop_new = [Fr,theta,tau0,lambda_vals(1),false,1,true];
                [xi_loop,y_loop]=construct_depo_wave(true,params_loop_new,true,xi_loop,y_loop,params_loop);
                params_loop = params_loop_new;
            end
            params_old = params_loop;
            xi_curr=xi_loop;
            y_curr=y_loop;
            for k=1:size(lambda_vals,2)
                lambda=lambda_vals(k);
                params_curr = [Fr,theta,tau0,lambda,false,1,true];
                try
                    [xi_curr,y_curr]=construct_depo_wave(true,params_curr,true,xi_curr,y_curr,params_old);
                 catch ME
                    try
                        lambda=lambda_vals(k+1);
                        params_curr = [Fr,theta,tau0,lambda,false,1,true];
                        [xi_curr,y_curr]=construct_depo_wave(true,params_curr,true,xi_curr,y_curr,params_old);
                        if (size(y_curr,1) == 6)
                            min_lambda(i,j) = lambda-10;
                            break
                        end
                    catch ME
                        fail_Fr(end+1) = Fr;
                        fail_theta(end+1) = theta;
                        fail_tau0(end+1) = tau0;
                    break
                end           
                params_old = params_curr;
                if (size(y_curr,1) == 6)
                    for l=1:round(lambda_vals(k)-lambda_vals(k-1))
                        lambda=lambda_vals(k)-l;
                        try
                            params_curr = [Fr,theta,tau0,lambda,false,1,true];
                            [xi_curr,y_curr]=construct_depo_wave(true,params_curr,true,xi_curr,y_curr,params_old);
                            params_old = params_curr;
                            if (size(y_curr,1) == 6)
                                min_lambda(i,j) = lambda+1;
                            break    
                            end
                            catch ME
                                min_lambda(i,j) = lambda+1;
                                break    
                        end
                    end
                elseif (k==size(lambda_vals,2))
                    min_lambda(i,j) = -1;
                end
            end
            if (mod(j,5)==0)
                fail_out = vertcat(fail_Fr,fail_theta,fail_tau0);
                save("Fr_tau0_fail.txt","fail_out","-ascii")
                save("Fr_tau0_min_lambda.txt","min_lambda","-ascii")
            end
        end
    end
    fail_out = vertcat(fail_Fr,fail_theta,fail_tau0);
    save("Fr_tau0_fail.txt","fail_out","-ascii")
    save("Fr_tau0_min_lambda.txt","min_lambda","-ascii")
end