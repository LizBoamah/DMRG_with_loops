function [left_env, middle_env,right_env] = contract_environments(M, M_, H, edge, N)
% middle_env = [];
    % Special case where edge is between the last and first sites
    if edge(1) == N && edge(2) == 1
        % Initialize left and right environments to identity
        left_env = 1;
        right_env = 1;
        % 
         middle_env = contract_middle_environment(M,H,M_,N);

    elseif edge(1) == 1
        % Contract environment from the right
        right_env = 1;
        
        
        for j = N:-1:edge(2)+1
            right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
        end
        left_env = 1;
        middle_env = 1;
    elseif edge(1) == N-1
        % Contract environment from the left
        left_env = 1;
        for j = 1:edge(1)-1
            left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
        end
        right_env = 1;
        middle_env = 1;
    else
        % Contract environment from the left
        left_env = 1;
        for j = 1:edge(1)-1
            left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
        end
        
        % Contract environment from the right
        right_env = 1;
        for j =  N:-1:edge(2)+1
            right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
        end
        middle_env =1;
    end
end
