function [newleft_env, newmiddle_env,newright_env] = new_contract_environments(M, M_, H, edge, N)
% middle_env = [];
    % Special case where edge is between the last and first sites
    if edge(1) == N && edge(2) == 1
        % Initialize left and right environments to identity
        newleft_env = 1;
        newright_env = 1;
        % 
        newmiddle_env = newcontract_middle_environment(M,H,M_,N);
        
    elseif edge(1) == 1
        % Contract environment from the right
        newright_env = 1;
        
        
        for j = N:-1:edge(2)+1
            newright_env = newcontract_right_environment(newright_env, M{j}, H{j}, M_{j});
        end
        newleft_env = 1;
        newmiddle_env = 1;
    elseif edge(1) == N-1
        % Contract environment from the left
        newleft_env = 1;
        for j = 1:edge(1)-1
            newleft_env = newcontract_left_environment(newleft_env, M{j}, H{j}, M_{j});
        end
        newright_env = 1;
        newmiddle_env = 1;
    else
        % Contract environment from the left
        newleft_env = 1;
        for j = 1:edge(1)-1
            newleft_env = newcontract_left_environment(newleft_env, M{j}, H{j}, M_{j});
        end
        
        % Contract environment from the right
        newright_env = 1;
        for j =  N:-1:edge(2)+1
            newright_env = newcontract_right_environment(newright_env, M{j}, H{j}, M_{j});
        end
        newmiddle_env =1;
    end
end
