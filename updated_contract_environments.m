function [updated_left_env, updated_middle_env, updated_right_env] = updated_contract_environments(M, M_dagger, H, edge, numTensors)
    % This function calculates the left, middle, and right environments
    % given the current state M, its conjugate M_dagger, the Hamiltonian H,
    % the current edge under consideration, and the total number of tensors.
    
    % Initialize middle_env to an empty array
    updated_middle_env = 1;
    
    % Special case where edge is between the last and first sites
    if edge(1) == numTensors && edge(2) == 1
        % Initialize left and right environments to identity
        updated_left_env = 1;
        updated_right_env = 1;
        
        % Compute middle environment
        updated_middle_env = tensorprod(M{2}, H{2}, 2,2);
        updated_middle_env = tensorprod(updated_middle_env, M_dagger{2}, 4,2);

        % Contract remaining tensors to build middle environment
        for i = 3:(numTensors-1)
            temp_env = tensorprod(M{i}, H{i}, 2,2);
            temp_env = tensorprod(temp_env, M_dagger{i}, 4,2);
            updated_middle_env = tensorprod(updated_middle_env, temp_env, [2, 4, 6],[1, 3, 5]);
             updated_middle_env = permute(updated_middle_env, [1 4 2 5 3 6]);
        end

    else  % General case
        % Initialize left environment if edge begins at the first site
        if edge(1) == 1
            updated_left_env = 1;
           
        else
            % Contract tensors to build left environment
            updated_left_env = tensorprod(M{1}, H{1}, 2,2);
            updated_left_env = tensorprod(updated_left_env, M_dagger{1}, 4,2);
          
        end

        % Continue contracting tensors for left environment up to the edge
        for i = 2:edge(1)-1
            temp_env = tensorprod(M{i}, H{i}, 2, 2);
            temp_env = tensorprod(temp_env, M_dagger{i}, 4, 2);
            updated_left_env = tensorprod(updated_left_env, temp_env, [2, 4 ,6], [1, 3, 5]);
             % updated_left_env = permute(updated_left_env,[1,4,2,5,3,6]);
             updated_left_env = permute(updated_left_env,[1,4,2,5,3,6]);
        end

        % Initialize right environment if edge ends at the last site
        if edge(2) == numTensors
            updated_right_env = 1;
            
        else
            % Contract tensors to build right environment
             updated_right_env = tensorprod(M{numTensors}, H{numTensors}, 2, 2);
             updated_right_env = tensorprod(updated_right_env, M_dagger{numTensors},4,2);
              updated_right_env = permute(updated_right_env,[1,2,3,6,4,5]);
             
        end

        % Continue contracting tensors for right environment up to the edge
        for i = (numTensors-1):-1:(edge(2)+1)
            temp_env = tensorprod(M{i}, H{i}, 2, 2);
            temp_env = tensorprod(temp_env, M_dagger{i}, 4, 2);
            updated_right_env = tensorprod(temp_env,updated_right_env, [2,4, 6], [1, 3, 5]);
            updated_right_env = permute(updated_right_env,[1,4,2,5,3,6]);
        end
    end
end

































































% function [left_env, middle_env, right_env] = contract_environments(M, M_dagger, H, edge, numTensors)
%     % Special case for edge = [N 1]
%     middle_env = [];
%     if edge(1) == numTensors && edge(2) == 1
%         % Initialize left_env, right_env and middle_env
%         left_env = 1;
%         right_env = 1;
% 
%         middle_env = tensorprod(M{2}, H{2}, 2,2);
%         middle_env = tensorprod(middle_env, M_dagger{2}, 4,2);
% 
%         % Contract the rest of the tensors from 3 to N-1
%         for i = 3:(numTensors-1)
%             temp_env = tensorprod(M{i}, H{i}, 2,2);
%             temp_env = tensorprod(temp_env, M_dagger{i}, 4,2);
%             middle_env = tensorprod(middle_env, temp_env,[2, 4, 5], [2,4,5]);
%             middle_env = permute(middle_env,[1 3 2 5 4 6]);
%         end
%         
%     else
%         if edge(1) == 1
%               left_env = 1;
%         else
%               % Contract the left environment
%                 left_env = tensorprod(M{1}, H{1}, 2,1);
%                 left_env = tensorprod(left_env, M_dagger{1},3 ,2);
%         end
%         % Contract the rest of the tensors up to the edge
%         for i = 2:edge(1)-1
%             temp_env = tensorprod(M{i}, H{i}, 2, 2);
%             temp_env = tensorprod(temp_env, M_dagger{i}, 4, 2);
%             left_env = tensorprod(left_env, temp_env, [1, 2 ,3], [1, 2, 3]);
%         end
%         if edge(2) == numTensors
%             right_env = 1;
%         else
%             % Contract the right environment
%             right_env = tensorprod(M{numTensors}, H{numTensors}, 2, 2);
%             right_env = tensorprod(right_env, M_dagger{numTensors}, 4,2);
%         end
%             % Contract the rest of the tensors from the end to the edge
%     for i = (numTensors-1):-1:(edge(2)+1)
%         temp_env = tensorprod(M{i}, H{i}, 2,2);
%         temp_env = tensorprod(temp_env, M_dagger{i}, 4,2);
%         right_env = tensorprod(right_env, temp_env,[1, 2, 3], [1, 2, 3]);
%     end
%     end
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 


























  




























