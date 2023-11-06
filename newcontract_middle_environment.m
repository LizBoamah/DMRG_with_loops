function newmiddle_env = newcontract_middle_environment(M,H,M_dagger,numTensors)
% Initialize middle_env to an empty array
        
        % Compute middle environment
        newmiddle_env = tensorprod(M{2}, H{2}, 2,2);
        newmiddle_env = tensorprod(newmiddle_env, M_dagger{2}, 4,2);

        % Contract remaining tensors to build middle environment
        for i = 3:(numTensors-1)
            temp_env = tensorprod(M{i}, H{i}, 2,2);
            temp_env = tensorprod(temp_env, M_dagger{i}, 4,2);
            newmiddle_env = tensorprod(newmiddle_env, temp_env, [2, 4, 6], [1, 3, 5]);
            newmiddle_env = permute(newmiddle_env, [1 4 3 5 2 6]);
        end
    end