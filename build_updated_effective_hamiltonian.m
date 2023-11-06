function updated_effective_H = build_updated_effective_hamiltonian(left_env, H_edge1, H_edge2, right_env, middle_env)
    % This function constructs the effective Hamiltonian for a given edge
    % based on the left, right, and middle environments.

    % Case 1: Special case where both left and right environments are scalars and equal to 1
    % and middle environment exists
    if isscalar(left_env) && all(left_env == 1) && isscalar(right_env) && all(right_env == 1) && exist('middle_env', 'var') && ~isempty(middle_env)
        
        % Contract middle environment with the Hamiltonian tensor at edge 1
        temp1 = tensorprod(H_edge2, middle_env, 4, 3);

        % Further contract with the Hamiltonian tensor at edge 2
        updated_effective_H = tensorprod(temp1, H_edge1, 6, 1);

         % updated_effective_H = permute(updated_effective_H, [1,3,5,7,9,2,4,8,6,10]);
           % updated_effective_H = permute(effective_H, [2,4,6,8,10,1,3,5,7,9,]);


        % Calculate dimensions for the subsequent reshape operation
        siz1 = (size(updated_effective_H, 1)*size(updated_effective_H, 2)*size(updated_effective_H, 3)*size(updated_effective_H, 4)*size(updated_effective_H, 5));
        siz2 = (size(updated_effective_H, 6)*size(updated_effective_H, 7)*size(updated_effective_H, 8)*size(updated_effective_H, 9)*size(updated_effective_H, 10));

    % Case 2: Special case where only the left environment is 1
    elseif left_env == 1

        % Contract left environment with Hamiltonian tensor at edge 1
        temp1 = tensorprod(left_env, H_edge1, 2, 1);

        % Contract right environment with Hamiltonian tensor at edge 2
         temp2 = tensorprod(right_env, H_edge2, 3, 4);
       
            % Contract the two resulting tensors
            updated_effective_H = tensorprod(temp1, temp2, 4, 6);
            % Permute the tensor for proper reshaping
                 updated_effective_H = permute(updated_effective_H, [1,3,7,5,9,2,4,8,6,10]);
               % updated_effective_H = permute(updated_effective_H, [1,3,5,7,9,2,4,6,8,10]);
           

        % Calculate dimensions for the subsequent reshape operation
        siz1 = (size(updated_effective_H, 1)*size(updated_effective_H, 2)*size(updated_effective_H, 3)*size(updated_effective_H, 4)*size(updated_effective_H, 5));
        siz2 = (size(updated_effective_H, 6)*size(updated_effective_H, 7)*size(updated_effective_H, 8)*size(updated_effective_H, 9)*size(updated_effective_H, 10));

    % Case 3: Special case where only the right environment is 1
    elseif right_env == 1

        % Contract left environment with Hamiltonian tensor at edge 1
        temp1 = tensorprod(left_env, H_edge1, 4, 1);

        % Contract right environment with Hamiltonian tensor at edge 2
        temp2 = tensorprod(right_env, H_edge2, 2, 4);

        % Contract the two resulting tensors
        updated_effective_H = tensorprod(temp1, temp2, 8, 2);
        % Permute the tensor for proper reshaping
            % updated_effective_H = permute(updated_effective_H, [1,2,3,6,9,4,5,7,8,10]);
            updated_effective_H = permute(updated_effective_H, [1,3,5,7,9,2,4,6,8,10]);

        % Calculate dimensions for the subsequent reshape operation
        siz1 = (size(updated_effective_H, 1)*size(updated_effective_H, 2)*size(updated_effective_H, 3)*size(updated_effective_H, 4)*size(updated_effective_H, 5));
        siz2 = (size(updated_effective_H, 6)*size(updated_effective_H, 7)*size(updated_effective_H, 8)*size(updated_effective_H, 9)*size(updated_effective_H, 10));

    % Case 4: General case (neither left nor right environments are 1)
    else

        % Contract left environment with Hamiltonian tensor at edge 1
        temp1 = tensorprod(left_env, H_edge1, 4, 1);

        % Contract right environment with Hamiltonian tensor at edge 2
        temp2 = tensorprod(right_env, H_edge2, 3, 4);

        % Contract the two resulting tensors
         % updated_effective_H  = tensorprod(temp1, temp2, [2,5,8], [1,3,5]);
         updated_effective_H  = tensorprod(temp1, temp2, [1,4,8], [2,5,6]);
         % updated_effective_H = tensorprod(temp1, temp2, 8, 6);


        % Permute the tensor for proper reshaping
           % updated_effective_H = permute(updated_effective_H, [1,3,5,7,9,2,4,6,8,10]);
            updated_effective_H = permute(updated_effective_H, [3,5,7,9,8,1,2,4,6,10]);
          
          % Calculate dimensions for the subsequent reshape operation

        siz1 = (size(updated_effective_H, 1)*size(updated_effective_H, 2)*size(updated_effective_H, 3)*size(updated_effective_H, 4)*size(updated_effective_H, 5));
        siz2 = (size(updated_effective_H, 6)*size(updated_effective_H, 7)*size(updated_effective_H, 8)*size(updated_effective_H, 9)*size(updated_effective_H, 10));
         % siz1 = (size(updated_effective_H, 1)*size(updated_effective_H, 2)*size(updated_effective_H, 3));
         % siz2 = (size(updated_effective_H, 4)*size(updated_effective_H, 5)*size(updated_effective_H, 6));


       
    end

    % Reshape the effective Hamiltonian
    updated_effective_H = reshape(updated_effective_H, siz1, siz2);
end
