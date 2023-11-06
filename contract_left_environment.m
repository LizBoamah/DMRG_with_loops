function contracted_env = contract_left_environment(env, MPS_tensor,mpo_tensor, M_)
% Contract the left environment, an MPS tensor, and an MPO tensor
temp = tensorprod(env, MPS_tensor, 1, 1, "NumDimensionsA", 3);
temp = tensorprod(temp, mpo_tensor,[1, 3], [1, 2], "NumDimensionsA",4 );
contracted_env = tensorprod(temp, M_, [1, 3], [1, 2], "NumDimensionsA",4);







% function contracted_env = contract_left_environment(env, MPS_tensor, mpo_tensor, M_)
%     % Check if env is a scalar
%     if isscalar(env)
%         temp = tensorprod(MPS_tensor, mpo_tensor,[1, 3], [1, 2], "NumDimensionsA",4 );
%         contracted_env = tensorprod(temp, M_, [1, 3], [1, 2], "NumDimensionsA",4);
%     else
%         % Contract the left environment, an MPS tensor, and an MPO tensor
%         temp = tensorprod(env, MPS_tensor, 1, 1, "NumDimensionsA", 3);
%         temp = tensorprod(temp, mpo_tensor,[1, 3], [1, 2], "NumDimensionsA",4 );
%         contracted_env = tensorprod(temp, M_, [1, 3], [1, 2], "NumDimensionsA",4);
%     end
% end


























%  function contracted_env = contract_left_environment(env, MPS_tensor,mpo_tensor, M_)
%             % Contract tensors in the left direction
%         M_perm = permute(M_, [1, 3, 2]);
%          M_reshaped = reshape(M_perm, [size(MPS_tensor, 1) * size(MPS_tensor, 3), size(MPS_tensor, 2)]);
%         
%         MPO_perm = permute(mpo_tensor, [3, 1, 2, 4]);
%          MPO_reshaped = reshape(MPO_perm, [size(mpo_tensor, 3), size(mpo_tensor, 1) * size(mpo_tensor, 2) * size(mpo_tensor, 4)]);
%         
%         temp1 = M_reshaped * MPO_reshaped;
% %         temp1 = M_perm * MPO_perm;
%         temp1 = temp1*env;
% 
%         temp2 = reshape(temp1, [size(MPS_tensor,1), size(MPS_tensor,2), size(mpo_tensor,3), size(mpo_tensor,4)]);
%         E1 = permute(temp2, [1, 3, 2, 4]);
%         temp3 = reshape(E1, [size(E1,1)*size(E1,2), size(E1,3)*size(E1,4)]) * reshape(M_, [size(M_,1), size(M_,2)*size(M_,3)]);
%         contracted_env = reshape(temp3, [size(E1,1), size(E1,2), size(M_,2), size(M_,3)]);
% end