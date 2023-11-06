function contracted_env = contract_right_environment(env, MPS_tensor,mpo_tensor, M_)
% Contract the right environment, an MPS tensor, and an MPO tensor
% Contract the left environment, an MPS tensor, and an MPO tensor
temp = tensorprod(env,MPS_tensor, 1, 3, "NumDimensionsA", 3);
temp = tensorprod(temp, mpo_tensor,[1, 4], [4, 2], "NumDimensionsA",4 );
contracted_env = tensorprod(temp, M_, [1, 4], [3, 2], "NumDimensionsA",4);

end