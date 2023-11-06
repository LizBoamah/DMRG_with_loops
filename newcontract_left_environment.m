function newcontracted_env = newcontract_left_environment(env, MPS_tensor,mpo_tensor, M_)
% Contract the left environment, an MPS tensor, and an MPO tensor
temp = tensorprod(env, MPS_tensor, 1, 1, "NumDimensionsA", 3);
temp = tensorprod(temp, mpo_tensor,[1, 3], [1, 2], "NumDimensionsA",4 );
newcontracted_env = tensorprod(temp, M_, [1, 3], [1, 2], "NumDimensionsA",4);
