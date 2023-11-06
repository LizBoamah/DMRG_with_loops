function H_mpo = mpo_to_hamiltonian(MPO)

    N = numel(MPO); % assuming MPO is a cell array
    if N < 2
        error('MPO must contain at least two elements.');
    end
    nextDim = 6; % Start with 6 as mentioned 6-D when i=3
    temp = tensorprod(MPO{1}, MPO{2}, 4, 1);
    for i = 3:N
        temp1 = tensorprod(temp, MPO{i}, nextDim, 1);
        temp = temp1; % Update temp to hold the result of the current iteration
        nextDim = nextDim + 2; % Increase dimensionality by 2 for the next iteration
    end   
    % Get the size of the original tensor
       sz = size(temp1);

    % Number of dimensions of the tensor
      numDims = length(sz);

    % Generate permutation order dynamically
       permOrder = [1:2:numDims, 2:2:numDims];
    

    % Permute the tensor
    B = permute(temp1, permOrder);

    
    % Calculate the products of the sizes of odd and even dimensions
    s1 = prod(sz(1:2:numDims));
    s2 = prod(sz(2:2:numDims));

    % Reshape the permuted tensor
    H_mpo = reshape(B, [s1, s2]);

end
