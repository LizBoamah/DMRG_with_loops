function trace_val = traceMPS(A)
    % A is a cell array representing the MPS, A{1}, A{2}, ..., A{N}
    N = length(A); % Number of sites
 % Start with the tensor at the first site
    contracted_tensor = A{1};

    % Contract tensors site by site
    for site = 2:N
        % Define the dimensions to be contracted for the tensors
        dimA = site + 1; % The third index for the first tensor, fourth for the second, and so on
        dimB = 1; % The first index for the site tensor

        % Contract current tensor with the next tensor in the MPS
        contracted_tensor = tensorprod(contracted_tensor, A{site}, dimA, dimB);
    end       

    % Contract the last indices to close the loop and get the trace
    trace_val = tensorprod(contracted_tensor, contracted_tensor, [1, N+2], [N+2, 1]);
end
  
