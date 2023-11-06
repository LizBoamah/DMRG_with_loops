function all_symmetric = is_MPS_symmetric(MPS)
    num_tensors = length(MPS);
    all_symmetric = true;  % Start with the assumption that all are symmetric

    for i = 1:num_tensors
        tensor = MPS{i};

        % Conjugate transpose over the physical index
        tensor_conj_transpose = conj(tensor);

        if ~isequal(tensor, tensor_conj_transpose)
            all_symmetric = false;
            return;
        end
    end
end
% 
% % Example usage:
% N = 10; % Number of sites
% d = 2;  % Physical dimension (e.g., spin-1/2)
% MPS = cell(1, N);
% 
% % Just an example to fill the MPS; this won't be symmetric
% for i = 1:N
%     MPS{i} = rand([2, d, 2]);
% end

% result = is_MPS_symmetric(MPS);
