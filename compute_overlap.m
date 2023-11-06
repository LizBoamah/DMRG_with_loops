function [overlap] = compute_overlap(M)
    % Create a conjugate version of the tensors
    M_conj = cellfun(@(tensor) conj(tensor), M, 'UniformOutput', false);

    % Contract the tensors over their physical indices
    contracted_tensors = cell(1, length(M));
    for m = 1:length(M)
        % Contract tensor m from M with tensor m from M_conj over physical index
        % Assuming physical index is at position 3
        contracted_tensors{m} = tensorprod(M{m}, M_conj{m}, 2); 
    end

    % Compute the overlap
    overlap = 1;
    for m = 1:length(contracted_tensors)
        overlap = overlap * sum(contracted_tensors{m}(:));
    end
 % overlap = 1;
 %    for i = 1:length(M)
 %        overlap = overlap * trace(M_dagger{i} * M{i});
 %    end
 % 

% function result = contract(T1, T2, index)
%     % A simplistic tensor contraction function over specified indices
%     % This will need to be made more robust for tensors with more than 3 indices, for example
%     result = sum(bsxfun(@times, T1, T2), index);
% end
