function M_normalized = normalize_mps(M, overlap_val)
    % Compute the normalization factor
    norm_factor = sqrt(overlap_val);
    
    % Normalize each tensor in M
    M_normalized = cellfun(@(tensor) tensor / norm_factor, M, 'UniformOutput', false);
end
