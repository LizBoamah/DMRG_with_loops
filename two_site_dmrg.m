function [dmrg,E_exact, energy_values,new_energy_values, E_DMRG_first,E_DMRG_second] = two_site_dmrg(N, mu, P,D, t, max_sweeps)
    % two_site_dmrg: Solves the ground state of a 1D Hubbard model using DMRG
    % 
    % Inputs:
    %   N          : Number of lattice sites
    %   mu         : Chemical potential
    %   D          : Maximum bond dimension
    %   U          : On-site interaction strength
    %   t          : Hopping parameter between adjacent sites
    %   max_sweeps : Maximum number of DMRG sweeps
    %
    % Outputs:
    %   dmrg           : Final Matrix Product State (MPS) approximation for the ground state
    %   E_exact        : Ground state energy from exact diagonalization
    %   energy_sweeps  : Vector of energies at each sweep
    %   coefficients   : Coefficients from the singular value decomposition (SVD)
    % 
    
   
    % Initialize tensor network state and edge set using example0 function.
    % This sets up the initial configuration for our many-body quantum system.
        [E_set,E_conn, mps] = example0(P,D);


    % Assuming you have a TensorNetwork object named 'TN'
        numTensors = size(mps.tensor, 1); % Get the number of tensors in the TensorNetwork
        M = cell(1,numTensors);
        % M_dagger = cell(1,numTensors);
        % for index = 1:numTensors
        %     M{index} = mps.tensor(index).value;
        % end
        for index = 1:numTensors
            % Tensor value extraction and initialization of its conjugate
            M{index} = mps.tensor(index).value;
            M_dagger{index} = conj(M{index});
            % Normalization of each MPS Tensor
            M{index} = M{index} / norm(M{index}(:));
         end
   
           % all_symmetric = is_MPS_symmetric(M);
         

     % Construct the MPO for the Hubbard Hamiltonian
            H = hubbard_mpo_site(t,mu, N);
            Hmatrix = mpo_to_hamiltonian(H);
            [~,E_exact] = exact_diagonalization(Hmatrix);

            

    for sweep = 1:1
        % Forward sweep
        % % Initialize the energy_sweeps vector  [Energy after forward sweep, Energy after backward sweep]
               energy_values = zeros(sweep,N-1); 
        for i = 1:size(E_set, 1)-1
            edge = E_set(i, :);
                [left_env,middle_env, right_env] = contract_environments(M, M_dagger, H, edge, numTensors);
                effective_H = build_effective_hamiltonian(left_env, H{edge(1)}, H{edge(2)}, right_env, middle_env);

               % Diagonalize the effective Hamiltonian
                [eig_vec, eig_val] = eig((effective_H + effective_H')/2);
                [energy_new1,ground_state_index] = min(diag(eig_val));
                psi = eig_vec(:, ground_state_index);


            % Store the energy after each tensor optimization during forward sweep
               energy_values(sweep, edge(1,2)) = energy_new1;

                % Reshape psi and perform SVD
                 psi_matrix = reshape(psi, [size(M{edge(1)}, 1) * size(M{edge(1)}, 2), size(M{edge(2)}, 2) * size(M{edge(2)}, 3)]);
                [U, S, V] = svd(psi_matrix, 'econ');
             

            % Truncate to maximum bond dimension D
                if size(U, 2) > D
                    U = U(:, 1:D);
                    S = S(1:D, 1:D);
                    V = V(:, 1:D);
                end
                % Update M and M_dagger using the ground state wavefunction
                % This typically involves reshaping psi and possibly performing an SVD to maintain the MPS form
                  if edge(1) == 1
                    % Update M{edge(1)} and M{edge(2)}
                    M{edge(1)} = reshape(U, [size(left_env,1), size(H{edge(1)}, 2), size(S, 1)]);
                    M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 2), size(right_env, 1)]);
                elseif edge(1) == N-1
                    M{edge(1)} = reshape(U , [size(left_env, 1), size(H{edge(1)}, 2), size(S, 1)]);
                    M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 3), size(right_env,2)]);
                  else
                       M{edge(1)} = reshape(U, [size(left_env, 1), size(H{edge(1)}, 2), size(S, 1)]);
                       M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 3), size(right_env, 3)]);
                  end
                    % Update the conjugate MPS tensor M_dagger
                    M_dagger{edge(1)} = conj(M{edge(1)});
                    M_dagger{edge(2)} = conj(M{edge(2)});
            
        end
        
        % Backward sweep
        for i = size(E_set,1)-1:-1:1
            edge = E_set(i, :);
            [left_env,middle_env, right_env] = contract_environments(M, M_dagger, H, edge, numTensors);
            effective_H = build_effective_hamiltonian(left_env, H{edge(1)}, H{edge(2)}, right_env, middle_env);
            
            % Diagonalize the effective Hamiltonian
            [eig_vec, eig_val] = eig((effective_H + effective_H')/2);
            energy_new2 = min(diag(eig_val));
            idx = find(diag(eig_val) == eig_val, 1);
            psi = eig_vec(:, idx);

            % Store the energy after each tensor optimization during backward sweep
            energy_values(sweep, edge(1,2)) = energy_new2;

            % Reshape psi and perform SVD
            psi_matrix = reshape(psi, [size(M{edge(1)}, 1) * size(M{edge(1)}, 2), size(M{edge(2)}, 2) * size(M{edge(2)}, 3)]);
            [U, S, V] = svd(psi_matrix, 'econ');
            
            % Truncate to maximum bond dimension D
            if size(U, 2) > D
                U = U(:, 1:D);
                S = S(1:D, 1:D);
                V = V(:, 1:D);
            end
    
            % Update M and M_dagger using the ground state wavefunction
             if edge(1) == 1
                    % Update M{edge(1)} and M{edge(2)}
                    M{edge(1)} = reshape(U*S, [size(left_env,1), size(H{edge(1)}, 2), size(S, 1)]);
                    M{edge(2)} = reshape(V', [size(S, 2), size(H{edge(2)}, 2), size(right_env, 1)]);
                elseif edge(1) == N-1
                    M{edge(1)} = reshape(U*S , [size(left_env, 1), size(H{edge(1)}, 2), size(S, 1)]);
                    M{edge(2)} = reshape(V', [size(S, 2), size(H{edge(2)}, 3), size(right_env,2)]);
                  else
                       M{edge(1)} = reshape(U*S, [size(left_env, 1), size(H{edge(1)}, 2), size(S, 1)]);
                       M{edge(2)} = reshape(V', [size(S, 2), size(H{edge(2)}, 3), size(right_env, 3)]);
             end
            
            % Update the conjugate MPS tensor M_dagger
            M_dagger{edge(1)} = conj(M{edge(1)});
            M_dagger{edge(2)} = conj(M{edge(2)});
        end
    end

 % Two-site update for the (N,1) edge after the first sweep

   % Loop over all connected edges from E_conn
for edge_idx = 1:size(E_conn,1)
    edge = E_conn(edge_idx,:);

    % Contract environments and build the effective Hamiltonian
    [left_env, middle_env, right_env] = contract_environments(M, M_dagger, H, edge, numTensors);
    effective_H = build_effective_hamiltonian(left_env, H{edge(1)}, H{edge(2)}, right_env, middle_env);

    % Diagonalize the effective Hamiltonian
    [eig_vec, eig_val] = eig((effective_H + effective_H')/2);
    energy_for_edge = min(diag(eig_val));
    idx = find(diag(eig_val) == energy_for_edge, 1);
    psi = eig_vec(:, idx);
    energy_values(sweep,edge(1,2))=energy_for_edge;


    % Reshape psi and perform SVD
    psi_matrix = reshape(psi, [size(M{edge(1)}, 1) * size(middle_env, 2), size(M{edge(2)}, 2) * size(middle_env, 3)]);
    [U, S, V] = svd(psi_matrix, 'econ');

    % Truncate to maximum bond dimension D
    if size(U, 2) > D
        U = U(:, 1:D);
        S = S(1:D, 1:D);
        V = V(:, 1:D);
    end
   % Update M and M_dagger for the current edge
     % Update M and M_dagger for the (N,1) edge
    M{edge(1)} = reshape(U*S, [ size(middle_env, 2), size(H{edge(1)}, 2),size(S, 1)]);
    M{edge(2)} = reshape(V',[ size(S,2),size(H{edge(2)}, 2) ,size(middle_env, 2)]);

    % Update the conjugate MPS tensor M_dagger
    M_dagger{edge(1)} = conj(M{edge(1)});
    M_dagger{edge(2)} = conj(M{edge(2)});
end

        % Continue with subsequent sweeps (from 2 to max_sweeps) with periodic MPS
        new_energy_values= zeros(max_sweeps,N);
    for sweep = 2:max_sweeps
        % forward and backward sweep code including the (N,1) edge
        % Forward sweep
        for i = 1:size(E_set, 1)
            edge = E_set(i, :);
            % contract the environments
            [updated_left_env, updated_middle_env, updated_right_env] = updated_contract_environments(M, M_dagger, H, edge, numTensors);
            % compute the  new effective hamiltonian
             updated_effective_H = build_updated_effective_hamiltonian(updated_left_env, H{edge(1)}, H{edge(2)}, updated_right_env, updated_middle_env);

                % Diagonalize the effective Hamiltonian
                 [eig_vec, eig_val] = eig((updated_effective_H+updated_effective_H')/2);
                 energy_new1 = sort(diag(eig_val),'ascend');
                 energy_new1 = energy_new1(1);
                 idx = find(diag(eig_val) == eig_val, 1);
                  psi = eig_vec(:, idx);
                    new_energy_values(sweep,:) = energy_new1;

                % Reshape psi and perform SVD

 % Calculate the new energy after completing a left-to-right sweep
       % Reshape psi and perform SVD for edge(1)=1
         if edge(1)== 1
             % Reshape psi and perform SVD for other values of edge
             psi_matrix = reshape(psi, [size(M{edge(1)}, 1) * size(M{edge(1)}, 2), size(M{edge(2)}, 2) * size(M{edge(2)}, 3)]);
        elseif edge(1) == N-1
            % Reshape psi and perform SVD for edge(1)=N-1
              psi_matrix = reshape(psi, [size(M{edge(1)}, 1) * size(M{edge(2)}, 2), size(M{edge(2)}, 2) * size(M{edge(2)}, 3)]);
      elseif edge(1) ==  N
           % Reshape psi and perform SVD for edge(1)=N
           psi_matrix = reshape(psi, [size(M{edge(1)}, 1) * size(M{edge(2)}, 2), size(M{edge(2)}, 2) * size(M{edge(2)}, 3)]);
        else
            % Reshape psi and perform SVD for other values of edge
               psi_matrix = reshape(psi, [size(M{edge(1)}, 1) * size(M{edge(1)}, 2), size(M{edge(2)}, 2) * size(M{edge(2)}, 3)]);
        end

                [U, S, V] = svd(psi_matrix, 'econ');
                % coefficients = diag(S);

            % Truncate to maximum bond dimension D
                if size(U, 2) > D
                    U = U(:, 1:D);
                    S = S(1:D, 1:D);
                    V = V(:, 1:D);
                end
                % Update M and M_dagger using the ground state wavefunction
                % This typically involves reshaping psi and possibly performing an SVD to maintain the MPS form
                  if edge(1) == 1
                    % Update M{edge} and M{edge+1}
                     M{edge(1)} = reshape(U, [size(updated_right_env,1), size(H{edge(1)}, 2), size(S, 1)]);
                     M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 3), size(updated_right_env, 5)]);
                elseif edge(1) == N-1
                    % Update M{edge} and M{edge+1}
                    M{edge(1)} = reshape(U , [size(updated_left_env, 1), size(H{edge(1)}, 2), size(S, 1)]);
                    M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 3), size(updated_left_env,2)]);
                  elseif edge(1) ==  N
                     M{edge(1)} = reshape(U, [size(updated_middle_env, 1), size(H{edge(1)}, 2), size(S, 1)]);
                     M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 3), size(updated_middle_env, 5)]);
                  else
                       M{edge(1)} = reshape(U, [size(updated_left_env, 2), size(H{edge(1)}, 2),size(S, 1)]);
                       M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 3),size(updated_right_env, 1)]);
                 end

                % Update the conjugate MPS tensor M_dagger
                M_dagger{edge(1)} = conj(M{edge(1)});
                M_dagger{edge(2)} = conj(M{edge(2)});

        end

        % Backward sweep
        for i = size(E_set,1):-1:1
              edge = E_set(i, :);
             [updated_left_env, updated_middle_env, updated_right_env] = updated_contract_environments(M, M_dagger, H, edge, numTensors);
             % effective_H = build_effective_hamiltonian(updated_left_env, H{edge(1)}, H{edge(2)}, updated_right_env, updated_middle_env);
             updated_effective_H = build_updated_effective_hamiltonian(updated_left_env, H{edge(1)}, H{edge(2)}, updated_right_env, updated_middle_env);

                % Diagonalize the effective Hamiltonian
                 [eig_vec, eig_val] = eig((updated_effective_H+updated_effective_H')/2);
                % [eig_vec, eig_val] = eigs(effective_H,1,'smallestreal');
                 energy_new1 = sort(diag(eig_val),'ascend');
                 energy_new1 = energy_new1(1);
                 idx = find(diag(eig_val) == eig_val, 1);
                  psi = eig_vec(:, idx);
                     new_energy_values(sweep,:) = energy_new1;

            % Reshape psi and perform SVD
       % Reshape psi and perform SVD for edge(1)=1
         if edge(1)== 1
             % Reshape psi and perform SVD for other values of edge
              psi_matrix = reshape(psi, [size(M{edge(1)}, 1) * size(M{edge(1)}, 2), size(M{edge(2)}, 2) * size(M{edge(2)}, 3)]);
            elseif edge(1) == N-1
                % Reshape psi and perform SVD for edge(1)=N-1
                  psi_matrix = reshape(psi, [size(M{edge(1)}, 2) * size(M{edge(1)}, 1), size(M{edge(2)}, 2) * size(M{edge(2)}, 3)]);
          elseif edge(1) ==  N
               % Reshape psi and perform SVD for edge(1)=N
               psi_matrix = reshape(psi, [size(M{edge(1)}, 1) * size(M{edge(2)}, 2), size(M{edge(2)}, 2) * size(M{edge(2)}, 3)]);
            else
                % Reshape psi and perform SVD for other values of edge
                   psi_matrix = reshape(psi, [size(M{edge(1)}, 1) * size(M{edge(1)}, 2), size(M{edge(2)}, 2) * size(M{edge(2)}, 3)]);
        end
            % psi_matrix = reshape(psi, [size(M{edge(1)}, 1) * size(M{edge(1)}, 2), size(M{edge(2)}, 2) * size(M{edge(2)}, 3)]);
            [U, S, V] = svd(psi_matrix, 'econ');

            % Truncate to maximum bond dimension D
            if size(U, 2) > D
                U = U(:, 1:D);
                S = S(1:D, 1:D);
                V = V(:, 1:D);
            end
              % Update M and M_dagger using the ground state wavefunction
                % This typically involves reshaping psi and possibly performing an SVD to maintain the MPS form
                  if edge(1) == 1
                    % Update M{s} and M{s+1}
                     M{edge(1)} = reshape(U*S, [size(updated_right_env,1), size(H{edge(1)}, 2), size(S, 1)]);
                     M{edge(2)} = reshape(V', [size(S, 2), size(H{edge(2)}, 3), size(updated_right_env, 5)]);
                elseif edge(1) == N-1
                    % Update M{s} and M{s+1}
                    M{edge(1)} = reshape(U*S , [size(updated_left_env, 2), size(H{edge(1)}, 2), size(S, 1)]);
                    M{edge(2)} = reshape(V', [size(S, 2), size(H{edge(2)}, 3), size(updated_left_env,1)]);
                  elseif edge(1) ==  N
                     M{edge(1)} = reshape(U*S, [size(updated_middle_env, 1), size(H{edge(1)}, 2), size(S, 1)]);
                     M{edge(2)} = reshape(V', [size(S, 2), size(H{edge(2)}, 3), size(updated_middle_env, 5)]);
                  else
                       % M{edge(1)} = reshape(U, [size(M{edge(1)}, 3), size(H{edge(1)}, 2), size(S, 1)]);
                       % M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 3), size(M{edge(2)}, 1)]);
                       M{edge(1)} = reshape(U*S, [size(updated_left_env, 2), size(H{edge(1)}, 2),size(S, 1)]);
                       M{edge(2)} = reshape(V', [size(S, 2), size(H{edge(2)}, 3),size(updated_right_env, 1)]);
                 end


           

            % Update the conjugate MPS tensor M_dagger
            M_dagger{edge(1)} = conj(M{edge(1)});
            M_dagger{edge(2)} = conj(M{edge(2)});
        end


   end


     dmrg = M;  % Final MPS state
     E_DMRG_first = energy_values(end);  % The last energy value after the 1st full sweep(Open boundary chain)
       E_DMRG_second = new_energy_values(end);  % The last energy value after the max full sweeps(periodic chain)

    end

    

















































































% [left_env,middle_env, right_env] = new_contract_environments(M, M_dagger, H, edge, numTensors);
    %             effective_H = build_effective_hamiltonian(left_env, H{edge(1)}, H{edge(2)}, right_env, middle_env);
    % 
    %            % Diagonalize the effective Hamiltonian
    %             [eig_vec, eig_val] = eig((effective_H + effective_H')/2);
    %             energy_new1 = min(diag(eig_val));
    %             idx = find(diag(eig_val) == eig_val, 1);
    %             psi = eig_vec(:, idx);
    % 
    %         % Store the energy after each tensor optimization during forward sweep
    %            energy_values(sweep, edge(1)) = energy_new1;
    % 
    %             % Reshape psi and perform SVD
    %              psi_matrix = reshape(psi, [size(M{edge(1)}, 1) * size(M{edge(1)}, 2), size(M{edge(2)}, 2) * size(M{edge(2)}, 3)]);
    %             [U, S, V] = svd(psi_matrix, 'econ');
    %             coefficients = diag(S);
    % 
    %         % Truncate to maximum bond dimension D
    %             if size(U, 2) > D
    %                 U = U(:, 1:D);
    %                 S = S(1:D, 1:D);
    %                 V = V(:, 1:D);
    %             end
    %             % Update M and M_dagger using the ground state wavefunction
    %             % This typically involves reshaping psi and possibly performing an SVD to maintain the MPS form
    %               if edge(1) == 1
    %                 % Update M{s} and M{s+1}
    %                 M{edge(1)} = reshape(U, [size(left_env,1), size(H{edge(1)}, 2), size(S, 1)]);
    %                 M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 2), size(right_env, 1)]);
    %             elseif edge(1) == N-1
    %                 % Update M{s} and M{s+1}
    %                 M{edge(1)} = reshape(U , [size(left_env, 1), size(H{edge(1)}, 2), size(S, 1)]);
    %                 M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 3), size(right_env,2)]);
    % 
    %               else
    %                    M{edge(1)} = reshape(U, [size(left_env, 1), size(H{edge(1)}, 2), size(S, 1)]);
    %                    M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 3), size(right_env, 3)]);
    %              end
    % 
    %             % Update the conjugate MPS tensor M_dagger
    %             M_dagger{edge(1)} = conj(M{edge(1)});
    %             M_dagger{edge(2)} = conj(M{edge(2)});
    % 
    %     end
    % 
    %     % Backward sweep
    %     for i = size(E_set,1):-1:1
    %         edge = E_set(i, :);
    %         [left_env,middle_env, right_env] = contract_environments(M, M_dagger, H, edge, numTensors);
    %         effective_H = build_effective_hamiltonian(left_env, H{edge(1)}, H{edge(2)}, right_env, middle_env);
    % 
    %         % Diagonalize the effective Hamiltonian
    %         [eig_vec, eig_val] = eig((effective_H + effective_H')/2);
    %         energy_new2 = min(diag(eig_val));
    %         idx = find(diag(eig_val) == eig_val, 1);
    %         psi = eig_vec(:, idx);
    % 
    %         % Store the energy after each tensor optimization during backward sweep
    %         energy_values(sweep, edge(1)) = energy_new2;
    % 
    %         % Reshape psi and perform SVD
    %         psi_matrix = reshape(psi, [size(M{edge(1)}, 1) * size(M{edge(1)}, 2), size(M{edge(2)}, 2) * size(M{edge(2)}, 3)]);
    %         [U, S, V] = svd(psi_matrix, 'econ');
    % 
    %         % Truncate to maximum bond dimension D
    %         if size(U, 2) > D
    %             U = U(:, 1:D);
    %             S = S(1:D, 1:D);
    %             V = V(:, 1:D);
    %         end
    % 
    %         % Update M and M_dagger using the ground state wavefunction
    %          if edge(1) == 1
    %                 % Update M{s} and M{s+1}
    %                 M{edge(1)} = reshape(U, [size(left_env,1), size(H{edge(1)}, 2), size(S, 1)]);
    %                 M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 2), size(right_env, 1)]);
    %             elseif edge(1) == N-1
    %                 % Update M{s} and M{s+1}
    %                 M{edge(1)} = reshape(U , [size(left_env, 1), size(H{edge(1)}, 2), size(S, 1)]);
    %                 M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 3), size(right_env,2)]);
    % 
    %               else
    %                    M{edge(1)} = reshape(U, [size(left_env, 1), size(H{edge(1)}, 2), size(S, 1)]);
    %                    M{edge(2)} = reshape(S*V', [size(S, 2), size(H{edge(2)}, 3), size(right_env, 3)]);
    %          end
    % 
    %         % Update the conjugate MPS tensor M_dagger
    %         M_dagger{edge(1)} = conj(M{edge(1)});
    %         M_dagger{edge(2)} = conj(M{edge(2)});
    %     end
    % end












