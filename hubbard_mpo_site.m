 %  function MPO = hubbard_mpo_site(U, t,mu, N)
% % % Constructs the MPO tensor for the Hubbard Hamiltonian in the second
% % % quantized form for a one-dimensional chain of N sites.
% % % Input:
% % %   U: onsite interaction strength
% %   t: hopping parameter
% %   N: number of sites
% %   d: bond dimension
% % Output:
% %   MPO: cell array containing MPO tensors
% 
% % Create MPO for Hubbard Hamiltonian
% MPO = cell(1, N);
% % d = 2;
% % D= 10;
% 
% % Identity and fermion creation/annihilation matrices
% I = eye(2);
% c = [0, 1; 0, 0]; % annihilation operator
% c_dag = c';       % creation operator
% n = c_dag * c;    % number operator
% ph = diag([1, -1]);
% 
% % Construct MPO tensors
% for i = 1:N
%     if i == 1
%         W = zeros(1, 2, 2, 4);
%         W(1, :, :, 1) = I;
%         W(1, :, :, 2) = -t*c;
%         W(1, :, :, 3) = t *c_dag;
%         W(1, :, :, 4) = mu *n;
%     elseif i == N
%         W = zeros(4,2,2,1);
%         W(1, :, :, 1) = mu*n;
%         W(2, :, :, 1) = c;
%         W(3, :, :, 1) = c_dag;
%         W(4, :, :, 1) = I;
%     else
%         W = zeros(4, 2, 2, 4);
%         W(1, :, :, 1) = I;
%         W(2, :, :, 2) = -t*c;
%         W(3, :, :, 3) = t *c_dag;
%         W(4, :, :, 4) =mu*n;
%         W(2, :, :, 4) = c;
%         W(3, :, :, 4) = c_dag;
%         W(4, :, :, 1) = I;
%     end
%     MPO{i} = W;
% end



 function MPO = hubbard_mpo_site(t,mu, N)
% Constructs the MPO tensor for the Hubbard Hamiltonian in the second
% quantized form for a one-dimensional chain of N sites.
% Input:
%   U: onsite interaction strength
%   t: hopping parameter
%   N: number of sites
%   d: bond dimension
% Output:
%   MPO: cell array containing MPO tensors
% 
% Create MPO for Hubbard Hamiltonian
MPO = cell(1, N);
% d = 2;
% D= 10;

% Identity and fermion creation/annihilation matrices
% Define local operators in the spinless fermion basis: |0>, |1>
I = eye(2);
c = [0, 1; 0, 0]; % annihilation operator
c_dag = c';       % creation operator
n = c_dag * c;    % number operator
ph = diag([1, -1]);

% Construct MPO tensors
for i = 1:N
    if i == 1
        W = zeros(1,2, 2, 4);
        W(:,:, :, 1) = I;
        W(:, :, :, 2) = t*c_dag*ph;
        W(:,:, :, 3) = -t*c*ph;
        W(:, :, :,4) = mu*n;
    elseif i == N
        W = zeros(4,2,2,1);
        W(1, :, :,:) = mu*n;
        W(2, :, :,:) =  c ;
        W(3, :, :,:) =  c_dag ;
        W(4, :, :,:) = I;
    else
        W = zeros(4, 2, 2, 4);
        W(1, :, :, 1) = I;
        W(1, :, :, 2) = t *c_dag*ph;
        W(1, :, :, 3) = -t*c*ph;
        W(1, :, :, 4) = mu*n;
        W(2, :, :, 4) = c;  
        W(3, :, :, 4) = c_dag;
        W(4, :, :, 4) = I;
    end
    MPO{i} = W;
 

end













