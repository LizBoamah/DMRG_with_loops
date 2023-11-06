function [E,E_conn,TN] = example0(P,D)
    % % Define the parameters of the system
    % P = 2;  % Physical leg 
    % D = 10; % Maximum bond dimension

    % Define the edge connections for long-range and nearest-neighbor interactions
     E_conn = [6 1]; % Long range interactions (periodic edge)
    E_unconn = [1 2; 2 3; 3 4; 4 5; 5 6;]; % Nearest neighbour interactions

    %  E_conn = [8 1]; % Long range interactions (periodic edge)
    % E_unconn = [1 2; 2 3; 3 4; 4 5; 5 6;6 7;7 8;]; % Nearest neighbour interactions

    % Combine the long-range and nearest-neighbor interactions
    E =[ E_unconn;E_conn];

    % Compute the maximum index in the edge set to determine the number of tensors
    M = max(max(E));

    % Compute the coordination number for each tensor in the network
    for m=1:M
        e(m) = size(find(E==m),1);
    end   

    % Construct auxiliary matrices to store the bond dimensions for each tensor
    for i = 1:M
        [r,~] = find(E == i); % Find all edges connected to tensor i
        bd{i} = [E(r,:); i, i]; % Store the edge and physical dimensions
    end

    % Initialize the node tensors A(m), m=1,...,M with random values
    A = [];
    for m=1:M
        if m == 1
           sizevec = [1,P,D*ones(1,ceil(e(m)/2))];  % Physical index in the middle and the first site
        elseif m == M
              sizevec = [D*ones(1,floor(e(m)/2)),1,P];  % Physical index in the middle and the last site
           
        else
            sizevec = [D*ones(1,floor(e(m)/2)), P, D*ones(1,ceil(e(m)/2))];  % Physical index in the middle and the rest of the sites
        end
        A = [A; Tensor(rand(sizevec), bd{m})]; % Create a tensor with random entries and attach it to the TensorNetwork
    end

    % Create the TensorNetwork object with the initialized node tensors and edge set
    TN = TensorNetwork(A, E);

    % Retrieve the tensor for the last site
    last_tensor_data = TN(end).tensor(end).value;  % Assuming the Tensor object has a 'data' property
    
    % Permute its dimensions
    last_tensor_permuted = permute(last_tensor_data, [1, 3, 2]);
    
    % Set the permuted tensor back to the Tensor Network
    TN(end).tensor(end).value = last_tensor_permuted;
    % Visualize the TensorNetwork (optional)
    view(TN)
end












% function [E,TN] = example0(P,D)
%     % P = 4;  
%     % D = 10; 
%     % Example: kite with M=5
%      % E = [1 2; 2 3 ;3 4;];
%      E = [1 2; 2 3;3 4;4 5;5 6;];
%     M = max(max(E)); 
% 
%     % compute coordination number of m-th node tensor
%     for m=1:M
%         e(m) = size(find(E==m),1);
%     end   
% 
%     % construct auxiliary matrices(the rows of bd{i} are the edges connected to the i-th node)
%     for i = 1:M
%         [r,~] = find(E == i); % Find all edges connected to this tensor 
%         bd{i} = [E(r,:); i, i]; % Add in the physical dimension
%     end
% 
%     % set up node tensors A(m), m=1,...,M
%         % set up node tensors A(m), m=1,...,M
%     A = [];
%     for m=1:M
%         if m == 1  % check if it is the last site
%             sizevec = [1, D*ones(1, e(m)-1), P];  % dimensions for the first site
%         elseif m == M
%             sizevec = [D*ones(1, e(m)-1),P,1];  % dimensions for the last site
%         else
%             sizevec = [D*ones(1, e(m)-1), D, P];  % dimensions for all other sites
%         end
% 
%         A = [A; Tensor(rand(sizevec), bd{m})];
%     end


    % A = [];
    % for m=1:M
    %     sizevec = [D*ones(1,e(m)),P];  % put physical index last, e.g. for e(m)=3, sizevec = [D D D P]
    %     A = [A; Tensor(rand(sizevec), bd{m})];
    % end
    % 
  