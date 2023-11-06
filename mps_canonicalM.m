function M = mps_canonicalM(mps, N,dir,s)

    
    % Initialize the output MPS
    M = cell(1, N);
    
       switch dir
            case 'left'
                % Sweep from left to right
                for l = 1:N
                    T = mps{l};
                    T = reshape(T, [size(T,1)*size(T,2), size(T,3)]);
                    [U,S,V] = svd(T,'econ');
                    M{l} = reshape(U,[size(U,1)/size(mps{l},2), size(mps{l},2), size(U,2)]);
                    if l < N
                        M{l+1} = tensorprod(S*V',mps{l+1},2, 1);
                    end
                end
    %             for l = 1:N
    %                 T = mps{l};
    %                 [U, S, V] = svd(reshape(T, [size(T, 1)*size(T, 2) ,size(T, 3)]), 'econ');
    %                 M{l} = reshape(U, [size(U,1)/size(mps{l},2), size(mps{l},2), size(U,2)]);
    %                 SV = S * V';
    %                 if l < N
    %                      M{l+1} = tensorprod(SV, mps{l+1}, 2, 1);
    %                 end
    %             end


        case 'right'
            % Sweep from right to left
            for l = N:-1:1
                T = mps{l};
                T = reshape(T, [size(T,1), size(T,2)*size(T,3)]);
                [U,S,V] = svd(T,'econ');
                M{l} = reshape(V',[size(V,2), size(mps{l},2), size(V,1)/size(mps{l},2)]);
                % note V is Hermitian conjugated (c.c. + transpose)
                US = U *S;
                if l > 1
                    M{l-1} = tensorprod(mps{l-1},US,3,1);
                end
            end

      case 'mixed'
%          if s <=1 ||  s == N
%                 error('s must lie between 1 and N')
%            else
            % Sweep from left to right up to site s
            for l = 1:s-1
                T = mps{l};
                [U, S, V] = svd(reshape(T, [size(T, 1)*size(T, 2) ,size(T, 3)]), 'econ');
                M{l} = reshape(U, [size(U,1)/size(mps{l},2), size(mps{l},2), size(U,2)]);
                SV = S * V';
                if l < N
                     M{l+1} = tensorprod(SV, mps{l+1}, 2, 1);
                end
            end
            
                % Sweep from right to left starting from the last site to site s+1
            for l = N:-1:(s+1)
                T = mps{l};
                T = reshape(T, [size(T,1), size(T,2)*size(T,3)]);
                [U,S,V] = svd(T,'econ');
                M{l} = reshape(V',[size(V,2), size(mps{l},2), size(V,1)/size(mps{l},2)]);
                % note V is Hermitian conjugated (c.c. + transpose)
                US = U *S;
                if l > 1
                   M{l-1} = tensorprod(mps{l-1},US,3,1);
                end
           end
        end

        M= mps;
end


































