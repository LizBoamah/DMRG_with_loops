classdef Tensor
    properties
        value {mustBeNumeric}
        dims (:,2) {mustBeInteger, mustBeSorted}
    end
    
    methods
        function obj = Tensor(v, d)
            % Check to make sure v and d match.
            if isscalar(v)
                assert(size(d,1) == 0, ['Tensor is a scalar, but has ', ...
                    int2str(size(d,1)), ' desciptions for dimensions.']);
            elseif iscolumn(v)
                assert(size(d,1) == 1, ['Tensor is a vector, but has ', ...
                    'descriptions for ', int2str(size(d,1)), ' dimensions.']);
            else
                assert(ndims(v) == size(d,1), ...
                    ['The tensor has ', int2str(ndims(v)), ...
                    ' dimensions, but there are desciptions for ', ...
                    int2str(size(d,1)), ' dimensions.']);
            end

            obj.value = v;
            obj.dims = sort(d,2);
        end

        function out = tensorprod(T1, T2, d1, d2)
            % Compute the tensor product.
            val = tensorprod(T1.value, T2.value, d1, d2);
            % Determine which dimensions are left.
            T1_remaining = true(size(T1.dims, 1), 1);
            T1_remaining(d1) = false;
            T2_remaining = true(size(T2.dims, 1), 1);
            T2_remaining(d2) = false;
            % Construct the output tensor.
            out = Tensor(val, ...
                vertcat(T1.dims(T1_remaining,:), T2.dims(T2_remaining,:)));
        end

        function n = num_dims(T)
            n =size(T.dims,1);
        end

        function Tf = flip_dims(T)
            % If T is a vector, do nothing.
            if isvector(T.value)
                Tf = T;
            else
            % Otherwise, reverse the order of the dimensions of the tensor.
                Tf = Tensor(permute(T.value, flip(1:ndims(T.value))), ...
                    flip(T.dims, 1));
            end
        end
        
        function idx = find_dim(T, d)
            assert(all(size(d) == [2,1]) || all(size(d) == [1,2]));
            [idx1, ~] = find(T.dims(:,1) == d(1));
            [idx2, ~] = find(T.dims(:,2) == d(2));
            idx = intersect(idx1, idx2);
        end
    end
end

function mustBeSorted(a)
    assert(issorted(a,2), 'Smaller number must come first in edges.');
end



% classdef Tensor
%     properties
%         value {mustBeNumeric}
%         dims (:,:) {mustBeInteger}
%     end
%     
%     methods
%         function obj = Tensor(v, d)
%             % Check to make sure v and d match.
%             if numel(v) == 1
%                 assert(isempty(d), ['Tensor is a scalar, but has ', ...
%                     int2str(size(d,2)), ' descriptions for dimensions.']);
%             else
%                 assert(ndims(v) == size(d,2), ...
%                     ['The tensor has ', int2str(ndims(v)), ...
%                     ' dimensions, but there are descriptions for ', ...
%                     int2str(size(d,2)), ' dimensions.']);
%             end
% 
%             obj.value = v;
%             obj.dims = sort(d,2);
%         end
% 
%         function out = tensorprod(T1, T2, d1, d2)
%             % Compute the tensor product.
%             val = tensorprod(T1.value, T2.value, d1, d2);
%             % Determine which dimensions are left.
%             T1_remaining = true(size(T1.dims, 2), 1);
%             T1_remaining(d1) = false;
%             T2_remaining = true(size(T2.dims, 2), 1);
%             T2_remaining(d2) = false;
%             % Construct the output tensor.
%             out = Tensor(val, ...
%                 [T1.dims(:,T1_remaining), T2.dims(:,T2_remaining)]);
%         end
% 
%         function n = num_dims(T)
%             n = size(T.dims,2);
%         end
% 
%         function Tf = flip_dims(T)
%             % If T is a scalar, do nothing.
%             if numel(T.value) == 1
%                 Tf = T;
%             else
%             % Otherwise, reverse the order of the dimensions of the tensor.
%                 Tf = Tensor(permute(T.value, flip(1:ndims(T.value))), ...
%                     fliplr(T.dims));
%             end
%         end
%         
%         function idx = find_dim(T, d)
%             assert(isvector(d));
%             idx = all(T.dims == d, 1);
%         end
%     end
% end
% 
% 
% 
% 
% 
% 
% 




























