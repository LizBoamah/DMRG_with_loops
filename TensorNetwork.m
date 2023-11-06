classdef TensorNetwork
    properties
        tensor(:,1) Tensor
        edges(:,2) {mustBeInteger}
    end

    methods
        function obj = TensorNetwork(nodetensors, edges)
            obj.tensor = nodetensors;
            obj.edges = edges;
        end

        function n = num_nodes(TN)
            n = size(TN.tensor, 1);
        end

        function view(TN)
            % Plots the tensor network.
            n = num_nodes(TN);
            figure;
            plot(1:n, zeros(1,n), '.b', 'MarkerSize', 40);
            hold on;
            axis off;
            for i = 1:size(TN.edges,1)
                if abs(TN.edges(i,1) - TN.edges(i,2)) == 1
                    plot(TN.edges(i,:), [0 0], 'b', 'LineWidth', 2);
                else
                    h = @(x) 1 - 1/x;
                    diff = TN.edges(i,2)-TN.edges(i,1);
                    x = linspace(TN.edges(i,1), TN.edges(i,2), 100);
                    plot(x, h(abs(diff)) * sin(pi*(x-x(1))/diff), ...
                        'b', 'LineWidth', 2);
                end
            end
            for i = 1:n
                text(i, 0, num2str(i), 'Color', 'white', ...
                    'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            end
        end
    end
end