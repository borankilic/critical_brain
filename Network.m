classdef Network < handle
    % Network class for handling graph structure and static measures
    % Follows composition design pattern - keeps structure separate from dynamics
    
    properties
        A           % adjacency matrix (sparse NxN)
        N           % number of nodes
        labels      % optional node labels (cell array)
        coords      % optional node coordinates (Nx2 or Nx3 matrix)
        directed = false  % whether network is directed
        weighted = false  % whether network has weighted edges
    end
    
    methods
        function obj = Network(A, varargin)
            % Constructor for Network class
            % Usage: net = Network(A, 'directed', true, 'labels', {...})
            
            if nargin == 0
                return; % allow empty construction
            end
            
            % Set adjacency matrix
            obj.A = sparse(A);
            obj.N = size(A, 1);
            
            % Check if weighted
            obj.weighted = any(obj.A(:) ~= 0 & obj.A(:) ~= 1);
            
            % Parse optional inputs
            p = inputParser;
            addParameter(p, 'directed', false, @islogical);
            addParameter(p, 'labels', {}, @iscell);
            addParameter(p, 'coords', [], @isnumeric);
            parse(p, varargin{:});
            
            obj.directed = p.Results.directed;
            obj.labels = p.Results.labels;
            obj.coords = p.Results.coords;
            
            % Validate inputs
            if size(obj.A, 1) ~= size(obj.A, 2)
                error('Adjacency matrix must be square');
            end
            
            if ~isempty(obj.labels) && length(obj.labels) ~= obj.N
                error('Number of labels must match number of nodes');
            end
            
            if ~isempty(obj.coords) && size(obj.coords, 1) ~= obj.N
                error('Number of coordinate rows must match number of nodes');
            end

            % Assign coordinates if not given
            % if isempty(obj.coords)
            %     G = graph(A);
            %     h = plot(G, 'Layout', 'force');
            %     obj.coords = [h.XData; h.YData]';
            % end

        end
        
        function loadFromFile(obj, filepath, varargin)
            % Load adjacency matrix from file
            % Supports .txt, .csv, .mat files
            
            [~, ~, ext] = fileparts(filepath);
            
            switch lower(ext)
                case '.txt'
                     if contains('CAPA', filepath)
                        % Read the data, skipping the header line
                        % The file uses multiple spaces/tabs as delimiters
                        data = readmatrix(filename, 'NumHeaderLines', 1);
                        
                        % Check if data was read correctly
                        if size(data, 2) < 3
                            error('Data not read correctly. File may have formatting issues.');
                        end
                        
                        % Extract NodeId1, NodeId2, and Weighted Connectivity (column 3)
                        nodeId1 = data(:, 1);
                        nodeId2 = data(:, 2);
                        weightedConnectivity = data(:, 3);
                        
                        % Initialize 800x800 matrix with zeros
                        weightMatrix = zeros(800, 800);
                        
                        % Fill the matrix with weighted connectivity values
                        for i = 1:length(nodeId1)
       
                            if nodeId1(i) >= 1001 && nodeId1(i) <= 1400
                                idx1 = nodeId1(i) - 1000;  % Maps 1001->1, 1400->400
                            elseif nodeId1(i) >= 2001 && nodeId1(i) <= 2400
                                idx1 = nodeId1(i) - 1600;  % Maps 2001->401, 2400->800
                            else
                                idx1 = 0;  % Invalid node
                            end

                            if nodeId2(i) >= 1001 && nodeId2(i) <= 1400
                                idx2 = nodeId2(i) - 1000;  % Maps 1001->1, 1400->400
                            elseif nodeId2(i) >= 2001 && nodeId2(i) <= 2400
                                idx2 = nodeId2(i) - 1600;  % Maps 2001->401, 2400->800
                            else
                                idx2 = 0;  % Invalid node
                            end
                            % Check for valid indices
                            if idx1 > 0 && idx2 > 0
                                % Fill upper triangle
                                weightMatrix(idx1, idx2) = weightedConnectivity(i);
                                
                                % Fill lower triangle (symmetric matrix)
                                weightMatrix(idx2, idx1) = weightedConnectivity(i);
                            end
                        end
                        
                        % Ensure diagonal is zero (should already be zero from initialization)
                        weightMatrix(logical(eye(800))) = 0;
                        
                        % Display matrix properties
                        fprintf('Matrix size: %dx%d\n', size(weightMatrix, 1), size(weightMatrix, 2));
                        fprintf('Number of non-zero elements: %d\n', nnz(weightMatrix));
                        fprintf('Matrix is symmetric: %d\n', issymmetric(weightMatrix));
                        fprintf('Diagonal sum (should be 0): %.10f\n', sum(diag(weightMatrix)));
                        
                        % Optional: Visualize the matrix
                        figure;
                        imagesc(weightMatrix);
                        colorbar;
                        title('Weight Matrix (800x800)');
                        xlabel('Node Index');
                        ylabel('Node Index');
                        axis square;
                        
                        % Optional: Save the matrix
                        save('weight_matrix.mat', 'weightMatrix');


                     end
                     
                case '.csv'
                    obj.A = csvread(filepath);
                case '.mat'
                    data = load(filepath);
                    fields = fieldnames(data);
                    obj.A = data.(fields{1}); % take first variable
                otherwise
                    error('Unsupported file format. Use .txt, .csv, or .mat');
            end
            
            % Update object properties
            obj.A = sparse(obj.A);
            obj.N = size(obj.A, 1);
            obj.weighted = any(obj.A(:) ~= 0 & obj.A(:) ~= 1);
            
            % Parse additional options
            p = inputParser;
            addParameter(p, 'directed', false, @islogical);
            parse(p, varargin{:});
            obj.directed = p.Results.directed;
        end
        
        function filteredMatrix = matrixfilter(A, t)
            % A: NxN adjacency matrix
            % t: threshold factor (e.g. 1 → mean - 1*std)
            % A_thresh: thresholded adjacency matrix
            
            mu = mean(A(:));        % mean of all entries
            sigma = std(A(:));      % standard deviation of all entries
            threshold = mu - t * sigma;
            
            filteredMatrix = A;           % copy matrix
            for i = 1:size(A,1)
                for j = 1:size(A,2)
                    if A(i,j) < threshold
                        filteredMatrix(i,j) = 0;
                    end
                end
            end
        end

        function normalize_weights(obj)
            row_sums = sum(obj.A, 2);
            row_sums(row_sums == 0) = 1;  % avoid divide by zero
            obj.A = obj.A ./ row_sums;
        end
        function L = laplacian(obj)
            % Compute graph Laplacian matrix
            d = obj.degree();
            L = spdiags(d, 0, obj.N, obj.N) - obj.A;
        end
        
        function d = degree(obj)
            % Compute node degrees
            if obj.directed
                d.in = sum(obj.A, 1)';  % in-degree
                d.out = sum(obj.A, 2);  % out-degree
                d.total = d.in + d.out;
            else
                d = sum(obj.A, 2);
            end
        end
        
        function C = clusteringCoefficients(obj)
            % Compute local clustering coefficients
            C = zeros(obj.N, 1);
            
            for i = 1:obj.N
                neighbors = find(obj.A(i, :));
                k = length(neighbors);
                
                if k < 2
                    C(i) = 0;
                    continue;
                end
                
                % Count triangles
                subgraph = obj.A(neighbors, neighbors);
                triangles = sum(sum(subgraph)) / 2;
                
                % Clustering coefficient
                if obj.directed
                    possible = k * (k - 1);
                else
                    possible = k * (k - 1) / 2;
                end
                
                C(i) = triangles / possible;
            end
        end
        
        function bc = betweennessCentrality(obj)
            % Compute betweenness centrality (simplified version)
            % For large networks, consider using specialized toolboxes
            
            bc = zeros(obj.N, 1);
            
            % Simple implementation - for better performance use specialized algorithms
            for s = 1:obj.N
                for t = s+1:obj.N
                    if s == t, continue; end
                    
                    % Find shortest paths using Dijkstra (simplified)
                    paths = obj.findShortestPaths(s, t);
                    if ~isempty(paths)
                        for p = 1:length(paths)
                            path = paths{p};
                            % Add to betweenness for intermediate nodes
                            for node = path(2:end-1)
                                bc(node) = bc(node) + 1/length(paths);
                            end
                        end
                    end
                end
            end
            
            % Normalize
            if obj.N > 2
                if obj.directed
                    bc = bc / ((obj.N-1)*(obj.N-2));
                else
                    bc = bc / ((obj.N-1)*(obj.N-2)/2);
                end
            end
        end
        
        function stats = getNetworkStats(obj)
            % Get comprehensive network statistics
            d = obj.degree();
            % Convert sparse matrix into full for stats calculatiom
            d_full = full(d);
            if isstruct(d_full)
                stats.degree_mean = mean(d_full.total);
                stats.degree_std = std(d_full.total);
            else
                stats.degree_mean = mean(d_full);
                stats.degree_std = std(d_full);
            end
            
            stats.density = nnz(obj.A) / (obj.N * (obj.N - 1));
            stats.num_nodes = obj.N;
            stats.num_edges = nnz(obj.A) / (1 + ~obj.directed); % adjust for undirected
            stats.directed = obj.directed;
            stats.weighted = obj.weighted;
            
            % Clustering
            C = obj.clusteringCoefficients();
            stats.clustering_mean = mean(C);
            stats.clustering_std = std(C);
        end
        
        function save(obj, filepath)
            % Save network to .mat file with metadata
            A = obj.A; %#ok<NASGU>
            N = obj.N; %#ok<NASGU>
            directed = obj.directed; %#ok<NASGU>
            weighted = obj.weighted; %#ok<NASGU>
            labels = obj.labels; %#ok<NASGU>
            coords = obj.coords; %#ok<NASGU>
            
            save(filepath, 'A', 'N', 'directed', 'weighted', 'labels', 'coords');
            fprintf('Network saved to %s\n', filepath);
        end
        
        



        function visualize(obj, varargin)
            % Create input parser with enhanced options
            p = inputParser;
            addParameter(p, 'min_linewidth', 0.5, @isnumeric);
            addParameter(p, 'max_linewidth', 4, @isnumeric);
            addParameter(p, 'node_size', 5, @isnumeric);
            addParameter(p, 'node_color', [0.2 0.4 0.8], @(x) isnumeric(x) && numel(x) == 3);
            addParameter(p, 'edge_color', [0.7 0.7 0.7], @(x) isnumeric(x) && numel(x) == 3);
            addParameter(p, 'layout', 'force', @(x) ismember(x, {'force', 'circle', 'layered', 'subspace', 'auto'}));
            addParameter(p, 'show_labels', true, @islogical);
            addParameter(p, 'font_size', 8, @isnumeric);
            addParameter(p, 'color_edges_by_weight', false, @islogical);
            addParameter(p, 'colormap', 'parula', @ischar);
            parse(p, varargin{:});
            
            % Get parsed parameters
            results = p.Results;
            
            % Create figure
            figure('Color', 'white', 'Position', [100, 100, 800, 600]);
            
            % Create graph object from adjacency matrix
            if isa(obj.A, 'logical') || all(obj.A(:) == 0 | obj.A(:) == 1)
                if ~isempty(obj.labels)
                    G = graph(obj.A, obj.labels);
                else 
                    G = graph(obj.A);
                end
            else
                if ~isempty(obj.labels)
                    G = graph(obj.A, obj.labels, 'omitselfloops');
                else 
                    G = graph(obj.A, 'omitselfloops');
                end
            end
            
            % Get edge weights for scaling
            edge_weights = G.Edges.Weight;
            if isempty(edge_weights)
                edge_weights = ones(numedges(G), 1);
            end
            
            % Calculate scaled line widths
            min_weight = min(edge_weights);
            max_weight = max(edge_weights);
            
            if max_weight > min_weight
                linewidths = results.min_linewidth + (edge_weights - min_weight) * ...
                            (results.max_linewidth - results.min_linewidth) / (max_weight - min_weight);
            else
                linewidths = repmat((results.min_linewidth + results.max_linewidth) / 2, numedges(G), 1);
            end
            
            % Choose layout method
            if ~isempty(obj.coords) && size(obj.coords, 1) == numnodes(G)
                % Use provided coordinates
                x = obj.coords(:, 1);
                y = obj.coords(:, 2);
                if size(obj.coords, 2) == 3
                    % 3D coordinates
                    h = plot(G, 'XData', x, 'YData', y, 'ZData', obj.coords(:, 3));
                else
                    % 2D coordinates
                    h = plot(G, 'XData', x, 'YData', y);
                end
            else
                % Use automatic layout
                switch results.layout
                    case 'force'
                        h = plot(G, 'Layout', 'force', 'Iterations', 100, ...
                                'UseGravity', true, 'WeightEffect', 'inverse');
                    case 'circle'
                        h = plot(G, 'Layout', 'circle');
                    case 'layered'
                        h = plot(G, 'Layout', 'layered');
                    case 'subspace'
                        h = plot(G, 'Layout', 'subspace', 'Dimension', 2);
                    otherwise
                        h = plot(G);
                end
            end
        
            % Customize appearance
            hold on;
            
            % Set node properties
            h.NodeColor = results.node_color;
            h.MarkerSize = results.node_size;
            h.Marker = 'o';
            h.NodeFontSize = results.font_size;
            h.NodeLabel = {};
            
            % Set edge properties
            h.LineWidth = linewidths;
            h.EdgeColor = results.edge_color;
            
            % Color edges by weight if requested
            if results.color_edges_by_weight && max_weight > min_weight
                colormap(results.colormap);
                normalized_weights = (edge_weights - min_weight) / (max_weight - min_weight);
                h.EdgeCData = normalized_weights;
                colorbar;
                title('Edge Weights');
            end
            
            % Add node labels if requested
            if results.show_labels && ~isempty(obj.labels)
                if size(obj.coords, 2) == 3
                    % 3D text labels
                    for i = 1:numnodes(G)
                        text(h.XData(i), h.YData(i), h.ZData(i), obj.labels{i}, ...
                            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                            'FontSize', results.font_size, 'BackgroundColor', 'white', ...
                            'Margin', 0.5, 'FontWeight', 'bold');
                    end
                else
                    % 2D text labels
                    for i = 1:numnodes(G)
                        text(h.XData(i), h.YData(i), obj.labels{i}, ...
                            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                            'FontSize', results.font_size, 'BackgroundColor', 'white', ...
                            'Margin', 0.5, 'FontWeight', 'bold');
                    end
                end
            end
            
            % Enhance plot appearance
            title(sprintf('Network Visualization (%d nodes, %d edges)', numnodes(G), numedges(G)));
            grid on;
            axis equal;
            set(gca, 'Box', 'on', 'FontSize', 10);
            
            % Add interactive features
            zoom on;
            rotate3d on;
            
            % Add legend for edge weights
            if max_weight > min_weight
                annotation('textbox', [0.02, 0.02, 0.2, 0.05], 'String', ...
                    sprintf('Edge weights: %.2f - %.2f', min_weight, max_weight), ...
                    'FitBoxToText', 'on', 'BackgroundColor', 'white');
            end
            
            hold off;
        end
    end
    methods (Access = private)
        function paths = findShortestPaths(obj, source, target)
            % Simple shortest path finder (for betweenness centrality)
            % This is a simplified implementation
            paths = {};
            
            % Use breadth-first search for unweighted graphs
            queue = {source};
            visited = false(obj.N, 1);
            parent = zeros(obj.N, 1);
            
            visited(source) = true;
            
            while ~isempty(queue)
                current = queue{1};
                queue(1) = [];
                
                if current == target
                    % Reconstruct path
                    path = target;
                    while parent(path(1)) ~= 0
                        path = [parent(path(1)), path]; %#ok<AGROW>
                    end
                    paths{end+1} = path; %#ok<AGROW>
                    return;
                end
                
                neighbors = find(obj.A(current, :));
                for neighbor = neighbors
                    if ~visited(neighbor)
                        visited(neighbor) = true;
                        parent(neighbor) = current;
                        queue{end+1} = neighbor; %#ok<AGROW>
                    end
                end
            end
        end
    end
    
    methods (Static)
        function net = generateRandomNetwork(N, p, varargin)
            % Generate random Erdős–Rényi network
            % Usage: net = Network.generateRandomNetwork(100, 0.1, 'directed', false)
            
            p_parser = inputParser;
            addParameter(p_parser, 'directed', false, @islogical);
            parse(p_parser, varargin{:});
            
            if p_parser.Results.directed
                A = rand(N, N) < p;
                A = A - diag(diag(A)); % remove self-loops
            else
                A = triu(rand(N, N) < p, 1);
                A = A + A'; % make symmetric
            end
            
            net = Network(A, 'directed', p_parser.Results.directed);
        end
        
        function net = generateSmallWorldNetwork(N, k, beta)
            % Generate Watts-Strogatz small-world network
            % N: number of nodes, k: initial degree, beta: rewiring probability
            
            % Start with regular ring lattice
            A = zeros(N, N);
            for i = 1:N
                for j = 1:k/2
                    target1 = mod(i + j - 1, N) + 1;
                    target2 = mod(i - j - 1, N) + 1;
                    A(i, target1) = 1;
                    A(i, target2) = 1;
                end
            end
            
            % Rewire edges
            [rows, cols] = find(triu(A, 1));
            for i = 1:length(rows)
                if rand < beta
                    % Remove old edge and add new random edge
                    A(rows(i), cols(i)) = 0;
                    A(cols(i), rows(i)) = 0;
                    
                    % Find new target (avoiding self-loops and existing edges)
                    possible_targets = setdiff(1:N, [rows(i), find(A(rows(i), :))]);
                    if ~isempty(possible_targets)
                        new_target = possible_targets(randi(length(possible_targets)));
                        A(rows(i), new_target) = 1;
                        A(new_target, rows(i)) = 1;
                    end
                end
            end
            
            net = Network(A);
        end
    end
end