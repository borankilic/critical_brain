r1 = 0.001;
r2 = 0.64;
threshold = 0.00;

T= 5000;
% 25x25 grid test
n = 800;
node_number = n ;  % 625 nodes

% Create grid connectivity (each node connected to 4 neighbors: up, down, left, right)

load("weight_matrix.mat");
weight_matrix = weightMatrix;

connectivity = (weight_matrix >0);
net = Network(weight_matrix);
net.normalize_weights;

for i = 1:10
    for j = 1:120
        params = GreenbergHastingsModel.getDefaultParams();
        params.r2 = r2;
        params.r1=r1;
        params.threshold=threshold;
        params.initial_excited = 0;
        ghModel = GreenbergHastingsModel(net, params);
        ghModel.run(T)
        largest_sizes = zeros(1, T);
        second_largest_sizes = zeros(1, T);
        for t = 1:T
            [~, number_of_clusters, cluster_sizes] = ghModel.find_cluster(t);
            
            if number_of_clusters > 0
                top_three = get_top_three_cluster_sizes(cluster_sizes);
                largest_sizes(t) = top_three(1);
                second_largest_sizes(t) = top_three(2);
                
                %disp(number_of_clusters);
                %disp(top_three);
            end
        end

        fid = fopen('experiment.txt', 'a');
        avg_largest = mean(largest_sizes);
        avg_second_largest = mean(second_largest_sizes);
        fprintf(fid, '%f %f %f %f\n', avg_largest, avg_second_largest, r2, threshold);
        fprintf(fid, '\n');  % blank line
        
        fclose(fid);

        threshold = 0.000 + j * 0.002;
    end
    r2 = 0.64 + i * 0.03; 
end



function top_three_sizes = get_top_three_cluster_sizes(cluster_sizes)
    if isempty(cluster_sizes)
        top_three_sizes = [0, 0, 0];
        return;
    end
    
    % Sort in descending order
    sorted_sizes = sort(cluster_sizes, 'descend');
    
    % Get top 3 (pad with zeros if fewer than 3 clusters)
    top_three_sizes = [sorted_sizes, zeros(1, 3)];
    top_three_sizes = top_three_sizes(1:3);
end