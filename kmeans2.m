function [ means, clusters ] = kmeans2( data, n )
%kmeans2 - Get matrix of means, matrix of clusters for n
%Implmentation by Nikolay Bogoychev, s1031254

distance = get_distance_matrix(data);

%Get the means:
[means_init] = get_means(distance, n);

%Get the inital cluster for the means. Use it later on in case we have 2
%points that could go into either cluster.
dummy_clusters = zeros(length(data(:,1)));

[clusters_init]=get_clusters(means_init, distance,dummy_clusters, n);

old_means = zeros(n, 3);
for i = 1:n
    old_means(i, :) = data(means_init(i), :);
end   

[means,clusters] = reclusterize(data,distance,clusters_init, old_means);

end

function [clusters] = get_clusters(means, distance,old_clusters, n )
%get_clusters Get the clusters
%mean is actually just an index here
%In case of even distances we'd use old clusters to assign the the point to
%the cluster it was previously in

clusters = zeros(length(distance(1, :)),1); %We check the distance for each point with all the means and
                %assign it to the smallest
for i = 1:(length(distance(1,:)))
    smallest = Inf; %Use Infinity so that we don't get a distance that's bigger than it
    smallest_index = 0; % Check which one it is
    for j = 1:n
        if distance(means(j), i) < smallest
            smallest = distance(means(j), i);
            smallest_index = j;
        elseif distance(means(j), i) == smallest
            smallest_index = old_clusters(i);
        end
    end
    clusters(i) = smallest_index;
end
end

function [ distance, new_means, old_means ] = recompute_mean(data,distance, clusters,old_means)
%recompute_mean Recompute the means for the clusters
n = length(old_means(:,1));
new_means = zeros(n,3);
%For each cluster, sum it up and divide it by the number of elements to get
%the new mean
for i = 1:n
    sum_clust = 0;
    for j = 1:(length(clusters))
        if clusters(j)==i
            sum_clust = sum_clust + data(j, :);
        end
    end
    new_mean = sum_clust/(sum(clusters ==i));
    new_means(i,:) = new_mean;
end

%Because we use the index of the distance matrix to calculate the distance
%between every point and the centers, after we create new centers we just
%extend the distance matrix with more entries which correspond to the d
%distances between the centers and the points.

%Get how long is the current list
extend_index = length(distance(:,1));
for j = 1:n %get how many elements exactly are we extending with
    for i = 1:(length(distance(1, :)))
        distance(extend_index + j,i) = sqrt(sum((data(i,: ) - new_means(j, :)).^2));
    end
end
end

function [ new_means, clusters ] = reclusterize( data, distance,old_clusters, old_means )
%reclusterize Reccursively compute the final cluster

%Recompute the new means. if they are the same as the old ones, terminate
[distance, new_means, old_means] = recompute_mean(data, distance, old_clusters, old_means);
if (all(new_means == old_means))
    done = 1; %indicate end
else
    done = 0;
end

%Get the new indexes for the means. 
%After the first n elements, the indexes of the next n elements and so on
%are just the last n rows of the distances matrix, as computed by the 
%recompute_means function

n = length(new_means(:,1)); 
means_indexes = zeros(n,1);
for i = 1:n
    means_indexes(i) = length(distance) - (n - i);
end

[clusters]=get_clusters(means_indexes, distance, old_clusters, n);

%recursion if we're not done.
if ~done
    [ new_means, clusters ] = reclusterize(data, distance, clusters, new_means);
end
end

function [all_means] = get_means(distance, n)
%get_means get the indeces of all means.
[first,second] = ind2sub(size(distance),find(distance==max(max(distance))));
%Get the index from the maximum distance. Note that the entries are
%duplicated so we have to pick first(1) and second(1) or first(2) and first(2).
%In case of more than one value that turns out to be the max, those would pick the
%first found the loop computes the rest of the means, up to n 

all_means = zeros(n,1); %Container for all the measns
for i = 1:n
    if i == 1
        all_means(i) = first(1);
    elseif i == 2
        all_means(i) = second(1);
    else
        max_dist_temp = 0; %Maximum distance
        next_index = 0; %Index of the next mean
        for j = 1:length(distance(1, :))
            temp = 1;
            for n = 1:i %Calculate the biggest distance
                if all_means(n) ~= 0
                    temp = temp*distance(all_means(n),j);
                end
            end
            temp = temp^2; %Get abs value
            temp = sqrt(temp); 
            if temp > max_dist_temp
                max_dist_temp = temp;
                next_index = j;
            end
        end
        all_means(n) = next_index;
    end
end  
end

function [ distance ] = get_distance_matrix( input_matrix)
%get_distance_matrix Get the distance matrix, given the data file
mat_length = (length(input_matrix));
distance = zeros(mat_length,mat_length);
for i = 1:mat_length
    for j = 1:mat_length
        distance(i,j) = sqrt(sum((input_matrix(i,: ) - input_matrix(j,: )).^2));
    end
end
end
