function similarityMatrix = getGIPsimilarity(adj_matrix)
    % weighting
    widthSum = 0;
    
    % calculate the constant r that controls the generational width of 
    % the Gaussian kernel
    rows = size(adj_matrix, 1);
    for i = 1:rows
        width = norm(adj_matrix(i, :)); % find the biparadigm of a row vector
        widthSum = widthSum + width^2;
    end
    r = rows / widthSum;

    % find each element of the Gaussian kernel matrix
    similarityMatrix = zeros(rows, rows);
    for i = 1:rows
        for j = i:rows
            ip = norm(adj_matrix(i, :) - adj_matrix(j, :))^2;
            similarityMatrix(i, j) = exp(-r * ip);
            similarityMatrix(j, i) = similarityMatrix(i, j);
        end
    end
end
