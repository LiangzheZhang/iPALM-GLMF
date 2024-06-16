function [y3,matrixs ]=crossValidation(Y,GM,GD,classifier,n,m,use_WKNKN,K,eta)
%crossValidation runs cross validation experiments
%
% INPUT:
%  Y:           matrix to be modified
%  GM:          microbe similarity matrix  
%  GD:          disease similarity matrix  
%  classifier:  algorithm to be used for Microbe-disease associations 
%               prediction
%  n:           the 'n' in "n-fold experiment"
%  m:           number of known associations in the database
%
% output:
%  y3:          Prediction matrix
%  matrixs:     Final predicted ranking matrix

    % prediction method
    pred_fn = str2func(['alg_' classifier]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% cross validation setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % parameters
    fprintf('n = %i\n',n);      
    y3 = pred_fn(Y,GM,GD,use_WKNKN,K,eta); % predict!
    total_iterations = 100;  % total number of iteration
    iteration = 0;  % current number of iteration
    matrix = zeros(m, total_iterations);
while iteration < total_iterations
    % Randomly divide the 1's in the Y matrix into 5 parts
    num_ones = nnz(Y);
    indices = find(Y == 1);
    indices = indices(randperm(length(indices)));
    part_size = floor(length(indices) / 5);

    % Define containers for tracking cumulative rankings 
    % and occurrences of indexed locations
    index_ranks = containers.Map('KeyType','double','ValueType','double');

    for i = 1:5
        if i == 5
            part_indices = indices((i - 1) * part_size + 1:end);
        else
            part_indices = indices((i - 1) * part_size + 1: i * part_size);
        end
        
        Y_modified = Y;
        for j = 1:length(part_indices)
            index = part_indices(j);
            [row, col] = ind2sub(size(Y), index);
            Y_modified(row, col) = 0;
        end

        G = pred_fn(Y,GM,GD,use_WKNKN,K,eta);
        G_modified = pred_fn(Y_modified,GM,GD,use_WKNKN,K,eta);

        % locate modified element 1
        modified_indices = find(Y ~= Y_modified);

        % get the corresponding G_modified value
        [modified_rows, modified_cols] = ind2sub(size(Y), modified_indices);
        G_modified_values = G_modified(sub2ind(size(G_modified), modified_rows, modified_cols));

        % get the G value corresponding to the element whose Y is 0
        zero_indices = find(Y == 0);
        [zero_rows, zero_cols] = ind2sub(size(Y), zero_indices);
        G_zero_values = G(sub2ind(size(G), zero_rows, zero_cols));

        % compare the ranking of G_modified_values and G_zero_values
        ranks = zeros(length(G_modified_values), 1);
        for k = 1:length(G_modified_values)
            rank = 1;
            for l = 1:length(G_zero_values)
                if G_modified_values(k) < G_zero_values(l)
                    rank = rank + 1;
                end
            end
            ranks(k) = rank;
        end

        % stores the rank of each index position
        for n = 1:length(modified_indices)
            index = modified_indices(n);
            index_ranks(index) = ranks(n);
        end
    end

    disp(['current iteration: ', num2str(iteration + 1), '/', num2str(total_iterations)]);

    % output prediction ranking matrices 
    sorted_ranks = sortrows(cell2mat(index_ranks.keys()), 1);
    ranks = cell2mat(index_ranks.values);
    matrix(:, iteration + 1) = ranks;
    matrixs = matrix';
    iteration = iteration + 1;
end   
end
    
   
        

   

 
    
    

  
