function Y=preprocess_WKNKN(Y,GM,GD,K,eta)
%preprocess_WKNKN preprocesses the interaction matrix Y by replacing each
%of the 0's (i.e. presumed non-interactions) with a continuous value
%between 0 and 1. For each 0, the K nearest known drugs are used to infer
%a value, the K nearest known targets are used to infer another value, and
%then the average of the two values is used to replace that 0.
%
%   Y=preprocess_WKNKN(Y,GM,GD,K,eta)
% INPUT:
%  Y:   matrix to be modified
%  Sd:  pairwise row similarities matrix
%  St:  pairwise column similarities matrix
%  K:   number of nearest known neighbors to use
%  eta: decay rate
%
% OUTPUT:
%  Y:   the modified matrix

    % decay values to be used in weighting similarities later
    eta = eta .^ (0:K-1);

    y2_new1 = zeros(size(Y));
    y2_new2 = zeros(size(Y));

    empty_rows = find(any(Y,2) == 0);   % get indices of empty rows
    empty_cols = find(any(Y)   == 0);   % get indices of empty columns

    % for each microbe i...
    for i=1:length(GM)
        microbe_sim = GM(i,:); % get similarities of microbe i to other 
                               % microbes
        microbe_sim(i) = 0;    % set self-similiraty to ZERO

        indices  = 1:length(GM);    % ignore similarities 
        microbe_sim(empty_rows) = [];  % to microbes of 
        indices(empty_rows) = [];   % empty rows

        [~,indx] = sort(microbe_sim,'descend');    % sort descendingly
        indx = indx(1:K);       % keep only similarities of K nearest neighbors
        indx = indices(indx);   % and their indices

        % computed profile of microbe i by using its similarities to its K
        % nearest neighbors weighted by the decay values from eta
        microbe_sim = GM(i,:);
        y2_new1(i,:) = (eta .* microbe_sim(indx)) * Y(indx,:) ./ sum(microbe_sim(indx));
    end

    % for each disease j...
    for j=1:length(GD)
        disease_sim = GD(j,:); % get similarities of target j to other diseases
        disease_sim(j) = 0;    % set self-similiraty to ZERO

        indices  = 1:length(GD);        % ignore similarities 
        disease_sim(empty_cols) = [];    % to diseases of
        indices(empty_cols) = [];       % empty columns

        [~,indx] = sort(disease_sim,'descend');  % sort descendingly
        indx = indx(1:K);       % keep only similarities of K nearest neighbors
        indx = indices(indx);   % and their indices

        % computed profile of disease j by using its similarities to its K
        % nearest neighbors weighted by the decay values from eta
        disease_sim = GD(j,:);
        y2_new2(:,j) = Y(:,indx) * (eta .* disease_sim(indx))' ./ sum(disease_sim(indx));
    end

    % average computed values of the modified 0's from the microbe and
    % disease sides while preserving the 1's that were already in Y 
    Y = max(Y,(y2_new1 + y2_new2)/2);

end