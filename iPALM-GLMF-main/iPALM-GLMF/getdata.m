function [Interaction,mSM,dSM]=getdata(path,dataset)
%getdata loads the the interaction information of a particular dataset
%along with the pairwise similarities between the involved microbes and
%diseases
%
% INPUT
%  path:        the path where all the files to be loaded are
%  dataset:     the dataset of interest ('hmdad', 'disbiome')
%
% OUTPUT
%  Interaction: the interaction matrix (diseases rows, microbes columns)
%  dSM:         diseases Similarity Matrix
%  mSM:         microbes Similarity Matrix

    % ================================================================

    %------------------%
    % Adjacency Matrix %
    %------------------%

    % extract the adjacency matrix...
    Interaction = importdata ([ path dataset '_interaction.mat.']);
    % now, 'Interaction' has the matrix WITHOUT the column headers and
    % the row labels.

    % disease Similarity Matrix %
    dSM = getGIPsimilarity(Interaction);  

    % microbe Similarity Matrix %
    mSM = getGIPsimilarity(Interaction');

end