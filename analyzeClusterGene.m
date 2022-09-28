function [meanExp, percenExp, logFoldChange, pFisher, oddRatio, pRanksum] ...
    = analyzeClusterGene(expressMat, clusterID, isLogFormat)
% function to compute each gene statistics in all clusters of the
% single-cell data
%
% inputs:
%   - expressMat: expression matrix, assume to be in log-normalization
% (isLogFormat = true), which is done by Seurat. The matrix is prefered to
% be in sparse format
% (https://www.mathworks.com/help/matlab/sparse-matrices.html) for memory
% efficiency. Each row associates to a gene, and each column associates
% to a cell.
%   - clusterID: vector of cluster number assigned for each cell. The
% clusterID should be consecutive integers from 1 to the number of cluster
% founded in the data.
%   - isLogFormat: indicate whether the expressMat is being log-(base 2)
% normalized. isLogFormat = false means it isn't in log-base 2.
%
% outputs:
%   - meanExp: table (matrix) of mean gene expression for each gene
% (associated to a row) in each cluster (associated to a column)
%   - percenExp: table (matrix) of how many percentage of cells in each
% cluster (associated to a column) expressing the gene (row).
%   - logFoldChange: table (matrix) of fold-change (in base 2 logarithm)
% of the gene (row) expression, between cells belonging to the cluster
% (column) and not belonging to the cluster.
%   - pFisher:  table (matrix) of p-Value done by Fisher's exact test,
% which can be use as a statistical metric indicating whether a gene (row)
% is a marker for a cluster (column)
%   - oddRatio: table (matrix) of odd-ratio done by Fisher's exact test
%   - pRanksum:  table (matrix) of p-Value done by Wilcoxon-Ranksum test,
% which can be use as a statistical metric indicating whether a gene (row)
% is a marker for a cluster (column)

if isLogFormat == true
    expressMat = 2.^expressMat - 1;
end

meanExp = zeros(size(expressMat, 1), max(clusterID));
percenExp = zeros(size(expressMat, 1), max(clusterID));
logFoldChange = zeros(size(expressMat, 1), max(clusterID));
pFisher = ones(size(expressMat, 1), max(clusterID));
oddRatio = ones(size(expressMat, 1), max(clusterID));
pRanksum = ones(size(expressMat, 1), max(clusterID));

for j = 1 : max(clusterID)
    clusterIndex = find(clusterID == j);
    nonClusterIndex = find(clusterID ~= j);
    meanExp(:, j) = mean( expressMat(:, clusterIndex)' )';
    percenExp(:, j) = sum( sign(expressMat(:, clusterIndex))' )' / length(clusterIndex);
    logFoldChange(:, j) = log2(mean( expressMat(:, clusterIndex)' )' ...
        ./ mean( expressMat(:, nonClusterIndex)' )');

    for i = 1 : size(expressMat, 1)
        count1 = length ( find( expressMat(i, nonClusterIndex) == 0 ) );
        count2 = length ( find( expressMat(i, clusterIndex) == 0 ) );
        count3 = length ( find( expressMat(i, nonClusterIndex) > 0 ) );
        count4 = length ( find( expressMat(i, clusterIndex) > 0 ) );
        countNum = [ count1, count2; count3, count4];

        [h,pFisher(i, j),stats] = fishertest(countNum,'Tail','right','Alpha',0.01);
        oddRatio(i, j) = stats.OddsRatio;

        pRanksum(i, j) = ranksum(expressMat(i, clusterIndex), expressMat(i, nonClusterIndex));
    end
end

end