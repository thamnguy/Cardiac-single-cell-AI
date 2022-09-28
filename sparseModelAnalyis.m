function [allCellScore, enrichTable, sparseModel, vioPlot] = ...
    sparseModelAnalyis(expression, cellGroup, positiveGroup, negativeGroup)
% this function perform geneset/pathway/ontology enrichment analysis from 
% a single cell data using sparse model analysis. It require knowing
% 'positive control' cells and 'negative control' cells.
%
% inputs:
%   - expression: the gene expression matrix for each cell. Each row
% correspond to a gene; and the matrix only contains genes belonging to a 
% specific pathway/geneset/gene ontology. Each column correspond to a cell.
% The matrix should contains cells that are 'positive control'
% (expected/known to actively express the pathway/geneset/ontology) and
% 'negative control (opposite of positive control)
%   - cellGroup: the cell group (in name) for each cell, including the positive
% control and negative control groups. This should be string (for easier format)
% and link to the 'positiveGroup' and 'negativeGroup'.
%   - positiveGroup: specify which group name is 'positive control' in the
% cellGroup
%   - negativeGroup: specify which group name is 'negative control' in the
% cellGroup
%
% The sparse analysis will build a model that optimally separate the
% 'positiveGroup' cells from the 'negativeGroup' cells. Then it applies the
% model on other cells. Each cell will have a 'score', which can be used to
% quantify the activity of geneset/pathway/ontology represented by the
% genes used in the expression matrix.
% 
% outputs:
%   - allCellScore: the score for each cell computed by the sparse model
%   - enrichTable: the result of pairwise comparison among all cell groups
%   - sparseModel: the sparse model, include all parameters, to be reused
% if necessary


posIndex = find(ismember(cellGroup, {positiveGroup}) == 1);
negIndex = find(ismember(cellGroup, {negativeGroup}) == 1);
trainMat = [expression(:, posIndex), expression(:, negIndex)];
trainLbl = [ones(length(posIndex), 1); -ones(length(negIndex), 1)];

sparseModel = fitclinear(trainMat', trainLbl, ...
        'Solver','sparsa','Regularization','lasso');

allCellScore = expression' * sparseModel.Beta +  sparseModel.Bias;
figure, vioPlot = violinplot(allCellScore, cellGroup); hold off;

uniqueGroup = unique(cellGroup);

compareNote = {'Overall, Kruskal-Wallis test'};
group1Mean = NaN;
group2Mean = NaN;
pVal = kruskalwallis(allCellScore, cellGroup, 'off');
for i = 1 : length(uniqueGroup)
    for j = (i+1) : length(uniqueGroup)
        compareNote = [ compareNote; [uniqueGroup{i}, ' vs ', uniqueGroup{j}] ];
        group1Mean = [ group1Mean; mean( allCellScore(ismember(cellGroup, uniqueGroup(i))) ) ];
        group2Mean = [ group2Mean; mean( allCellScore(ismember(cellGroup, uniqueGroup(j))) ) ];
        pVal = [ pVal ; ranksum( allCellScore(ismember(cellGroup, uniqueGroup(i))), ...
           allCellScore(ismember(cellGroup, uniqueGroup(j))) ) ];
    end
end

enrichTable = table(compareNote, group1Mean, group2Mean, pVal);

end