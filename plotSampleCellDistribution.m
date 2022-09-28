function [allFigs] = plotSampleCellDistribution (sampleID, coordinate, clusterID)
% plot the sample cell distribution, include:
%   - localization of a group cell (specified by sampleID) on UMAP. Each
% group has it owns figure
%   - stacked bar, proportion of each cell cluster in each group

uniqueSampleID = unique(sampleID);

allFigs = cell( length(unique(sampleID)) + 1, 1 );
for i = 1 : length(uniqueSampleID)
    allFigs{i} = figure;
    index1 = find( ~ismember(sampleID, uniqueSampleID(i)) );
    scatter(coordinate(index1, 1), coordinate(index1, 2), 15, [192, 192, 192]/255, '.');
    hold on
    index1 = find( ismember(sampleID, uniqueSampleID(i)) );
    scatter(coordinate(index1, 1), coordinate(index1, 2), 30, 'r', '.');
    legend( 'Other sample', uniqueSampleID{i} );
    xlabel('umap 1'); ylabel('umap 2');
    % title([ marker,  ' > ', num2str(threshold)]);
    set(gca,'FontSize',16)
    % xlim([-13 13])
    % ylim([-20 16])
    hold off
    set(gca, 'Box', 'off');
end

allFigs{ length(unique(sampleID)) + 1 } = figure;
countCell = zeros(max(clusterID), length(uniqueSampleID));
for i = 1 : max(clusterID)
    for j = 1 : length(uniqueSampleID)
        countCell(i, j) = length( find( clusterID==i & ismember(sampleID, uniqueSampleID(j)) ) );
    end
end

ratioMat = zeros(size(countCell));
for i = 1 :  length(uniqueSampleID)
    ratioMat(:, i) = countCell(:, i) / sum(countCell(:, i));
end

bar(ratioMat','stacked');
legendLbl = cell(max(clusterID), 1);
for i = 1 : max(clusterID)
    legendLbl{i} = ['Cluster', num2str(i)];
end
legend(legendLbl);
xticklabels(uniqueSampleID);