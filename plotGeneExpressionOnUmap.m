function [newFig] = plotGeneExpressionOnUmap(markerExpression, coordinate, threshold)
% plot a specific gene expression on each cell on the top of a 2D cell
% input:
% - markerExpression: vector of gene expression in all cells
% - coordinate: 2D coordiate of each cells (i.e, from UMAP or TSNE plots)
% - threshold: the minimum expression threshold that will be colored. If
% the expression in the cell is less than the threshold, the cell will
% appear as a light grey dot.

figure
hold on
index1 = find(markerExpression <= threshold);
scatter(coordinate(index1, 1), coordinate(index1, 2), 15, [192, 192, 192]/255, '.');
index1 = find(markerExpression > threshold);
scatter(coordinate(index1, 1), coordinate(index1, 2), 15, markerExpression(index1), '.');
colorbar;
colormap('jet');
xlabel('umap 1'); ylabel('umap 2');
set(gca,'FontSize',16)
set(gca, 'linewidth', 1.5, 'XColor', 'k', 'Ycolor', 'k', 'TickDir', 'out', 'Box', 'off')
hold off

newFig = gca;
end