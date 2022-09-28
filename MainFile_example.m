%% Compute the autoencoder to cluster all cells (cell type identification

% get the raw transcript count
x = importdata('AdjustedCount.txt');
expressMat = spconvert(x);
save AdjustedCount.mat expressMat -mat; % optional, save the result for future usage
sumExpress = full( sum( expressMat, 2 ) ); % total expression (raw transcript
% count of each gene), which we will use to select the highest
% expressing genes to compute the autoencoder

% get the normalized expression matrix
x = importdata('NormalizedExpression.txt');
expressMat = spconvert(x);
save NormalizedExpression.mat expressMat -mat; % optional, save the result for future usage
autoencoderMat = full( expressMat( find(sumExpress > 40000) , : ) );
% this example use the genes with total expression of 40000 UMI and more to
% build the autoencoder. This threshold can change case-by-case. Less
% threshold means more genes will be used, and requires longer computing
% time

autoenc = trainAutoencoder(autoencoderMat, 'UseGPU', true);
% if the autoencoderMat is small enough to fit in the graphic card (GPU)
% memory, then set 'UseGPU' to true will drastically reduce computing time.
% Otherwise, we must set it to false.

cellEncode = encode(autoenc, autoencoderMat); %use autoencoder to embed the cells

save autoencoderNet.mat autoenc -mat; % optional, save the autoencoder for reusing
save cellEncode.mat cellEncode -mat;% optional, save the cell embeding for reusing

%% compute UMAP and clustering
% load the cell sample ID, which is also the cell groups in this analysis
sampleID = readtable('sampleID.csv');
sampleID = table2cell(sampleID(:, 1));
for i = 1 : length(sampleID)
    tempStr = split(sampleID{i}, '_'); % the cell ID format is <groupID>_<cell barcode>
    % so we split it by the underscore (_)
    sampleID{i} = tempStr{1};
end

addpath(genpath([pwd, '/MatlabUmap'])); % load the MatlabUmap, or simply
% all folders and subfolders in MatlabUmap into Matlab Path.
% See https://www.mathworks.com/help/matlab/ref/addpath.html

[umapCoor, umapStruct, umapCluster] = run_umap(cellEncode');
% umapCoor: umap the 2D embedding (coordinate)
% umapStruct: umap data structure, usually not being used
% umapCluster: the umap Matlab implementation can also perform cell
% clustering. In some cases, this result is useful.

save umapCoor.mat umapCoor -mat; % optional, save the umap 2D embedding for reusing
save umapCluster.mat umapCluster -mat; % optional, save the umap cluster result for reusing

% density-based clustering (dbscan) is more prefered to cluster the cells
% see more at https://www.mathworks.com/help/stats/dbscan-clustering.html
kD = pdist2(umapCoor,umapCoor,'euc','Smallest',50);
figure, plot(sort(kD(end,:))); % the 'curve turning area' in this figure tells how to set parameter in dbscan
clusterID = dbscan(umapCoor, 0.6, 50); % 0.6 seems to be at the end of the turning area. Other option is choosing the turning point.
figure, gscatter(umapCoor(:, 1), umapCoor(:, 2), clusterID);

% compute the statisitics for each gene in each cluster
geneTable = readtable('filterGeneList.csv');
gene = table2cell(geneTable(:, 2)); % load the gene list
[meanExp, percenExp, logFoldChange, pFisher, oddRatio, pRanksum] = ...
    analyzeClusterGene(expressMat, clusterID, true);
% so far the expression matrix is still under log base 2 normalization,
% done by Seurat. If NOT, change 'true' to 'false'
save ('Stat_NormalizedExpression.mat', 'clusterID', 'pFisher', 'meanExp', ...
    'percenExp', 'oddRatio', 'pRanksum', 'logFoldChange'); %optional, save the gene stats in each cluster

% plot RYR2, a cardiomyocyte marker. high RYR2 expression (i.e. > 8)
% implies a cardiomyocyte cluster
markerExpression = 2.^( expressMat( find(ismember(gene, 'RYR2')==1) , : ) );
plotGeneExpressionOnUmap(markerExpression, umapCoor, 8);
% it suggests that clusters 1, 3, 7, 8, 9 are cardiomyocytes.

%plot the distribution of cells
allFigs = plotSampleCellDistribution (sampleID, umapCoor, clusterID);

%% ontology enrichment analysis analyze cell cycle G1S transition
% import the cell cycle G1S (S phase) marker list. The AI sparse analysis 
% only use these marker genes to differentiate between fetal (highly
% proliferative/G1S) and CTL-P56 (non-proliferative) cardiomyocytes.
G1SGeneList = importdata('G1S_Gene.txt');
G1SGeneIndex = find( ismember(gene, G1SGeneList) == 1 );

% get the expression matrix containing only cardiomyocyte and S phase
% marker
cmCluster = [1, 3, 7, 8, 9]; % we saw before, clusters 1, 2, 7, 8, 9 are cardiomyocytes
cmCellIndex = find( ismember(clusterID, cmCluster) == 1 );
G1S_cm_expression = expressMat (G1SGeneIndex, cmCellIndex);

% also get the cardiomyocyte cell group (sample ID)
cmCellGroup = sampleID(cmCellIndex);

% then perform sparse model analysis
[allCellScore, enrichTable, sparseModel, vioPlot] = ...
    sparseModelAnalyis(G1S_cm_expression, cmCellGroup, 'Fetal', 'CTL-P56');
save 'G1S_sparseModelResult.mat' allCellScore enrichTable sparseModel -mat;
% optional, save the G1S sparse analysis result.

%% utilities in the pipeline

