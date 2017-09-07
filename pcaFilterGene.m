function pcaFilterGene(dataDir, sessid)
%  pcaFilterGene(dataDir,sessid)
for s = 1:length(sessid)
    % load data
    donor = sessid{s};
    geneFile = fullfile(dataDir,donor,'gene','Gene.mat');
    gene = load(geneFile);      
    express_val = gene.expression.value;
   
    % do filtering
    fprintf('Filter gene data for %s(probe,sample):(%d,%d)\n', donor,size(express_val));   
    [coeff, score] = pca(zscore(express_val));
    coeff(:,1) = 0;
    express_val = score*coeff'; 
    gene.expression.value = express_val; 
    
    % save filtered data   
    outFile = fullfile(dataDir,donor,'gene','GeneFiltered.mat');
    save(outFile,'-struct', 'gene');
end