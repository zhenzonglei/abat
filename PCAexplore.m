clear
close all
projDir = '/sni-storage/kalanit/users/zhenzl/AllenHumanBrainGeneExpression';
addpath(fullfile(projDir,'analysis','code'));
imgDir = fullfile(projDir,'image');

dataDir = fullfile(projDir, 'data');
sessid = {'H0351.1009', 'H0351.1012', 'H0351.1015', 'H0351.1016','H0351.2001','H0351.2002'};
nSubj = length(sessid);

figure('units','normalized','outerposition',[0 0 1 1],'name','PCA Variance');
for s = 1:nSubj
    gene = load(fullfile(dataDir,sessid{s},'gene','Gene.mat'));
    express_val = gene.expression.value;
    probe_id  = gene.expression.probe_id;
    
    [coeff, score, latent] = pca(zscore(express_val));
    latent = latent/sum(latent(:)) * 100;
    
    
    subplot(3,6,s),bar(latent(1:100))
    xlabel('PC');ylabel('Explained variance(%)');title( [sessid{s},':','With 1st PC']);
    axis square
    
    subplot(3,6,6+s),bar(latent(2:100));title('Without 1st PC');
    xlabel('PC');ylabel('Explained variance(%)')
    axis square
    
    subplot(3,6,12+s),plot(score(:,1:10));title('1st 10 scores');
    xlabel('Probe');ylabel('Score')
    axis square
end
fig = gcf;
fig.PaperPositionMode = 'auto';
imgName = fullfile(projDir,'image','pca.png');
print('-dpng','-r0',imgName);


close all
figure('units','normalized','outerposition',[0 0 1 1],'name','sample correlation');
for s = 1:nSubj
    gene = load(fullfile(dataDir,sessid{s},'gene','Gene.mat'));
    express_val = gene.expression.value;
    probe_id  = gene.expression.probe_id;
    for  i = 1:3
        I = randi(size(express_val,2),1,2);
        ax = subplot(3,6,(i-1)*6+s);
        ax = scatter(ax,express_val(:,I(1)), express_val(:,I(2)));
        lsline
        axis square
        title(sessid{s});
        xlabel(sprintf('Sample %d',I(1)));
        ylabel(sprintf('Sample %d',I(2)));
        
    end
end
fig = gcf;
fig.PaperPositionMode = 'auto';
imgName = fullfile(projDir,'image','sampleCorrScatter.png');
print('-dpng','-r0',imgName);




close all
figure('units','normalized','outerposition',[0 0 1 1],'name','sample distribution');
for s = 1:nSubj
    gene = load(fullfile(dataDir,sessid{s},'gene','Gene.mat'));
    express_val = gene.expression.value;
    probe_id  = gene.expression.probe_id;
    for  i = 1:3
        I = randi(size(express_val,2),1,1);
        ax = subplot(3,6,(i-1)*6+s);
        histogram(express_val(:,I));
        axis square
        title(sessid{s});
        xlabel(sprintf('Sample %d',I))
    end
end



fig = gcf;
fig.PaperPositionMode = 'auto';
imgName = fullfile(projDir,'image','sampleDistribution.png');
print('-dpng','-r0',imgName);













close all
figure('units','normalized','outerposition',[0 0 1 1],'name','Correlation Matrix');
for s = 1:nSubj
    gene = load(fullfile(dataDir,sessid{s},'gene','Gene.mat'));
    express_val = gene.expression.value;
    
    % correlation mat before filter
    ax1 = subplot(2,6,s);
    imagesc(corr(express_val));
    axis square
    title(sessid{s});
    colorbar('southoutside')
    
    % do pca
    [coeff, score] = pca(zscore(express_val));
    coeff(:,1) = 0;
    express_val = score*coeff';
    
    % correlation after filter
    ax2 = subplot(2,6,6+s);
    imagesc(corr(express_val));
    axis square
    title(sessid{s});
    colorbar('southoutside')
end

fig = gcf;
fig.PaperPositionMode = 'auto';
imgName = fullfile(projDir,'image','CorrMat.png');
print('-dpng','-r0',imgName);





close all
figure('units','normalized','outerposition',[0 0 1 1],'name','clustering');
for s = 1:nSubj
    gene = load(fullfile(dataDir,sessid{s},'gene','Gene.mat'));
    express_val = gene.expression.value;
    
    % cluster before filter
    ax1 = subplot(2,6,s);
    Y = pdist(express_val','correlation');
    Z = linkage(Y,'average');
    [H, T] = dendrogram(Z);
    title(sessid{s});
    set(gca,'xtick',[])
    
    
    % do pca
    [coeff, score] = pca(zscore(express_val));
    coeff(:,1) = 0;
    express_val = score*coeff';
    
    % cluster after filter
    ax2 = subplot(2,6,6+s);
    Y = pdist(express_val','correlation');
    Z = linkage(Y,'average');
    [H, T] = dendrogram(Z);
    title(sessid{s});
    set(gca,'xtick',[])
end

fig = gcf;
fig.PaperPositionMode = 'auto';
imgName = fullfile(projDir,'image','clustering.png');
print('-dpng','-r0',imgName);


figure, imagesc(zscore(express_val(1:100,:)));












% makeCytoAtlas()
% mniRoiFile = fullfile(projDir,'data','cytoAtlas','cytoMPM_thr25_2mm.nii.gz');
% roiId = 85;
% expr_val = extractRoiGeneExpr(dataDir,sessid(1),mniRoiFile,roiId);
% figure, plot(expr_val)
