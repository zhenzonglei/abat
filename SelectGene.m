% feature selection and cluster
clear
close all
projDir = '/sni-storage/kalanit/users/zhenzl/AllenHumanBrainGeneExpression';
addpath(fullfile(projDir,'analysis','code'));
imgDir = fullfile(projDir,'image');

dataDir = fullfile(projDir, 'data');
sessid = {'H0351.1009', 'H0351.1012', 'H0351.1015', 'H0351.1016','H0351.2001','H0351.2002'};
nSubj = length(sessid);


%%  Get roi gene
load(fullfile(projDir,'data','cytoAtlas','cytoGene.mat'),'roi');
nProb = 58692;
nRoi = 100;

for s = 1:nSubj
    figure('units','normalized','outerposition',[0 0 1 1],'name','Roi clustering');
    genePat = zeros(nRoi,nProb);
    for r = 1:nRoi
        genePat(r,:) = mean(roi(s,r).expr_val,2)';
    end
    idx = ~all(isnan(genePat),2);
    genePat = genePat(idx,:);
    label =  extractfield(roi(s,:), 'label');
    label = label(idx);
    
    
    ax1 = subplot(1,2,1);
    imagesc(corr(genePat'));
    axis square
    title(sessid{s});
    colorbar('southoutside')
    set(gca,'Ytick',1:length(label),'YtickLabel',label,'FontSize',8, 'TickDir','out','box','off');
    
    ax2 = subplot(1,2,2);
    Y = pdist(genePat,'correlation');
    Z = linkage(Y,'average');
    [H, T] = dendrogram(Z,0, 'Orientation', 'right','Labels',label);
    title(sessid{s});
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    imgName = fullfile(projDir,'image',[sessid{s},'_roiClustering.png']);
    print('-dpng','-r0',imgName);
end




