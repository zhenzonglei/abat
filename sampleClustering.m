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
rid = [81:93];% ,100]; % roi id
nRoi = length(rid);


for s = 1:nSubj
    figure('units','normalized','outerposition',[0 0 1 1],'name','Roi clustering');
    genePat = zeros(nRoi,nProb);
    label =  extractfield(roi(s,rid), 'label');
    
    for r = 1:nRoi
        genePat(r,:) = mean(roi(s,rid(r)).expr_val,2)';
    end
    idx = ~all(isnan(genePat),2);
    genePat = genePat(idx,:);
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
    
%     fig = gcf;
%     fig.PaperPositionMode = 'auto';
%     imgName = fullfile(projDir,'image',[sessid{s},'_visual_roiClustering.png']);
%     print('-dpng','-r0',imgName);
end


%% top gene
close all
toppct = 1;
figure('units','normalized','outerposition',[0 0 1 1],'name',sprintf('Roi clustering-top %d percent',toppct));
for s = 1:nSubj
    genePat = zeros(nRoi,nProb);
    label =  extractfield(roi(s,rid), 'label');
    
    for r = 1:nRoi
        genePat(r,:) = mean(roi(s,rid(r)).expr_val,2)';
    end
    idx = ~all(isnan(genePat),2);
    genePat = genePat(idx,:);
    label = label(idx);
    
    for i = 1:length(label)
        k = strfind(label{i},'-');
        label{i} = label{i}(k+1:end);
    end
   % label{end} = 'wThalamus';
    
    
    geneSd = std(genePat);
    % geneMean = mean(genePat);
    % geneCv = geneSd./geneMean;
    I = geneSd > prctile(geneSd,100-toppct);
    %     subplot(3,6,s); bar(geneCv(I))
    %     subplot(3,6,6+s); bar(geneMean(I))
    %     subplot(3,6,12+s); bar(geneSd(I))
    
    genePat = genePat(:,I);
  
    load(fullfile(dataDir, sessid{s},'gene','Gene.mat'),'probe');
    gene_symbol = probe.gene_symbol(I);
    
    
    fid = fopen(fullfile(dataDir, sessid{s},'gene',sprintf('top%dGene_sd.txt',toppct)),'w+');
    [Ycv,Icv] = sort(geneSd(I),2,'descend');
    symbol = gene_symbol(Icv);
    for i = 1:length(symbol)
        fprintf(fid,'%s, %.2f\n',symbol{i},Ycv(i));
    end
    fclose(fid);
    
    ax1 = subplot(2,6,s);
    imagesc(corr(genePat'));
    axis square
    title(sessid{s});
    colorbar('southoutside')
    set(gca,'Ytick',1:length(label),'YtickLabel',label,'FontSize',8, 'TickDir','out','box','off');
    
    ax2 = subplot(2,6,6+s);
    Y = pdist(genePat,'correlation');
    Z = linkage(Y,'average');
    [H, T] = dendrogram(Z,0, 'Orientation', 'right','Labels',label);
    title(sessid{s});
end

fig = gcf;
fig.PaperPositionMode = 'auto';
imgName = fullfile(projDir,'image',sprintf('visual_top%dgene_roiClustering.png',toppct));
print('-dpng','-r0',imgName);


%% bottom gene
botpct = 1;
figure('units','normalized','outerposition',[0 0 1 1],'name',sprintf('Roi clustering-bot %d percent',botpct));
for s = 1:nSubj
    genePat = zeros(nRoi,nProb);
    label =  extractfield(roi(s,rid), 'label');
    
    for r = 1:nRoi
        genePat(r,:) = mean(roi(s,rid(r)).expr_val,2)';
    end
    idx = ~all(isnan(genePat),2);
    genePat = genePat(idx,:);
    label = label(idx);
    
    for i = 1:length(label)
        k = strfind(label{i},'-');
        label{i} = label{i}(k+1:end);
    end
%     label{end} = 'wThalamus';
    
    geneStd = std(genePat);
    I = geneStd < prctile(geneStd,botpct);
    genePat = genePat(:,I);
    
    ax1 = subplot(2,6,s);
    imagesc(corr(genePat'));
    axis square
    title(sessid{s});
    colorbar('southoutside')
    set(gca,'Ytick',1:length(label),'YtickLabel',label,'FontSize',8, 'TickDir','out','box','off');
    
    ax2 = subplot(2,6,6+s);
    Y = pdist(genePat,'correlation');
    Z = linkage(Y,'average');
    [H, T] = dendrogram(Z,0, 'Orientation', 'right','Labels',label);
    title(sessid{s});    
end
fig = gcf;
fig.PaperPositionMode = 'auto';
imgName = fullfile(projDir,'image',sprintf('visual_bot%dgene_roiClustering.png',botpct));
print('-dpng','-r0',imgName);


%% random gene
randpct = 1;
figure('units','normalized','outerposition',[0 0 1 1],'name',sprintf('Roi clustering-rand %d percent',randpct));
for s = 1:nSubj
    genePat = zeros(nRoi,nProb);
    label =  extractfield(roi(s,rid), 'label');
    
    for r = 1:nRoi
        genePat(r,:) = mean(roi(s,rid(r)).expr_val,2)';
    end
    idx = ~all(isnan(genePat),2);
    genePat = genePat(idx,:);
    label = label(idx);
    
    for i = 1:length(label)
        k = strfind(label{i},'-');
        label{i} = label{i}(k+1:end);
    end
%     label{end} = 'wThalamus';
    I = randi(nProb,round(nProb*randpct/100),1);
    genePat = genePat(:,I);
    
    ax1 = subplot(2,6,s);
    imagesc(corr(genePat'));
    axis square
    title(sessid{s});
    colorbar('southoutside')
    set(gca,'Ytick',1:length(label),'YtickLabel',label,'FontSize',8, 'TickDir','out','box','off');
    
    ax2 = subplot(2,6,6+s);
    Y = pdist(genePat,'correlation');
    Z = linkage(Y,'average');
    [H, T] = dendrogram(Z,0, 'Orientation', 'right','Labels',label);
    title(sessid{s});
    
%     load(fullfile(dataDir, sessid{s},'gene','Gene.mat'),'probe');
%     %     length(unique(probe.gene_symbol(I)))
%     A = probe.gene_symbol(I);
%     fid = fopen(fullfile(dataDir, sessid{s},'gene',sprintf('top%dGene.txt',botpct)),'w+');
%     for i = 1:length(A)
%         fprintf(fid,'%s\n',A{i});
%     end
%     fclose(fid);
    
    %     figure
    %     [idx,label] = grp2idx(sort(A));
    %     histogram(idx,unique(idx));
    %     N = histcounts(idx,unique(idx));
    %     label(N>1)
    %     set(gca,'xTickLabel',label)
    
end
fig = gcf;
fig.PaperPositionMode = 'auto';
imgName = fullfile(projDir,'image',sprintf('visual_rand%dgene_roiClustering.png',randpct));
print('-dpng','-r0',imgName);





% %% top gene frequency
% close all
% toppct = 1;
% figure('units','normalized','outerposition',[0 0 1 1],'name',sprintf('Roi clustering-top %d percent',toppct));
% for s = 1:nSubj
%     genePat = zeros(nRoi,nProb);
%     label =  extractfield(roi(s,rid), 'label');
%     
%     for r = 1:nRoi
%         genePat(r,:) = mean(roi(s,rid(r)).expr_val,2)';
%     end
%     idx = ~all(isnan(genePat),2);
%     genePat = genePat(idx,:);
%     label = label(idx);
%     
%     for i = 1:length(label)
%         k = strfind(label{i},'-');
%         label{i} = label{i}(k+1:end);
%     end
%     label{end} = 'wThalamus';
%     
%     
%     geneSd = std(genePat);
%     % geneMean = mean(genePat);
%     % geneCv = geneSd./geneMean;
%     I = geneSd > prctile(geneSd,100-toppct);
%     %     subplot(3,6,s); bar(geneCv(I))
%     %     subplot(3,6,6+s); bar(geneMean(I))
%     %     subplot(3,6,12+s); bar(geneSd(I))
%     
%     subplot(2,6,s); bar(geneSd(I))
%     genePat = genePat(:,I);
%     
%     
%     load(fullfile(dataDir, sessid{s},'gene','Gene.mat'),'probe');
%     gene_symbol = probe.gene_symbol(I);
%     
%     
%     fid = fopen(fullfile(dataDir, sessid{s},'gene',sprintf('top%dGene_sd.txt',toppct)),'w+');
%     [Ycv,Icv] = sort(geneSd(I),2,'descend');
%     symbol = gene_symbol(Icv);
%     for i = 1:length(symbol)
%         fprintf(fid,'%s, %.2f\n',symbol{i},Ycv(i));
%     end
%     fclose(fid);
%     %
%     
%     
%     
%     %     ax1 = subplot(2,6,s);
%     %     imagesc(corr(genePat'));
%     %     axis square
%     %     title(sessid{s});
%     %     colorbar('southoutside')
%     %     set(gca,'Ytick',1:length(label),'YtickLabel',label,'FontSize',8, 'TickDir','out','box','off');
%     %
%     %     ax2 = subplot(2,6,6+s);
%     %     Y = pdist(genePat,'correlation');
%     %     Z = linkage(Y,'average');
%     %     [H, T] = dendrogram(Z,0, 'Orientation', 'right','Labels',label);
%     %     title(sessid{s});
%     
%     
%     
%     [idx,symbol] = grp2idx(sort(gene_symbol));
%     subplot(2,6,6+s); histogram(idx,unique(idx));
%     N = histcounts(idx,length(unique(idx)))';
%     
%     [Yf,If]= sort(N,1,'descend');
%     fid = fopen(fullfile(dataDir, sessid{s},'gene',sprintf('top%dGene_freq.txt',toppct)),'w+');
%     symbol = symbol(If);
%     for i = 1:length(symbol)
%         fprintf(fid,'%s %d\n',symbol{i},Yf(i));
%     end
%     fclose(fid);
% % xx = genePat(:,If)
%     
% end
% 
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% imgName = fullfile(projDir,'image',sprintf('visual_top%dgene_Count.png',toppct));
% print('-dpng','-r0',imgName);
% 
% 
% 
% figure, bar(genePat(:,4))
