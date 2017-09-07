clear
close all
projDir = '/sni-storage/kalanit/users/zhenzl/AllenHumanBrainGeneExpression';
addpath(fullfile(projDir,'analysis','code'));
imgDir = fullfile(projDir,'image');

dataDir = fullfile(projDir, 'data');
sessid = {'H0351.1009', 'H0351.1012', 'H0351.1015', 'H0351.1016','H0351.2001','H0351.2002'};
nSubj = length(sessid);

% pcaFilterGene(dataDir, sessid)
 
%%  Get roi gene
mniRoiFile = fullfile(projDir,'data','cytoAtlas','cytoMPM_thr25_2mm.nii.gz');
labelFile = fullfile(projDir,'data','cytoAtlas','cytoLabel.txt');
roiId = 1:100;
geneFile ='GeneFiltered.mat';
roi = extractRoiGeneExpr(dataDir,sessid,roiId,geneFile,mniRoiFile,labelFile);
save(fullfile(projDir,'data','cytoAtlas','cytoGene.mat'),'roi');


%% plot sample number for each donor
sampleNum = zeros(100,6);
for s = 1:6
    for i = 1:100
        sampleNum(i,s) = size(roi(s,i).expr_val,2);
    end
end

label =  extractfield(roi(1,:), 'label');
figure('units','normalized','outerposition',[0 0 1 1],'name','Sample number');
imagesc(sampleNum)
axis square
colorbar
set(gca,'Ytick',1:100,'YtickLabel',label,'FontSize',8, 'TickDir','out','box','off');
set(gca,'Xtick',1:6,'XtickLabel', sessid);
sum(sampleNum>0)
sum(sum(sampleNum,2)>0)
fig = gcf;
fig.PaperPositionMode = 'auto';
imgName = fullfile(imgDir, 'donorSampleNumImg.png');
print('-dpng','-r0',imgName);



figure('units','normalized','outerposition',[0 0 1 1],'name','sampled number');
for s = 1:6
    subplot(2,3,s)
    label =  extractfield(roi(s,:), 'label');
    barh(sampleNum(:,s));
    axis square
    % set(gca,'Ytick',1:100,'YtickLabel',label,'FontSize',8, 'TickDir','out','box','off');
    % xlabel('Number of samples')
    set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
    axis off
    title(sessid{s})
end
fig = gcf;
fig.PaperPositionMode = 'auto';
imgName = fullfile(imgDir,'donorSampleNumBar.png');
print('-dpng','-r0',imgName);

%% plot sample number for pooled samples 
label =  extractfield(roi(1,:), 'label');
figure('units','normalized','outerposition',[0 0 1 1],'name','Pooled sampled number');
barh(sum(sampleNum,2));
axis square
set(gca,'Ytick',1:100,'YtickLabel',label,'FontSize',8, 'TickDir','out','box','off');
xlabel('Number of samples')
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
fig = gcf;
fig.PaperPositionMode = 'auto';
imgName = fullfile(imgDir, 'sampleNumPooled.png');
print('-dpng','-r0',imgName);


