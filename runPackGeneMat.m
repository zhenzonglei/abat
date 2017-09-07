clear
close all
projDir = '/sni-storage/kalanit/users/zhenzl/AllenHumanBrainGeneExpression';
addpath(fullfile(projDir,'analysis','code'));

dataDir = fullfile(projDir, 'data');
sessid = {'H0351.1009', 'H0351.1012', 'H0351.1015', 'H0351.1016','H0351.2001','H0351.2002'};
nSubj = length(sessid);
for s = 1:length(sessid)
    packGeneData2Mat(dataDir, sessid(s))
end
