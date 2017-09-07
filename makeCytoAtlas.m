function makeCytoAtlas(atlasDir)
if nargin < 1, atlasDir = '/sni-storage/kalanit/users/zhenzl/AllenHumanBrainGeneExpression/data/Anatomy_Toolbox_22c/PMaps'; end


cwd = pwd;
cd(atlasDir);

%%  concate PMaps from anatomy toolbox
label = extractfield(dir('*.nii'),'name');
outvol = 'cytoPM.nii.gz';
invol = [];
fid = fopen('cytoLabel.txt','w+');
for i = 1:length(label)
    invol = [invol,' ',label{i}];
    
    [~,name,] = fileparts(label{i});
    fprintf(fid,'%d %s\n',i, strrep(name,'_','-'));
end
fclose(fid);

mergePM = sprintf('fslmerge -t %s %s',outvol,invol);
system(mergePM);

%% Pad zeros volume to 4D cytoPM
fslmerge -t cytoPMZeroPad.nii.gz AnatZerosMask.nii.gz cytoPM.nii.gz

%% Make cyto MPM
fslmaths cytoPMZeroPad.nii.gz -thr 0.25 -Tmaxn cytoMPM_thr25.nii.gz

%% swap orientation from RAS to LAS
fslswapdim cytoMPM_thr25.nii.gz -x y z  MNI152_cytoMPM_thr25.nii.gz
fslorient -swaporient MNI152_cytoMPM_thr25.nii.gz

%% resample MPM to FSL 
fslDir = getenv('FSLDIR');
refVol = fullfile(fslDir, '/data/standard/MNI152_T1_2mm.nii.gz');

flirt -in MNI152_cytoMPM_thr25.nii.gz -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz  -out MNI152_cytoMPM_thr25_2mm_new.nii.gz -usesqform -interp nearestneighbour



%% make yeo atlas 
% from LIA to LAS
fslswapdim Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz x z -y yeo7.nii.gz
% resample to 2mm
flirt -in yeo7_1mm.nii.gz -ref  $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz -out yeo7_2mm.nii.gz -applyxfm -usesqform  -interp nearestneighbour