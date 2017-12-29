function subj  = extractAnnotationBrainMeasNew(datadir,sessid,run,annotationfile,hemi,mode)
% roi = extractAnnotationBrainMeas(dataDir,sessid,run,annotationFile,hemi,mode)

if nargin < 6, mode = 'mni'; end
if nargin < 5, hemi = 'lh'; end

if strcmp(mode,'mni')
    mni =  true;
elseif strcmp(mode,'native')
    mni = false;
else
    error('Wrong mode, should be mni or native');
end


subjdir = getenv('SUBJECTS_DIR');
for s = 1:length(sessid)
    if mni
        annotation = fullfile(subjdir,'fsaverage5','label',[hemi,'.',annotationfile]);
    else
        annotation = fullfile(subjdir, sessid{s},'label',[hemi,'.',annotationfile]);
    end
    
    % read annotation file
    [~,label,ct] = read_annotation(annotation);
    
    % read surface for geometry ref
    L = unique(label); L(L==0) =[];
    nlabel = length(L);
    
    subj(s).name = sessid{s};
    subj(s).hemi = hemi;
    for i = 1:nlabel
        subj(s).label(i,1) = ct.struct_names(ct.table(:,5) == L(i));
    end
    
    
    roi =[];
    for r = 1:length(run)
        mrFile = fullfile(datadir, sessid{s},'bold',['gump.fs5.fc.', hemi], ['pr',run{r}], ...
            'res',['res-',run{r},'.nii.gz']);
        
        mri = MRIread(mrFile);
        meas =[];
        for i = 1:nlabel
            idx = label==L(i);
            meas(i,:) = squeeze(mean(mri.vol(:,idx,:),2));
        end
        roi{r} = meas;
    end
    subj(s).roi = roi;
end