function roi = sample2SurfAnnotation(dataDir,sessid,geneFile,annotationFile,surfFile,dist_thr,hemi,mode)
% roi = sample2SurfAnnotation(dataDir,sessid,geneFile,annotationFile,mode)
% geneFile, allen gene data file.
% annotationFile, annotation file name
% surfFile, surface file
% mode, mni or native

if nargin < 8, mode = 'mni'; end
if nargin < 7, hemi = 'lh'; end
if nargin < 6, dist_thr = 4; end
if nargin < 5, surfFile = 'pial'; end

if strcmp(mode,'mni')
    mni =  true; 
elseif strcmp(mode,'native')
    mni = false;
else
    error('Wrong mode, should be mni or native');
end

subjdir = getenv('SUBJECTS_DIR');
for s = 1:length(sessid)
    donor = sessid{s};
    if mni
        annotation = fullfile(subjdir,'fsaverage','label',[hemi,'.',annotationFile]);
        surf = fullfile(subjdir, 'fsaverage','surf',[hemi,'.',surfFile]);
    else 
        annotation = fullfile(subjdir, donor,'label',[hemi,'.',annotationFile]);
        surf = fullfile(subjdir, donor,'surf',[hemi,'.',annotationFile]);
    end
    
    % read annotation file
    [~,label,ct] = read_annotation(annotation);
    [label_coords,~] = read_surf(surf);
    
    % read surface for geometry ref
    L = unique(label); L(L==0) =[];
    nlabel = length(L);
    
    % read allen gene file
    gene = load(fullfile(dataDir,donor,'gene',geneFile));
    sample_coords = gene.sample.mni_xyz;
    if strcmp(hemi,'lh')
        x = sample_coords(:,1) < 0 ;
    else
        x = sample_coords(:,1) > 0 ;
    end
    
    
    for i = 1:nlabel
        D = pdist2(sample_coords,label_coords(label == L(i),:)); 
        Is = any( D < dist_thr,2) & x;
              
        %% get expression data
        roi(s,i).donor = donor;
        roi(s,i).hemi = hemi;
        roi(s,i).call = gene.call.value(Is,:);
        roi(s,i).label_name = ct.struct_names(ct.table(:,5) == L(i));
        % roi(s,i).label_id = L(i);
        roi(s,i).expr_val = gene.expression.value(Is,:);
        roi(s,i).sample_coords = sample_coords(Is,:);
        roi(s,i).meas = [];
    end
end