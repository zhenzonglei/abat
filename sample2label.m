function roi = sample2SurfLabel(dataDir,sessid,roiId,geneFile,mniRoiFile,labelFile)
% roi = mapSample2SurfLabel(dataDir,sessid,roiId,geneFile,mniRoiFile)
% roiId is vector.
% read roi volume
% if nargin < 5, mniRoiFile = 'cytoMPM_thr25_2mm.nii.gz'; end
% if nargin < 4, geneFile = 'GeneFiletered.mat';

fid  = fopen(labelFile);
label = textscan(fid,'%d %s');
label = label{2};
fclose(fid);
 
    
ni = niftiRead(mniRoiFile);
roi = [];
for s = 1:length(sessid)
    donor = sessid{s};
    gene = load(fullfile(dataDir,donor,'gene',geneFile));
    
    % transform samples' mni coords to img coords
    sample_mni_coords = gene.sample.mni_xyz;
    sample_img_coords = round(mrAnatXformCoords(ni.qto_ijk, sample_mni_coords));    
    
    for i = 1:length(roiId)
        % get img coors for a roi
        [I,J,K] =  ind2sub(size(ni.data),find(ni.data == roiId(i)));
        roi_img_coords = [I,J,K];
        
        % find overlaping voxels
        [C,Is,~] = intersect(sample_img_coords,roi_img_coords,'rows');
        
        % get expression data
        % roi(s,i).donor = donor;
        roi(s,i).call = gene.call.value(Is,:);
        roi(s,i).label = label{i};
        roi(s,i).expr_val = gene.expression.value(Is,:);
        roi(s,i).sample_img_coords = C;
        roi(s,i).sample_mni_coords = sample_mni_coords(Is,:);
        roi(s,i).meas = [];
    end
end
