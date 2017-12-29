function makeSampleProbMap(dataDir,sessid,mniFile,radius)
% makeSampleProbMap(dataDir,sessid,mniFile,radius)
if nargin < 4, radius = 5;end

 sph = makeSphRoiCoords([0 0 0],radius);





mni = niftiRead(mniFile);
for s = 1:length(sessid)
    donor = sessid{s};
    annotFile = fullfile(dataDir,donor,'gene','SampleAnnot.csv');
    sample = readSampleAnnot(annotFile);
   
    % transform samples' mni coords to img coords
    sample_mni_coords = sample.mni_xyz;
    sample_img_coords = round(mrAnatXformCoords(mni.qto_ijk, sample_mni_coords));
    
    
    pm = zeros(mni.dim);
    for i = 1:size(sample_mni_coords,1)
        sphere_coords = [sample_mni_coords(i,1)+sph(:,1), ...
            sample_mni_coords(i,2) + sph(:,2),...
            sample_mni_coords(i,3)+sph(:,3)];
                        
        % remove bad coords that lie outside the boundaries of the volume
        bad = sphere_coords(:,1)<1 | ...
            sphere_coords(:,2)<1 | ...
            sphere_coords(:,3)<1 | ...
            sphere_coords(:,1)>mni.dim(1) | ...
            sphere_coords(:,2)>mni.dim(2) | ...
            sphere_coords(:,3)>mni.dim(3);
        sphere_coords(bad,:) = [];
      
        
        % make a 3D image
        vol = zeros(mni.dim);
        vol(sub2ind(mni.dim, sphere_coords(:,1), sphere_coords(:,2), sphere_coords(:,3))) = 1;
        
        % add up
        pm = pm + vol;    
    end
    
    
    % change img data
    mni.data = pm;
    mni.cal_min = min(pm(:));
    mni.cal_max = max(pm(:));
    
    % write the nifti file
    mni.fname = fullfile(dataDir,donor,'gene', 'sampleProbMap.nii.gz');
    niftiWrite(mni);
end
