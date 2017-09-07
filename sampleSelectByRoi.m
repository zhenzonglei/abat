function roiaha= sampleSelectByRoi(aha,roiNi,roiId,roiName)

% transform samples' mni coords to img coords
sample_ijk = round(mrAnatXformCoords(roiNi.qto_ijk,aha.sample.mni_xyz));
for i = 1:length(roiId)
    fprintf('Select samples for (%s, %s)\n',aha.donor,roiName{i});
    % get img coors for a roi
    [I,J,K] =  ind2sub(size(roiNi.data),find(roiNi.data == roiId(i)));
    roi_ijk = [I,J,K];
    
    % find overlaping voxels
    [~,index,~] = intersect(sample_ijk,roi_ijk,'rows');
    
    % update aha.name
    roiaha(i).name = roiName{i};
    
    % update aha.samples
    roiaha(i).sample.stru_id = aha.sample.stru_id(index);
    roiaha(i).sample.slab_num = aha.sample.slab_num(index);
    roiaha(i).sample.well_id=aha.sample. well_id(index);
    roiaha(i).sample.slab_type = aha.sample.slab_type(index);
    roiaha(i).sample.stru_acronym = aha.sample.stru_acronym(index);
    roiaha(i).sample.stru_name = aha.sample.stru_name(index);
    roiaha(i).sample.polygon_id = aha.sample.polygon_id(index);
    roiaha(i).sample.nat_ijk = aha.sample.nat_ijk(index,:);
    roiaha(i).sample.mni_xyz = aha.sample.mni_xyz(index,:);
    
    % update aha.expression
    roiaha(i).expression.value = aha.expression.value(:,index);
    
    % update aha.call
    roiaha(i).call.value = aha.call.value(:,index);
end

