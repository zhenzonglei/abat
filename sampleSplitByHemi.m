function split_aha = sampleSplitByHemi(aha)
% split sample by hemisphere

hemi = {'left','right'};
stru_name = aha.sample.stru_name;

for h = 1:length(hemi)
    index = not(cellfun('isempty',strfind(stru_name, hemi{h})));
    
    % update aha.samples
    split_aha(h).sample.stru_id = aha.sample.stru_id(index);
    split_aha(h).sample.slab_num = aha.sample.slab_num(index);
    split_aha(h).sample.well_id=aha.sample. well_id(index);
    split_aha(h).sample.slab_type = aha.sample.slab_type(index);
    split_aha(h).sample.stru_acronym = aha.sample.stru_acronym(index);
    split_aha(h).sample.stru_name = aha.sample.stru_name(index);
    split_aha(h).sample.polygon_id = aha.sample.polygon_id(index);
    split_aha(h).sample.nat_ijk = aha.sample.nat_ijk(index,:);
    split_aha(h).sample.mni_xyz = aha.sample.mni_xyz(index,:);
    
    
    % update aha.expression
    split_aha(h).expression.value = aha.expression.value(:,index);
    
    
    % update aha.call
    split_aha(h).call.value = aha.call.value(:,index);
end



