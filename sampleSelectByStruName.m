function aha = sampleSelectByStruName(aha,struName)
% select sample based on brain structure
if nargin < 2, struName = 'CX'; end

switch struName
    case 'CX' % cerebral cortex
        cx = aha.sample.stru_id < 4250;        
        ins = aha.sample.stru_id >= 4270 & aha.sample.stru_id < 4280;
        tp =  aha.sample.stru_id >= 4800 & aha.sample.stru_id < 5000;
        index = cx | ins | tp;
        
    case 'THA' % Thalamus
        index =  aha.sample.stru_id >= 4400 & aha.sample.stru_id < 4540;
        
    case 'HYPO' % Hypothalamus
        index =  aha.sample.stru_id >= 4540 & aha.sample.stru_id < 4700;
        index = index | aha.sample.stru_id > 13000;
        
    case 'BG' %  Basal Forebrain and Basal Ganglia
        index = aha.sample.stru_id >= 4280 & aha.sample.stru_id < 4340;
        
    case 'AC' % Amygdaloid Complex
        index = aha.sample.stru_id >= 4340  &  aha.sample.stru_id < 4400;
        
    case 'HIP' % Hippocampus
        index = aha.sample.stru_id >= 4250 &  aha.sample.stru_id < 4270;
        
    case 'CB' % Cerebellum
        index = aha.sample.stru_id >= 4700 & aha.sample.stru_id < 4800;
        
    case 'MB' % Midbrain, Pons, and Medulla
        index =  aha.sample.stru_id >= 9000 & aha.sample.stru_id < 13000;
        cc = aha.sample.stru_id >= 9220 & aha.sample.stru_id < 9230;
        cg = aha.sample.stru_id >= 9240 & aha.sample.stru_id < 9250;
        wm  = cc | cg;
        index = index &(~wm);
        
    case 'WM' % White matter and Fiber
        cc = aha.sample.stru_id >= 9220 & aha.sample.stru_id < 9230;
        cg = aha.sample.stru_id >= 9240 & aha.sample.stru_id < 9250;
        index = cc | cg;
        
    otherwise
        error('Wrong structure name');
end


% display selected samples
stru_id = aha.sample.stru_id(index);
stru_name = aha.sample.stru_name(index);
[stru_id,I] = sort(stru_id);
stru_name = stru_name(I);
fprintf('(%s,%s): %d samples\n',aha.donor,struName,length(stru_name));
for i = 1:length(stru_id)
    fprintf('%d: %s\n',stru_id(i), stru_name{i});
end
% figure, histogram(aha.sample.stru_id,100 )


% update aha.samples
aha.sample.stru_id = aha.sample.stru_id(index);
aha.sample.slab_num = aha.sample.slab_num(index);
aha.sample.well_id=aha.sample. well_id(index);
aha.sample.slab_type = aha.sample.slab_type(index);
aha.sample.stru_acronym = aha.sample.stru_acronym(index);
aha.sample.stru_name = aha.sample.stru_name(index);
aha.sample.polygon_id = aha.sample.polygon_id(index);
aha.sample.nat_ijk = aha.sample.nat_ijk(index,:);
aha.sample.mni_xyz = aha.sample.mni_xyz(index,:);


% update aha.expression
aha.expression.value = aha.expression.value(:,index);


% update aha.call
aha.call.value = aha.call.value(:,index);



