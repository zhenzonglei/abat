function pr =  poolRoi(roi,group,group_name)
% pool roi according to the group
% roi, nDonor x nRoi structure array
% pr, nDonorxnGroup structure array
if nargin < 3, group_name = {'VP','DP'}; end
if nargin < 2, group={[81:84,88,92], [87,89,90,91,93]}; end
nDonor = size(roi,1);
nGroup = length(group);

pr(nDonor,nGroup).label =[];
pr(nDonor,nGroup).expr_val = [];
pr(nDonor,nGroup).sample_img_coords = [];
pr(nDonor,nGroup).sample_mni_coords = [];
pr(nDonor,nGroup).group = [];


if isfield(roi(1,1),'meas')
    hasMeas = true;
    pr(nDonor,nGroup).meas = [];
else
    hasMeas = false;
end



for d = 1:nDonor
    for g = 1:nGroup
        G = group{g};
        pr(d,g).label= group_name{g};
        for i = 1:length(G) 
            pr(d,g).expr_val = [pr(d,g).expr_val; roi(d,G(i)).expr_val];
            pr(d,g).sample_img_coords = [pr(d,g).sample_img_coords; roi(d,G(i)).sample_img_coords];
            pr(d,g).sample_mni_coords = [pr(d,g).sample_mni_coords; roi(d,G(i)).sample_mni_coords];
            
            if hasMeas
                pr(d,g).meas =  [pr(d,g).meas; roi(d,G(i)).meas]; 
            end
        end
    end
end