function pd = poolDonor(roi,group)
% Pool donor into N group for each roi
% nDonor x nRoi structure array
% group, vector indicate the group id for each donor: 1,2,..N
% pd, ngroupxnRoi structure array

[nDonor,nRoi] = size(roi);
if nargin < 2, group = ones(1,nDonor); end
ngroup = max(group);

pd(ngroup,nRoi).call = [];
pd(ngroup,nRoi).label = [];
pd(ngroup,nRoi).expr_val = [];
pd(ngroup,nRoi).sample_img_coords = [];
pd(ngroup,nRoi).sample_mni_coords = [];

if isfield(roi(1,1),'meas')
    hasMeas = true;
    pd(ngroup,nRoi).meas = [];
else
    hasMeas = false;
end

for g = 1:ngroup
    donor = find(group==g);
    for r = 1:nRoi
        for d = donor
            pd(g,r).expr_val = [pd(g,r).expr_val; roi(d,r).expr_val];
            pd(g,r).sample_img_coords = [pd(g,r).sample_img_coords; roi(d,r).sample_img_coords];
            pd(g,r).sample_mni_coords = [pd(g,r).sample_mni_coords; roi(d,r).sample_mni_coords];
            if hasMeas
                pd(g,r).meas = [pd(g,r).meas; roi(d,r).meas];
            end
        end
    end
end