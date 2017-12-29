function roi =  extractRoiBrainMeas(roi,mniFile,radius)
% roi =  extractRoiBrainMeas(roi,mniFile,radius)
if nargin < 3, radius = 3;end
sph = makeSphRoiCoords([0 0 0 ],radius);
[nDonor,nRoi]= size(roi);
mni = niftiRead(mniFile);
meas = reshape(mni.data,[],size(mni.data,4));
for s = 1:nDonor
    for i =  1:nRoi
        nSample = size(roi(s,i).sample_img_coords,1);
        roi(s,i).meas = zeros(nSample, size(mni.data,4));
        for c = 1:nSample
            sample_img_coords =  roi(s,i).sample_img_coords(c,:);
            
            sphere_coords = [sample_img_coords(1)+sph(:,1), ...
                sample_img_coords(2) + sph(:,2),...
                sample_img_coords(3)+sph(:,3)];
            
            % remove bad coords that lie outside the boundaries of the volume
            bad = sphere_coords(:,1)<1 | ...
                sphere_coords(:,2)<1 | ...
                sphere_coords(:,3)<1 | ...
                sphere_coords(:,1)>mni.dim(1) | ...
                sphere_coords(:,2)>mni.dim(2) | ...
                sphere_coords(:,3)>mni.dim(3);
            sphere_coords(bad,:) = [];
            
            
            I = sub2ind(mni.dim(1:3), sphere_coords(:,1), sphere_coords(:,2), sphere_coords(:,3));
            roi(s,i).meas(c,:) = mean(meas(I,:),1);
        end
    end
end