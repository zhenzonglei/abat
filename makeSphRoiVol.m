function vol = makeSphRoiVol(mask,centerCoord, radius)
% coords = dtiBuildSphereCoords(centerCoord, radius)
%
% centerCoord: 1x3 coordinate defining the sphere center.
% radius: scalar defining the radius, in voxel units
%
% If centerCoord is a 1x2, then we do a 2d version.
%
% RETURNS:
%  coords: an Nx3 list of the coordinates that intersect the sphere.
%  (Or Nx2 if the 2d option is selected by passing in a 1x2.)
%
% HISTORY:
% 2003.12.01 RFD (bob@white.stanford.edu) wrote it.
% 2006.09.12 RFD: added 2d option.

if(radius>0)
    if(size(centerCoord,2)==2)
        [X,Y] = meshgrid([-radius:+radius],...
            [-radius:+radius]);
        dSq = X.^2+Y.^2;
        keep = dSq(:) < radius.^2;
        coords = [X(keep)+centerCoord(1), Y(keep)+centerCoord(2)];
    else
        [X,Y,Z] = meshgrid([-radius:+radius],...
            [-radius:+radius],...
            [-radius:+radius]);
        dSq = X.^2+Y.^2+Z.^2;
        keep = dSq(:) < radius.^2;
        coords = [X(keep)+centerCoord(1), Y(keep)+centerCoord(2), Z(keep)+centerCoord(3)];
    end
else
    coords = centerCoord;
end







% neighbor list
offsets = mk_neighbor('sphere',radius);
nvox_in_full_sphere = size(offsets,1);


mask_dims = size(mask);
mask_coords = [];
[mask_coords(:,1) mask_coords(:,2) mask_coords(:,3)] = ind2sub(mask_dims,find(mask));
nvox_active_in_mask = size(mask_coords,1);

slinx = [];
for v = 1:nvox_active_in_mask
    vox_coords_tiled = repmat(mask_coords(v,:), nvox_in_full_sphere, 1);
    sphere_coords = vox_coords_tiled + offsets;

    % winnow voxels that lie outside the boundaries of the volume
    badvoxels = find( sphere_coords(:,1)<1 | ...
        sphere_coords(:,2)<1 | ...
        sphere_coords(:,3)<1 | ...
        sphere_coords(:,1)>mask_dims(1) | ...
        sphere_coords(:,2)>mask_dims(2) | ...
        sphere_coords(:,3)>mask_dims(3));
    sphere_coords(badvoxels,:) = [];

    sphere_idx_mask = sub2ind(mask_dims, ...
        sphere_coords(:,1), ...
        sphere_coords(:,2), ...
        sphere_coords(:,3));

    % winnow voxels that aren't in the MASK.if voxel in mask,
    % element of sphere_idx_included_in_mask = 1
    sphere_idx_included_in_mask = mask(sphere_idx_mask);

    % how many voxels in this sphere, after throwing away
    % out-of-bounds voxels and voxels excluded by the mask.
    idx_pat = sphere_idx_mask(logical(sphere_idx_included_in_mask));

    % index for this searchlight
    slinx{v} = idx_pat;
end














return;