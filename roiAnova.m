function sig = roiAnova(roi)
% do anova to test difference between rois

nRoi = length(roi);
nGene = size(roi(1).expr_val,2);


sample_label = {};
nSample = zeros(nRoi,1);
for i = 1:nRoi
    nSample(i) = size(roi(i).sample_img_coords,1);
    L = cell(1,nSample(i));
    L(:) =  {roi(i).label};
    sample_label = [sample_label, L];
end


sig = nan(nGene,1);
GE = nan(sum(nSample),1);
for g = 1:nGene
    s = 1;
    for i = 1:nRoi
        t = s + nSample(i)-1;
        GE(s:t) = roi(i).expr_val(:,g);
        s = t;
    end
    sig(g) = -log10(anova1(GE, sample_label,'off'));
end
