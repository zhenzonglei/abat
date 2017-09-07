function aha = geneSelectByCall(aha,callRatio)
% select gene based on PACall
if nargin < 2, callRatio = 0.05; end

call = aha.call.value;
index = mean(call,2) > callRatio; 
fprintf('%s, callRatio= %.2f, %.2f%% genes left\n',aha.donor, callRatio,100*sum(index)/length(index));


 
%% update aha.expression
aha.expression.probe_id = aha.expression.probe_id(index,:);
aha.expression.value = aha.expression.value(index,:);


%% update aha.call
aha.call.probe_id = aha.call.probe_id(index,:);
aha.call.value = aha.call.value(index,:);

%% update aha.probe
aha.probe.id = aha.probe.id(index);
aha.probe.name = aha.probe.name(index);
aha.probe.gene_id = aha.probe.gene_id(index);
aha.probe.gene_symbol = aha.probe.gene_symbol(index);
aha.probe.gene_name = aha.probe.gene_name(index);
aha.probe.entrez_id = aha.probe.entrez_id(index);
aha.probe.chromosome = aha.probe.chromosome(index);

