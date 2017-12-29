function extractEntrezGene(dataDir,sessid,raw_file)
% Extract the probes which are in entrez
% dataDir, Dir for session
% sessid, donor Id
% raw_file, raw ABA data file

if nargin < 3, raw_file='raw.mat'; end

% only the genes with entrez id and correct symbol are kept.
gene = load(fullfile(dataDir,sessid{1},'gene',raw_file));
entrez = gene.probe.entrez_id~=0;
na = strcmp(gene.probe.gene_symbol,'na');
entrez(na) = false;
clear gene;

for s = 1:length(sessid)
    donor = sessid{s};
    fprintf('Extract Entrez done for %s(sample,probe):', donor);
    load(fullfile(dataDir,donor,'gene',raw_file));
    
    % call
    call.probe_id = call.probe_id(entrez);
    call.value = call.value(:,entrez);

    % expre
    expression.probe_id = expression.probe_id(entrez);
    expression.value = expression.value(:,entrez);
   
    % probe
    probe.id = probe.id(entrez);
    probe.name = probe.name(entrez);
    probe.gene_id = probe.gene_id(entrez);
    probe.gene_symbol = probe.gene_symbol(entrez);
    probe.gene_name = probe.gene_name(entrez);
    probe.entrez_id = probe.entrez_id(entrez);
    probe.chromosome = probe.chromosome(entrez);
    
    % save raw entrez file
    outFile = fullfile(dataDir,sessid{s},'gene','raw_entrez.mat');
    save(outFile,'donor','sample', 'call','expression', 'probe')
    fprintf('(%d,%d)\n', size(expression.value));
end