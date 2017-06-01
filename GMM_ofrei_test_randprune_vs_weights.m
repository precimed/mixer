t = 'H:\NORSTORE\MMIL';                       if exist(t, 'dir'), mmil_gwas_path = t; end;
t = '/work/users/oleksanf/NORSTORE/MMIL';     if exist(t, 'dir'), mmil_gwas_path = t; end;
t = '/space/syn03/1/data/GWAS';               if exist(t, 'dir'), mmil_gwas_path = t; end;
if ~exist('mmil_gwas_path', 'var'), error('Unable to locate data folder'); end;

use_2m_template = true;
nprune = 10000;

if use_2m_template,
    if ~exist('LDr2_p1sparse', 'var'), load(fullfile(mmil_gwas_path, 'old_25_SNPs_processed_summary/GWAS_Annot/ldmat_2m_p1.mat')); end;
    if ~exist('LDr2_p2sparse', 'var'), load(fullfile(mmil_gwas_path, 'old_25_SNPs_processed_summary/GWAS_Annot/ldmat_2m_p2.mat')); end;
    if ~exist('LDr2_p8sparse', 'var'), load(fullfile(mmil_gwas_path, 'old_25_SNPs_processed_summary/GWAS_Annot/ldmat_2m_p8.mat')); end;
    result_file = 'w_ld_2m.mat';
else
    if ~exist('LDr2_p1sparse', 'var'), load(fullfile(mmil_gwas_path, 'new_9m_SNPs_processed_summary/GWAS_Annot/ldmat_9m_p1.mat')); end;
    if ~exist('LDr2_p2sparse', 'var'), load(fullfile(mmil_gwas_path, 'new_9m_SNPs_processed_summary/GWAS_Annot/ldmat_9m_p2.mat')); end;
    if ~exist('LDr2_p8sparse', 'var'), load(fullfile(mmil_gwas_path, 'new_9m_SNPs_processed_summary/GWAS_Annot/ldmat_9m_p8.mat')); end;
    result_file = 'w_ld_9m.mat';
end

LDmats = {LDr2_p1sparse, LDr2_p2sparse, LDr2_p8sparse};

defvec = ~isnan(mafvec) & (mafvec > 0.005); defidx = find(defvec);

w_randomprune = zeros(length(defvec), length(LDmats));
for LDmati = 1:length(LDmats)
    LDmat = LDmats{LDmati};
    hits = zeros(size(defvec));
    w_ld = full(sum(LDmat, 2));
    
    % For each SNP $i$ calculate an ww_ld(i) = an average TLD of the SNPs in LD with $i$.
    % If ww_ld(i) is much larger than the LD of $i$ then we expect that $i$
    % will have a higher chance to be picked by random pruning.
    %ww_ld = zeros(size(w_ld));for i=1:length(w_ld), ww_ld(i) = mean(w_ld(find(LDmat(:, i)))); end

    for prune_index = 1:nprune
        si = defidx(randperm(length(defidx)));
        prunevec = false(size(defvec));
        for i = 1:length(si)
          if ~prunevec(si(i))
             prunevec(find(LDmat(:,si(i)))) = true; 
             prunevec(si(i)) = false;
          end
        end
        hits(~prunevec) = hits(~prunevec) + 1;
        fprintf('Matrix %i of %i, iterations %i of %i\n', LDmati, length(LDmats), prune_index, nprune);
    end
    
    w_randomprune(:, LDmati) = (hits ./ nprune);
end

w_ld_p1 = w_randomprune(:, 1);
w_ld_p2 = w_randomprune(:, 2);
w_ld_p8 = w_randomprune(:, 3);

save(result_file, 'w_ld_p1', 'w_ld_p2', 'w_ld_p8');
 

%k=sum(hits)/nprune/sum(1./w_ld);
%for i=0:nprune, b(i+1)=mean(1./w_ld(hits==i)); end; plot((0:nprune)/nprune, b); 
%xlabel('w_j = #\{t : j\in J_t\}/T');
%ylabel('mean(1 / w_j)');
