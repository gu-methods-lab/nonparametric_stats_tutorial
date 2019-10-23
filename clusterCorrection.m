function results = clusterCorrection(surf_ds,vo,fo,v_inf)

% Load surface neighborhood structure
load('/Users/srd49/Untitled Folder/sample_data/cluster_nbrhood.mat');
nbrhood_mat = cosmo_convert_neighborhood(cluster_nbrhood,'matrix');

% Load subcortical masks to remove meaningless vertices
subcorticalMsk_dir = '/Users/srd49/Untitled Folder/sample_data';
fid = fopen(fullfile(subcorticalMsk_dir,'lh_subcortical_64.1D.roi'));
lhMask = textscan(fid,'%f%f');
fclose(fid);
fid = fopen(fullfile(subcorticalMsk_dir,'rh_subcortical_64.1D.roi'));
rhMask = textscan(fid,'%f%f');
fclose(fid);
subcortMsk = [lhMask{1}+1;(rhMask{1}+1+size(v_inf,1)/2)];


%Get p-values Data
%parfor p = 1:size(surf_ds.samples,2) % iterate through each permutation
%    % for a given permutation calculate the p-value at each surface node. This p-value represents the percentage of values in the null distribution that are greater than the observed statistic
%    pVal(:,p) = sum(surf_ds.samples(:,p)<=surf_ds.samples(:,1:end),2)/(size(surf_ds.samples,2));
%end
%% remove effect of useless vertices (e.g. subcortical vertices)
%pVal(subcortMsk,:) = 1;
%pVal(isnan(surf_ds.samples)) = 1;
load('/Users/srd49/Untitled Folder/sample_data/p_val.mat'); %pre-computed to save time


% Identify clusters in observed data and build null distribution of maximum cluster sizes
for p = 1:size(surf_ds.samples,2)
    tmp = cosmo_clusterize(double(pVal(:,p)<=.001)',nbrhood_mat); % calls on cosmo_clusterize in the CosmoMVPA toolbox
    if ~isempty(tmp) % if there are clusters
        for c = 1:length(tmp)
            tmp2(c) = sum(surf_ds.samples(tmp{c},p)); % get the mass of each cluster (sum of the statistics within the cluster)
        end
        if p==1 % if it's our actual observed data then keep all the clusters and their corresponding cluster masses
            clusters  =tmp;
            clusMass = tmp2;
        else % if it's a permuted distribution, the identify the mass of the larges cluster
            clusPerm{p-1} = tmp;
            [clusMassPerm(p-1),idx] = max(tmp2); % keeps track of mass of the largest cluster
        end
    else
        if p==1
            clusters  =[];
            clusMass = 0;
        else
            clusPerm{p-1} = [];
            clusMassPerm(p-1) = 0;
        end
        
    end
    clear tmp2 tmp
end

% Identify the significance of each of our observed clusters by comparing it to the null distribution of largest clusters
for c = 1:length(clusMass)
    p_clus(c) = sum(clusMassPerm>clusMass(c))/length(clusMassPerm);
end

sigIdx = find(p_clus<=.05);
sig_p = p_clus(sigIdx);
sigClus = clusters(sigIdx);
results.sigIdx = sigIdx;
results.sig_p = sig_p;
results.sigClus = sigClus;
results.p_val = pVal;