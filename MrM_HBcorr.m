function [Pvalues_corr, mask_corr] = MrM_HBcorr(Pvalues, alpha)
% MrM_HBcorr performs correction for multiple comparison
% based on S. Holm 1979;

% p-values to be corrected inputed in the (Pvalues) matrice (2D or 3D),
% and outputed in [data_corr, mask_corr]
% with data_corr being the corrected p-values
% and mask_corr being the mask a significancy (1 correspond to p=<alpha)

% Use as: MrM_HBcorr(cfg,data)
% The configuration can have the following parameters
% cfg.alpha = selected threshold (by Default = 0.05)


if nargin<2,
    alpha=0.05;
elseif alpha<=0,
    error('Alpha must be greater than 0.');
elseif alpha>=1,
    error('Alpha must be less than 1.');
end

s = size(Pvalues);
%depending on the number of dim, reorganise the input p matrice
if isvector(Pvalues)
    if size(Pvalues,1)>1
        Pvalues = Pvalues';
    end
    [Pvalues_1D_sort, ind] = sort(Pvalues);
else
                Pvalues_1D = reshape (Pvalues, 1, prod(s));
    [Pvalues_1D_sort, ind] = sort(Pvalues_1D); 
end

% [dummy, unsort_ids]=sort(ind); %indices to return sorted_p to pvalues order

%compute the correction factor
m_fact = length(Pvalues_1D_sort):-1:1; 
%correction
Pvalues_corr_1D_sort = Pvalues_1D_sort.*m_fact;
%adjust the values to ensure the order
Pvalues_corr_1D_sort(2:length(Pvalues_1D_sort))=max([Pvalues_corr_1D_sort(1:(length(Pvalues_1D_sort)-1)); Pvalues_corr_1D_sort(2:length(Pvalues_1D_sort))]); 

%re-sort in the original order
Pvalues_corr_1D(ind) = Pvalues_corr_1D_sort; 
% Pvalues_corr_1D = Pvalues_corr_1D_sort(unsort_ids);

%reshape with orginal dimensions
Pvalues_corr = reshape(Pvalues_corr_1D, s);
%create the binary mask of significance
mask_corr = Pvalues_corr < alpha;
end
