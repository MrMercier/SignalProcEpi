function [data] = MrM_stats_Pcorrection(cfg, data)
% stats_Pthreshold performs a basic thresholding of significativity
% using different FDR methods or
% based on consecutives significant values along time (using the p-values)
% output create both thresholded p-values and binary mask

% Data should be the output of either:
% MrM_stats_erp_baseline.m / MrM_stats_erp_unpaired.m / MrM_stats_erp_unpaired_sum.m
% MrM_stats_tf_baseline.m / MrM_stats_tf_unpaired.m / MrM_stats_tf_unpaired_sum.m
%
% Use as: MrM_stats_Pcorrection(cfg,data)
% The configuration can have the following parameters
% cfg.channel        = label of the selected channel {''}     
% cfg.time           = selected time period (extra is zeroed in output)
% cfg.freq           = selected frequencies (extra is zeroed in output)
% cfg.parameterField = name of the field on whicth the thresholding will be
%                      performed
%                      'pow_p' 'pci_p' or 'erp_p'
% cfg.method         = 'Bonferroni' / 'Bonferroni-Holm'
%                      'Benjamini & Hochberg' / 'Benjamini & Yekutieli'
%                      'Benjamini, Krieger, & Yekutieli'
% cfg.p              = p value
%                      [0 1] (ex: 0.05)
% cfg.nbtails        = number of tails used in the stats (1 or 2)
%                      selected p value will be modiied in consequence
%                      (generally ERP and power are 2 tailed,
%                      for phase it could be 1 or 2)
% cfg.duration       = number of consecutive time point that should be
%                      within the specified P-value range
%                      [ms]


%time index
switch isfield(cfg, 'time')
    case 1
        if length(cfg.time) == 2
            id_time =[find(round(data.time*1000) == round(cfg.time(1)*1000)) find(round(data.time*1000) == round(cfg.time(2)*1000))];            
        elseif strcmp(cfg.time, 'all');
            id_time = [1 length(data.time)];
        else
            error('time period wrongly stated');
        end
    case 0
        error('no time of interest defined');
end

%freq index
if strcmp(cfg.parameterField, 'pow_p') | strcmp(cfg.parameterField, 'pci_p')
        if length(cfg.freq) == 2
            id_freq =[find(data.freq == cfg.freq(1)) find(data.freq == cfg.freq(2))];            
        elseif strcmp(cfg.freq, 'all');
            id_freq = [1 length(data.freq)];
        else
            error('frequencies wrongly stated');
        end
end

%channels selection
if ischar(cfg.channel)
    id_ch = strmatch(cfg.channel, data.label, 'exact');
elseif ~isfield (cfg, 'channel')
    error('no channel was specified');
else
    error('channel field wrongly specified');
end

%check for the selected p-value and correct it depending on nb of tails
switch isfield(cfg, 'p')
    case 1
        if cfg.nbtails == 2
            cfg.p = cfg.p/2;
        elseif strcmp(cfg.parameterField, 'erp_p') ...
                || strcmp(cfg.parameterField, 'pow_p') ...
                || (strcmp(cfg.parameterField, 'pci_p') && length(data.cfg.pci_prob) == 2)
            cfg.p = cfg.p/2;
        else            
        end            
    case 0
        error('no p range specified');
end

% data subset selection
if strcmp(cfg.parameterField, 'erp_p')
    tmp.p    = data.erp_p(id_ch,id_time(1):id_time(2));
elseif strcmp(cfg.parameterField, 'pow_p')
    tmp.p    = data.pow_p(id_ch,id_freq(1):id_freq(2),id_time(1):id_time(2));
elseif strcmp(cfg.parameterField, 'pci_p')
    tmp.p    = data.pci_p(id_ch,id_freq(1):id_freq(2),id_time(1):id_time(2));
end
%_____________________
% Correction methods:
if isfield(cfg, 'method')
    if strcmp(cfg.method, 'Bonferroni')
        cfg.p = cfg.p / numel(tmp.p);
        tmp.mask   = tmp.p < cfg.p;
        tmp.p      = tmp.p .* numel(tmp.p);
    elseif strcmp(cfg.method, 'Bonferroni-Holm')
        [tmp.p, tmp.mask] = MrM_HBcorr(tmp.p, cfg.alpha);
%          [tmp.p, tmp.mask] = bonf_holm(tmp.p, cfg.alpha);        
    elseif strcmp(cfg.method, 'Benjamini & Hochberg')
        [tmp.mask, cfg.crit_p, tmp.p ] = fdr_bh(tmp.p, cfg.alpha, 'pdep', cfg.report);
    elseif strcmp(cfg.method, 'Benjamini & Yekutieli')
        [tmp.mask, cfg.crit_p, tmp.p ] = fdr_bh(tmp.p, cfg.alpha, 'dep', cfg.report);
    elseif strcmp(cfg.method, 'Benjamini, Krieger, & Yekutieli')
        [tmp.mask, cfg.crit_p]= fdr_bky(tmp.p,cfg.alpha,cfg.report);
    else
    end
else
%     % p-value threshold:
%     tmp.mask = tmp.p < cfg.p;
end
% p-value threshold:
tmp.mask = tmp.p < cfg.p;


%______________________
% Duration threshold:
%define time index length base on the sampling
switch isfield(cfg, 'duration')
    case 1
        Tstep   = data.time(2)*1000 - data.time(1)*1000;
        Twindow = round(cfg.duration/Tstep);
        if rem(cfg.duration, round(Tstep)) ~= 0;
            error('duration threshold define do not match the step resolution of your data (check the time field)');
        else

        end
        %scan for consecutive significant time points
        for ch = 1:length(id_ch)
            %time domain case
            if ~isfield(data, 'freq')
                tmp.mask_Dt = zeros(1,length(tmp.mask));
                for t = 1:size(tmp.mask,2)- (Twindow-1)
                    if all(tmp.mask(1, t:t+(Twindow-1))==1)
                        tmp.mask_Dt(1, t:t+(Twindow-1)) = 1;
                    else
                    end
                end
            %time-frequency space case
            elseif isfield(data, 'freq')
                tmp.mask_Dt = zeros(1,size(tmp.mask,2), size(tmp.mask,3));
                for f = 1:length(data.freq)
                    for t = 1:size(tmp.mask,3)- (Twindow-1)
                        if all(tmp.mask(1, f, t:t+(Twindow-1))==1)
                            tmp.mask_Dt(1, f, t:t+(Twindow-1)) = 1;
                        else
                        end
                    end
                end
            end
        end
        tmp.mask = tmp.mask_Dt;
    case 0
end
%_______________
%prepare output
if strcmp(cfg.parameterField, 'erp_p')
    mask_out                                 = zeros(size(data.erp_p));
    mask_out(id_ch, id_time(1):id_time(2))   = tmp.mask;
    data.erp_mask_threshold                  = mask_out;
    Pcorr_out                                = zeros(size(data.erp_p));    
    Pcorr_out(id_ch, id_time(1):id_time(2))  = tmp.p;
    data.erp_p_threshold                     = Pcorr_out;
elseif strcmp(cfg.parameterField, 'pow_p')
    mask_out                                                        = zeros(size(data.pow_p));
    mask_out(id_ch, id_freq(1):id_freq(2), id_time(1):id_time(2))    = tmp.mask;    
    data.pow_mask_threshold                                         = mask_out;
    Pcorr_out                                                       = zeros(size(data.pow_p));
    Pcorr_out(id_ch, id_freq(1):id_freq(2), id_time(1):id_time(2))  = tmp.p;    
    data.pow_p_threshold                                            = Pcorr_out;
elseif strcmp(cfg.parameterField, 'pci_p')
    mask_out                                                       = zeros(size(data.pci_p));
    mask_out(id_ch, id_freq(1):id_freq(2), id_time(1):id_time(2))  = tmp.mask;     
    data.pci_mask_threshold                                        = mask_out;
    Pcorr_out                                                      = zeros(size(data.pci_p));    
    Pcorr_out(id_ch, id_freq(1):id_freq(2), id_time(1):id_time(2))                     = tmp.p;     
    data.pci_p_threshold                                           = Pcorr_out;
end
data.cfg_correction = cfg;
end
