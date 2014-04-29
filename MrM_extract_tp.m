function [lat] = MrM_extract_tp(cfg, data)
% extract_tp extract the time point that see a variable of interest
% reaching a specified value
% this could be a p-value or based on the mask
% the nested use of stats_Pthreshold function require corresponding config

% Data should be the output of either:
% stats_erp_baseline.m / stats_erp_unpaired.m / stats_erp_unpaired_sum.m
% stats_tf_baseline.m / stats_tf_unpaired.m / stats_tf_unpaired_sum.m

% Use as:
% extract_tp(cfg,data)
% The configuration can have the following parameters
% cfg.channel        = selected channel (string format)
% cfg.time           = time period over which the stats will be performed
%                     [sec sec] or 'all'
% cfg.freq           = selected freq
%                     [Hz Hz] or 'all'
% cfg.var            = variable of interest
%                     'erp' 'pow' 'plf'
% cfg.type           = data type
%                     'Baseline' or 'Difference'
% cfg.toi            = 'first', 'absAmpl', 'maxAmpl', 'minAmpl'
%                      (first onset or absolute maximum or maximum or minimum in amplitude)
%                      future implementation mean or median
% cfg.correction     = 'no' or 'yes' in this case check stats_Pthreshold
%                      function for related config options and requirement

%channels selection
if strcmp(cfg.var, 'erp')
    if strcmp(cfg.channel, 'all')
        id_ch = 1:length(data.label);
    elseif ischar(cfg.channel)
        for c=1:length(cfg.channel)
            id_ch = strmatch(cfg.channel, data.label, 'exact');
        end
    elseif ~isfield (cfg, 'channel')
        error('no channel was specified');
    else
        error('channel field wrongly specified');
    end
else strcmp(cfg.var, 'pow') || strcmp(cfg.var, 'plf');
    id_ch = 1;
end

%time index
switch isfield(cfg, 'time')
    case 1
        if length(cfg.time) == 2
            id_time =[find(round(data.time*1000) == round(cfg.time(1)*1000)) find(round(data.time*1000) == round(cfg.time(2)*1000))];
        elseif strcmp(cfg.time, 'all');
            id_time = [data.time{1}(1) data.time{1}(end)];
        else
            error('time period wrongly stated');
        end
    case 0
        error('no time of interest defined');
end
%frequency index
switch isfield(cfg, 'freq')
    case 1
        if length(cfg.freq) == 2
            id_freq = [find(data.freq == cfg.freq(1)) find(data.freq == cfg.freq(2))];
        elseif strcmp(cfg.freq, 'all');
            id_freq = [1 length(data.freq)];
        else
            error('frequencies wrongly stated');
        end
    case 0
end

if strcmp(cfg.correction, 'yes') && strcmp(cfg.type, 'Baseline')
        data = MrM_stats_Pcorrection(cfg, data);
        if strcmp(cfg.var, 'erp');
        tmp = data.erp_mask_threshold .* data.erp;
        elseif strcmp(cfg.var, 'pow');
        tmp = squeeze(data.pow_mask_threshold) .* squeeze(data.powspctrm);    
        elseif strcmp(cfg.var, 'plf');
        tmp = squeeze(data.plf_mask_threshold) .* squeeze(data.plfspctrm);
        end
elseif strcmp(cfg.correction, 'no') && strcmp(cfg.type, 'Baseline')
        if strcmp(cfg.var, 'erp');
        tmp = data.erp_mask .* data.erp;
        elseif strcmp(cfg.var, 'pow');
        tmp = squeeze(data.pow_mask) .* squeeze(data.powspctrm);    
        elseif strcmp(cfg.var, 'plf');
        tmp = squeeze(data.plf_mask) .* squeeze(data.plfspctrm);
        end
elseif strcmp(cfg.correction, 'yes') && strcmp(cfg.type, 'Difference')
        data = MrM_stats_Pcorrection(cfg, data);
        if strcmp(cfg.var, 'erp');
        tmp = data.erp_mask_threshold .* data.diff;
        elseif strcmp(cfg.var, 'pow');
        tmp = squeeze(data.pow_mask_threshold) .* squeeze(data.diff_pow);    
        elseif strcmp(cfg.var, 'plf');
        tmp = squeeze(data.plf_mask_threshold) .* squeeze(data.diff_plf);
        end
elseif strcmp(cfg.correction, 'no') && strcmp(cfg.type, 'Difference')
        if strcmp(cfg.var, 'erp');
        tmp = data.erp_mask .* data.diff;
        elseif strcmp(cfg.var, 'pow');
        tmp = squeeze(data.pow_mask) .* squeeze(data.diff_pow);    
        elseif strcmp(cfg.var, 'plf');
        tmp = squeeze(data.plf_mask) .* squeeze(data.diff_plf);
        end        
end

for c = length(id_ch)
if strcmp(cfg.toi, 'first')
    if strcmp(cfg.var, 'erp')
        lat{c,1} = data.label(id_ch);
        if length(find(tmp(id_ch, id_time(1):id_time(2)),1, 'first')) ==1;
            [e,t]    = find(tmp(id_ch, id_time(1):id_time(2)),1,'first');
            lat{c,2} = data.time(t+id_time(1)); % no -1
        else
            lat{c,2} = NaN;
        end
    elseif strcmp(cfg.var, 'pow') || strcmp(cfg.var, 'plf')
        lat{c,1} = num2str(data.cfg.channel);
        if length(find(tmp(id_freq(1):id_freq(2), id_time(1):id_time(2)),1, 'first')) == 1;
            [f,t]    = find(tmp(id_freq(1):id_freq(2), id_time(1):id_time(2)),1, 'first');
            lat{c,2} = data.time(t+id_time(1)-1);
            lat{c,3} = data.freq(f+id_freq(1)-1);            
        else
            lat{c,2} = NaN;
            lat{c,3} = NaN;
        end
    end
elseif strcmp(cfg.toi, 'absAmpl')
    if strcmp(cfg.var, 'erp')
        lat{c,1} = data.label(id_ch);
        if any(find(tmp(id_ch, id_time(1):id_time(2))));
            [ampl,t]    = max(abs(tmp(id_ch, id_time(1):id_time(2))));
            lat{c,2}    = data.time(t+id_time(1)); %no -1
        else
            lat{c,2} = NaN;
        end
    elseif strcmp(cfg.var, 'pow') || strcmp(cfg.var, 'plf')
        lat{c,1} = num2str(data.cfg.channel);
        if any(find(tmp(id_freq(1):id_freq(2), id_time(1):id_time(2))));
            [f,t]       = find(abs(tmp(id_freq(1):id_freq(2), id_time(1):id_time(2))) == max(max(abs(tmp(id_freq(1):id_freq(2), id_time(1):id_time(2))))));
            lat{c,2}    = data.time(t+id_time(1)-1);
            lat{c,3}    = data.freq(f+id_freq(1)-1);
        else
            lat{c,2} = NaN;
            lat{c,3} = NaN;            
        end
    end    
elseif strcmp(cfg.toi, 'maxAmpl')
    if strcmp(cfg.var, 'erp')
        lat{c,1} = data.label(id_ch);
        if any(find(tmp(id_ch, id_time(1):id_time(2))));
            [ampl,t]    = max(tmp(id_ch, id_time(1):id_time(2)));
            lat{c,2}    = data.time(t+id_time(1)); %no -1
        else
            lat{c,2} = NaN;
        end
    elseif strcmp(cfg.var, 'pow') || strcmp(cfg.var, 'plf')
        lat{c,1} = num2str(data.cfg.channel);
        if any(find(tmp(id_freq(1):id_freq(2), id_time(1):id_time(2))));
            [f,t]       = find(tmp(id_freq(1):id_freq(2), id_time(1):id_time(2)) == max(max(tmp(id_freq(1):id_freq(2), id_time(1):id_time(2)))));
            lat{c,2}    = data.time(t+id_time(1)-1);
            lat{c,3}    = data.freq(f+id_freq(1)-1);
        else
            lat{c,2} = NaN;
            lat{c,3} = NaN;            
        end
    end
elseif strcmp(cfg.toi, 'minAmpl')
    if strcmp(cfg.var, 'erp')
        lat{c,1} = data.label(id_ch);
        if any(find(tmp(id_ch, id_time(1):id_time(2))));
            [ampl,t]    = min(tmp(id_ch, id_time(1):id_time(2)));
            lat{c,2}    = data.time(t+id_time(1)-1);
        else
            lat{c,2} = NaN;
        end
    elseif strcmp(cfg.var, 'pow') || strcmp(cfg.var, 'plf')
        lat{c,1} = num2str(data.cfg.channel);
        if any(find(tmp(id_freq(1):id_freq(2), id_time(1):id_time(2))));
            [f,t]       = find(tmp(id_freq(1):id_freq(2), id_time(1):id_time(2)) == min(min(tmp(id_freq(1):id_freq(2), id_time(1):id_time(2)))));
            lat{c,2}    = data.time(t+id_time(1)-1);
        else
            lat{c,2} = NaN;
            lat{c,3} = NaN;              
        end
    end
else
    error('toi wrongly stated');
end
end
end

