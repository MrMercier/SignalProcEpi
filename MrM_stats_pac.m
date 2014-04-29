function [stats] = MrM_stats_pac(cfg,data)

% MrM_stats_pac compute for a series of trials phase amplitude coupling
% and performs randomisation statistics (in the trials dimension)
% Data  should be the output of ft_preprocess
% processing steps are:
% - band-pass to get low frequency modulation
% - Hilbert transform to get low freq phase
%
% -band-pass to get high frequency modulation
% -Hilbert transform to get high freq amplitude
% -band-pass to get low freq modulation of the high freq ampl
% -Hilbert transform to get the phase of low freq modulation of the high freq amplthen
%
% compute Phase Amplitude coopling based on PLV
% (compare Low&High frequencies angle difference consistency across trial)

%use as: MrM_stats_pac(cfg, data)
% Configuration need the following parameters
% cfg.channel       = selected channel
% cfg.trig          = trigger
% cfg.time          = time period over which the stats will be performed
%                     [sec sec] or 'all'
% cfg.Lowfreq       = selected Low frequency band
% cfg.Highfreq      = selected High frequency band
%                     [Hz Hz]
% cfg.bpfiltord     = filter order
% cfg.pac_prob      = statistic threshold (one tail since pac index range is [0,1])
%                     [95]
% cfg.iteration     = number of random permutation
% cfg.distribution  = 'yes' save the random distributions
% cfg.savetmp       = 'yes' save intermediate temporary matrices
% cfg.baseline      = good question!!!

%time index
switch isfield(cfg, 'time')
    case 1
        if length(cfg.time) == 2
            id_time =[find(round(data.time{1}*1000) == round(cfg.time(1)*1000)) find(round(data.time{1}*1000) == round(cfg.time(2)*1000))];
        else
            error('time period wrongly stated');
        end
    case 0
        error('no time of interest defined');
end

%band pass to get low frequency modulation
cfg_tmp                 = [];
cfg_tmp.channel         = cfg.channel;
cfg_tmp.bpfilter        = 'yes';
cfg_tmp.bpfreq          = cfg.Lowfreq;
cfg_tmp.bpfiltord       = cfg.bpfiltord;
cfg_tmp.keeptrials      = 'yes';
cfg_tmp.trials          = find(data.trialinfo(:,1) == cfg.trig);
stats.Low_Freq_Bp       = ft_preprocessing(cfg_tmp,data);
% Low_Freq_Bp             = ft_preprocessing(cfg,data);
% if strcmp(cfg.savetmp, 'yes')
%     save('Low_Freq_Bp.mat','Low_Freq_Bp','-v7.3');
%     eval(['save(''Low_Freq_Bp_' cfg.channel '.mat'', ''Low_Freq_Bp'');' ]);
% end
%Hilbert transform to get low freq phase
cfg_tmp             = [];
cfg_tmp.hilbert     = 'angle';
stats.Low_Freq_Bp_Ht  = ft_preprocessing(cfg_tmp, stats.Low_Freq_Bp);
% if strcmp(cfg.savetmp, 'yes')
%     save('Low_Freq_Bp_Ht.mat','Low_Freq_Bp_Ht','-v7.3');
% end
% clear 'Low_Freq_Bp';

%band pass to get high frequency modulation
cfg_tmp                 = [];
cfg_tmp.channel         = cfg.channel;
cfg_tmp.bpfilter        = 'yes';
cfg_tmp.bpfreq          = cfg.Highfreq;
cfg_tmp.bpfiltord       = cfg.bpfiltord;
cfg_tmp.keeptrials      = 'yes';
cfg_tmp.trials          = find(data.trialinfo(:,1) == cfg.trig);
stats.High_Freq_Bp      = ft_preprocessing(cfg_tmp, data);
% if strcmp(cfg.savetmp, 'yes')
%     save('High_Freq_Bp.mat','High_Freq_Bp','-v7.3');
% end
%Hilbert transform to get high freq amplitude
cfg_tmp               = [];
cfg_tmp.hilbert       = 'abs';
stats.High_Freq_Bp_Ht = ft_preprocessing(cfg_tmp, stats.High_Freq_Bp);
% if strcmp(cfg.savetmp, 'yes')
%     save('High_Freq_Bp_Ht.mat','High_Freq_Bp_Ht','-v7.3');
% end
% clear 'High_Freq_Bp';
%band pass to get low freq modulation of the high freq ampl
cfg_tmp                 = [];
cfg_tmp.bpfilter        = 'yes';
cfg_tmp.bpfreq          = cfg.Lowfreq;
cfg_tmp.bpfiltord       = cfg.bpfiltord;
cfg_tmp.keeptrials      = 'yes';
stats.High_Freq_Bp_Ht_Bp= ft_preprocessing(cfg_tmp, stats.High_Freq_Bp_Ht);
% if strcmp(cfg.savetmp, 'yes')
%     save('High_Freq_Bp_Ht_Bp.mat','High_Freq_Bp_Ht_Bp','-v7.3');
% end
% clear 'High_Freq_Bp_Ht';
%Hilbert transform to get the phase of low freq modulation of the high freq ampl
cfg_tmp                     = [];
cfg_tmp.hilbert             = 'angle';
stats.High_Freq_Bp_Ht_Bp_Ht = ft_preprocessing(cfg_tmp, stats.High_Freq_Bp_Ht_Bp);
% if strcmp(cfg.savetmp, 'yes')
%     save('High_Freq_Bp_Ht_Bp_Ht.mat','High_Freq_Bp_Ht_Bp_Ht','-v7.3');
% end
% clear 'High_Freq_Bp_Ht_Bp';

% data selection
stats.Low_freq_Ph      = reshape(cell2mat(stats.Low_Freq_Bp_Ht.trial),size(stats.Low_Freq_Bp_Ht.trial{1},1),size(stats.Low_Freq_Bp_Ht.trial{1},2), size(stats.Low_Freq_Bp_Ht.trial,2));
stats.Low_freq_Ph      = stats.Low_freq_Ph(:,id_time(1):id_time(2),:);
% clear 'Low_Freq_Bp_Ht';
stats.High_freq_AmplPh = reshape(cell2mat(stats.High_Freq_Bp_Ht_Bp_Ht.trial),size(stats.High_Freq_Bp_Ht_Bp_Ht.trial{1},1),size(stats.High_Freq_Bp_Ht_Bp_Ht.trial{1},2), size(stats.High_Freq_Bp_Ht_Bp_Ht.trial,2));
stats.High_freq_AmplPh = stats.High_freq_AmplPh(:,id_time(1):id_time(2),:);
% clear 'High_Freq_Bp_Ht_Bp_Ht';

%compute Ph AMpl coopling based on PLV ...
%(to compare difference consistency across trial)
stats.pac = abs(mean(exp(1i * (stats.Low_freq_Ph - stats.High_freq_AmplPh )),3));

%stats
stats.cfg = cfg;
stats.pac_Rd = nan(size(stats.pac,1),size(stats.pac,2), cfg.iteration);
for i=1:cfg.iteration
    stats.pac_Rd(:,:,i) = abs(mean(exp(1i * ...
        (stats.Low_freq_Ph(:,:,randperm(size(stats.Low_freq_Ph,3))) ...
        - stats.High_freq_AmplPh(:,:,randperm(size(stats.High_freq_AmplPh,3)))))...
        ,3));
     disp(['Randomization #' num2str(i) ' over ' num2str(cfg.iteration)]);
end
% compute threshold
stats.pac_thrshld = prctile(stats.pac_Rd,cfg.pac_prob,3);
% build mask
stats.pac_mask = stats.pac >= stats.pac_thrshld;
% compute statistics
stats.pac_p = nan(size(stats.pac,1), size(stats.pac,2));
    for e = 1:size(stats.pac,1)
        for t = 1:size(stats.pac,2)
        stats.pac_p(e,t) = length(find(stats.pac_Rd(e,t,:) > stats.pac(e,t)))/size(stats.pac_Rd,3);
        end
    end
disp('Randomization done');
end

%% check distribution
% e = 15;  
% t = 900;
% figure; 
%     steps_edges = linspace(min(stats.pac_Rd(e,t,:)), max(stats.pac_Rd(e,t,:)), 101);
%     bar(steps_edges, squeeze(histc(stats.pac_Rd(e,t,:),steps_edges)), 'k');
%     hold on;
%     thr = stats.pac_thrshld(e,t);
%     line([thr thr], [0 max(histc(stats.pac_Rd(e,t,:), steps_edges))], 'Color', 'r');    
%     dat = stats.pac(e,t);
%     line([dat dat], [0 max(histc(stats.pac_Rd(e,t,:), steps_edges))], 'Color', 'b');    
    