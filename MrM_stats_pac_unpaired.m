function [stats] = MrM_stats_pac_unpaired(cfg,data)

% MrM_stats_pac_unpaired compared phase amplitude coupling for two conditions
% using unpaired randomisation statistics (in the trials dimension)
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
% cfg.trig1         = trigger
% cfg.trig2         = trigger
% cfg.time          = time period over which the stats will be performed
%                     [sec sec] or 'all'
% cfg.Lowfreq       = selected Low frequency band
% cfg.Highfreq      = selected High frequency band
%                     [Hz Hz]
% cfg.bpfiltord     = filter order
% cfg.pac_prob      = statistic threshold (two tails)
%                     [2.5 97.5]
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

%band-pass to get LOW frequency modulation
cfg_tmp                 = [];
cfg_tmp.channel         = cfg.channel;
cfg_tmp.bpfilter        = 'yes';
cfg_tmp.bpfreq          = cfg.Lowfreq;
cfg_tmp.bpfiltord       = cfg.bpfiltord;
cfg_tmp.keeptrials      = 'yes';
cfg_tmp.trials          = find(data.trialinfo(:,1) == cfg.trig1);
stats.trig1_Low_Freq_Bp = ft_preprocessing(cfg_tmp,data);
cfg_tmp.trials          = find(data.trialinfo(:,1) == cfg.trig2);
stats.trig2_Low_Freq_Bp = ft_preprocessing(cfg_tmp,data);
%Hilbert transform to get low freq phase
cfg_tmp                     = [];
cfg_tmp.hilbert             = 'angle';
stats.trig1_Low_Freq_Bp_Ht  = ft_preprocessing(cfg_tmp, stats.trig1_Low_Freq_Bp);
stats.trig2_Low_Freq_Bp_Ht  = ft_preprocessing(cfg_tmp, stats.trig2_Low_Freq_Bp);

%band-pass to get HIGH frequency modulation
cfg_tmp                  = [];
cfg_tmp.channel          = cfg.channel;
cfg_tmp.bpfilter         = 'yes';
cfg_tmp.bpfreq           = cfg.Highfreq;
cfg_tmp.bpfiltord        = cfg.bpfiltord;
cfg_tmp.keeptrials       = 'yes';
cfg_tmp.trials           = find(data.trialinfo(:,1) == cfg.trig1);
stats.trig1_High_Freq_Bp = ft_preprocessing(cfg_tmp, data);
cfg_tmp.trials           = find(data.trialinfo(:,1) == cfg.trig2);
stats.trig2_High_Freq_Bp = ft_preprocessing(cfg_tmp, data);
%Hilbert transform to get high freq amplitude
cfg_tmp               = [];
cfg_tmp.hilbert       = 'abs';
stats.trig1_High_Freq_Bp_Ht = ft_preprocessing(cfg_tmp, stats.trig1_High_Freq_Bp);
stats.trig2_High_Freq_Bp_Ht = ft_preprocessing(cfg_tmp, stats.trig2_High_Freq_Bp);
%band pass to get low freq modulation of the high freq ampl
cfg_tmp                 = [];
cfg_tmp.bpfilter        = 'yes';
cfg_tmp.bpfreq          = cfg.Lowfreq;
cfg_tmp.bpfiltord       = cfg.bpfiltord;
cfg_tmp.keeptrials      = 'yes';
stats.trig1_High_Freq_Bp_Ht_Bp= ft_preprocessing(cfg_tmp, stats.trig1_High_Freq_Bp_Ht);
stats.trig2_High_Freq_Bp_Ht_Bp= ft_preprocessing(cfg_tmp, stats.trig2_High_Freq_Bp_Ht);
%Hilbert transform to get the phase of low freq modulation of the high freq ampl
cfg_tmp                     = [];
cfg_tmp.hilbert             = 'angle';
stats.trig1_High_Freq_Bp_Ht_Bp_Ht = ft_preprocessing(cfg_tmp, stats.trig1_High_Freq_Bp_Ht_Bp);
stats.trig2_High_Freq_Bp_Ht_Bp_Ht = ft_preprocessing(cfg_tmp, stats.trig2_High_Freq_Bp_Ht_Bp);

% data selection
stats.trig1_Low_freq_Ph      = reshape(cell2mat(stats.trig1_Low_Freq_Bp_Ht.trial),size(stats.trig1_Low_Freq_Bp_Ht.trial{1},1),size(stats.trig1_Low_Freq_Bp_Ht.trial{1},2), size(stats.trig1_Low_Freq_Bp_Ht.trial,2));
stats.trig1_Low_freq_Ph      = stats.trig1_Low_freq_Ph(:,id_time(1):id_time(2),:);
stats.trig2_Low_freq_Ph      = reshape(cell2mat(stats.trig2_Low_Freq_Bp_Ht.trial),size(stats.trig2_Low_Freq_Bp_Ht.trial{1},1),size(stats.trig2_Low_Freq_Bp_Ht.trial{1},2), size(stats.trig2_Low_Freq_Bp_Ht.trial,2));
stats.trig2_Low_freq_Ph      = stats.trig2_Low_freq_Ph(:,id_time(1):id_time(2),:);

stats.trig1_High_freq_AmplPh = reshape(cell2mat(stats.trig1_High_Freq_Bp_Ht_Bp_Ht.trial),size(stats.trig1_High_Freq_Bp_Ht_Bp_Ht.trial{1},1),size(stats.trig1_High_Freq_Bp_Ht_Bp_Ht.trial{1},2), size(stats.trig1_High_Freq_Bp_Ht_Bp_Ht.trial,2));
stats.trig1_High_freq_AmplPh = stats.trig1_High_freq_AmplPh(:,id_time(1):id_time(2),:);
stats.trig2_High_freq_AmplPh = reshape(cell2mat(stats.trig2_High_Freq_Bp_Ht_Bp_Ht.trial),size(stats.trig2_High_Freq_Bp_Ht_Bp_Ht.trial{1},1),size(stats.trig2_High_Freq_Bp_Ht_Bp_Ht.trial{1},2), size(stats.trig2_High_Freq_Bp_Ht_Bp_Ht.trial,2));
stats.trig2_High_freq_AmplPh = stats.trig2_High_freq_AmplPh(:,id_time(1):id_time(2),:);

%compute Ph AMpl coopling based on PLV ...
%(to compare difference consistency across trial)
stats.trig1_pac = abs(mean(exp(1i * (stats.trig1_Low_freq_Ph - stats.trig1_High_freq_AmplPh )),3));
stats.trig2_pac = abs(mean(exp(1i * (stats.trig2_Low_freq_Ph - stats.trig2_High_freq_AmplPh )),3));
stats.diff_pac  = stats.trig1_pac -stats.trig2_pac;

%stats
stats.cfg = cfg;
stats.Rd_diff_pac = nan(size(stats.diff_pac,1),size(stats.diff_pac,2), cfg.iteration);
All_trig_Low_freq_Ph      = cat(3,stats.trig1_Low_freq_Ph, stats.trig2_Low_freq_Ph);
All_trig_High_freq_AmplPh = cat(3,stats.trig1_High_freq_AmplPh, stats.trig2_High_freq_AmplPh);
nb_tr_1 = size(stats.trig1_Low_freq_Ph,3);

for i=1:cfg.iteration
    Rd_tr = randperm(size(All_trig_Low_freq_Ph,3));
    All_trig_Low_freq_Ph       = All_trig_Low_freq_Ph(:,:,Rd_tr);
    All_trig_High_freq_AmplPh  = All_trig_High_freq_AmplPh(:,:,Rd_tr);
    
    Rd_pool1_pac = abs(mean(exp(1i * ( ...
        All_trig_Low_freq_Ph(:,:,1:nb_tr_1) ...
        - All_trig_High_freq_AmplPh(:,:,1:nb_tr_1) ...
                                          )),3));
    Rd_pool2_pac = abs(mean(exp(1i * ( ...
        All_trig_Low_freq_Ph(:,:,1+nb_tr_1 :end) ...
        - All_trig_High_freq_AmplPh(:,:,1+nb_tr_1 :end) ...
                                          )),3));
    stats.Rd_diff_pac(:,:,i)  = Rd_pool1_pac - Rd_pool2_pac;    
    disp(['Randomization #' num2str(i) ' over ' num2str(cfg.iteration)]);
end
% compute threshold
stats.pac_thrshld = prctile(stats.Rd_diff_pac,cfg.pac_prob,3);
% build mask
stats.pac_mask = (stats.diff_pac <= stats.pac_thrshld(:,:,1)) | (stats.diff_pac >= stats.pac_thrshld(:,:,2));
% compute statistics
stats.pac_p = nan(size(stats.diff_pac,1), size(stats.diff_pac,2));
    for e = 1:size(stats.diff_pac,1)
        for t = 1:size(stats.diff_pac,2)
           if stats.diff_pac(e,t) >= mean(stats.Rd_diff_pac(e,t,:),3);
           stats.pac_p(e,t) = length(find(stats.Rd_diff_pac(e,t,:) > stats.diff_pac(e,t)))/size(stats.Rd_diff_pac,3);
           elseif stats.diff_pac(e,t) <= mean(stats.Rd_diff_pac(e,t,:),3);
           stats.pac_p(e,t) = length(find(stats.Rd_diff_pac(e,t,:) < stats.diff_pac(e,t)))/size(stats.Rd_diff_pac,3);
           end
        end
    end
disp('Randomization done');
end

% %% check distribution
% e = 12;  
% t = 671;
% figure; 
%     steps_edges = linspace(min(stats.Rd_diff_pac(e,t,:)), max(stats.Rd_diff_pac(e,t,:)), 101);
%     bar(steps_edges, squeeze(histc(stats.Rd_diff_pac(e,t,:),steps_edges)), 'k');
%     hold on;
%     thr_low  = stats.pac_thrshld(e,t,1);
%     thr_high = stats.pac_thrshld(e,t,2);
%     line([thr_low thr_low], [0 max(histc(stats.Rd_diff_pac(e,t,:), steps_edges))], 'Color', 'r', 'LineStyle', ':');
%     line([thr_high thr_high], [0 max(histc(stats.Rd_diff_pac(e,t,:), steps_edges))], 'Color', 'r', 'LineStyle', ':');    
%     dat = stats.diff_pac(e,t);
%     line([dat dat], [0 max(histc(stats.Rd_diff_pac(e,t,:), steps_edges))], 'Color', 'b');    
% %     