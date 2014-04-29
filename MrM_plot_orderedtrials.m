function fig = MrM_plot_orderedtrials(cfg, data)
% plot single trials ordered as a function of a given value
%
% cfg.channel       = cell containing the selected channel
% cfg.time          = time period to be plot
%                     [sec sec] or 'all'
% cfg.freq          = frequency(ies) to be plot [Hz] or [Hz Hz]
%                     (if a range is specified the mean will be computed
% cfg.parameter     = 'Amplitude', 'Power' or 'Phase'
% cfg.sort          = 'RT', 'Amplitude, 'Power', 'Phase' or any matrix
%                     RTs must be present in trial info,
%                     external matrice should match in dimension
% cfg.sort_tf       = time point at which Amplitude/Power/Phase is used to sort the trial TODO !!!
% cfg.event         = specify triggers corresponding to selected trials
%                     (refer to trialinfo field)
% cfg.response      = specify trigger code for response, used to plot RT
% cfg.baseline      = 'no' or time period used to baseline correct the data
%                     [sec sec] (do not apply on phase)
% cfg.baselinetype  = 'relative' or 'absolute' or 'relchange' ('complexdiff' not  ...)
% cfg.zlim          = [Zmin Zmax] 'minmax' (default for Amplitude and Power)
% cfg.smooth        = number of trials used for the smoothing
% cfg.size          = width x height ex: [1024 768]
% cfg.title         = figure title
% cfg.colormap      = colormap (jet, hsv ... )
% cfg.parent        = panel{1...n};

%channel index
switch isfield(cfg, 'channel')
    case 1
        id_ch = strcmp(data.label, cfg.channel);
    case 0
        error('channel not found');
end

%time index
switch isfield(cfg, 'time')
    case 1
        if length(cfg.time) == 2 && strcmp(cfg.parameter, 'Amplitude')
            id_time =[find(round(data.time{1}*1000) == round(cfg.time(1)*1000)) find(round(data.time{1}*1000) == round(cfg.time(2)*1000))];
            id_time0 =find(round(data.time{1}*1000) == 0)-id_time(1)+1;
            t_step = (data.time{1}(2)-data.time{1}(1)) * 1000;
            if strcmp(cfg.sort, 'Amplitude');
                id_sort_tf = find(round(data.time{1}*1000) == round(cfg.sort_tf*1000))-id_time(1)+1;
            end         
        elseif length(cfg.time) == 2 && (strcmp(cfg.parameter,'Power') || strcmp(cfg.parameter,'Phase'))
            id_time =[find(round(data.time*1000) == round(cfg.time(1)*1000)) find(round(data.time*1000) == round(cfg.time(2)*1000))];
            id_time0 =find(round(data.time*1000) == 0)-id_time(1)+1;
            t_step = (data.time(2)-data.time(1)) * 1000;
            if strcmp(cfg.sort, 'Power') || strcmp(cfg.sort, 'Phase');
                id_sort_tf = find(round(data.time*1000) == round(cfg.sort_tf*1000))-id_time(1)+1;
            end
        elseif strcmp(cfg.time, 'all');
            id_time = [1 length(data.time)];
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
        elseif length(cfg.freq) == 1
            id_freq = find(data.freq == cfg.freq(1));
        elseif strcmp(cfg.freq, 'all');
            id_freq = [1 length(data.freq)];
        end
    case 0
end

%Baseline index
switch isfield(cfg, 'baseline');
    case 1
        if isnumeric(cfg.baseline) && length(cfg.baseline) == 2 && strcmp(cfg.parameter, 'Amplitude');
            if round(1000*cfg.baseline(1))>= round(1000*data.time{1}(id_time(1))) && round(1000*cfg.baseline(2))<= round(1000*data.time{1}(id_time(2)))
                id_Bl = [find(round(data.time{1}*1000) == round(cfg.baseline(1)*1000)) find(round(data.time{1}*1000) == round(cfg.baseline(2)*1000))];
            else
                error('baseline period should be included into the selected time period');
            end
        elseif isnumeric(cfg.baseline) && length(cfg.baseline) == 2 && strcmp(cfg.parameter,'Power');
            if round(1000*cfg.baseline(1))>= round(1000*data.time(id_time(1))) && round(1000*cfg.baseline(2))<= round(1000*data.time(id_time(2)))
                id_Bl = [find(round(data.time*1000) == round(cfg.baseline(1)*1000)) find(round(data.time*1000) == round(cfg.baseline(2)*1000))];
            else
                error('baseline period should be included into the selected time period');
            end
        else
            error('baseline wrongly stated');
        end
    case 0
end

% data selection and compute requested parameter

if strcmp(cfg.parameter, 'Amplitude')
    
    %data selection
    dataS = reshape(cell2mat(data.trial), size(data.trial{1}, 1), size(data.trial{1}, 2), size(data.trial,2));
    dataS = dataS(:,:,data.trialinfo(:,1) == cfg.event);
    dat_tmp =  squeeze(dataS(id_ch,id_time(1):id_time(2),:));
    
    if ~strcmp(cfg.baseline, 'no')
        %compute the average amplitude during the baseline period
        dat_ampl_Bl = squeeze(nanmean(dataS(id_ch,id_Bl(1):id_Bl(2),:),2)); bug
        %compute baseline correction
        if strcmp(cfg.baselinetype,'absolute')
            dat_tmp = dat_tmp - repmat(dat_ampl_Bl, [1 size(dat_tmp,1)])';
        elseif strcmp(cfg.baselinetype,'relative')
            dat_tmp = dat_tmp ./ repmat(dat_ampl_Bl, [1 size(dat_tmp,1)])';
        elseif strcmp(cfg.baselinetype,'relchange')
            dat_tmp = (dat_tmp - repmat(dat_ampl_Bl, [1 size(dat_tmp,1)])')./ repmat(dat_ampl_Bl, [1 size(dat_tmp,1)])';
        end
    end
    dat_tmp = dat_tmp';
    
elseif strcmp(cfg.parameter, 'Power')
    dat_tmp = abs(squeeze(data.fourierspctrm(:, id_ch,id_freq(1):id_freq(2),id_time(1):id_time(2)))) .^2;
    
    if ~strcmp(cfg.baseline, 'no')
        %compute the average power during the baseline period
        dat_pow_Bl = squeeze(nanmean(abs(squeeze(data.fourierspctrm(:,id_ch,id_freq(1):id_freq(2),id_Bl(1):id_Bl(2)))) .^2,3));
        %compute baseline correction
        if strcmp(cfg.baselinetype,'absolute')
            dat_tmp = dat_tmp - repmat(dat_pow_Bl, [1 1 size(dat_tmp,3)]);
        elseif strcmp(cfg.baselinetype,'relative')
            dat_tmp = dat_tmp ./ repmat(dat_pow_Bl, [1 1 size(dat_tmp,3)]);
        elseif strcmp(cfg.baselinetype,'relchange')
            dat_tmp = (dat_tmp - repmat(dat_pow_Bl, [1 1 size(dat_tmp,3)]))./ repmat(dat_pow_Bl, [1 1 size(dat_tmp,3)]);
        end
        
        %average across freq if necessary
        if length(cfg.freq) == 2
            dat_tmp = squeeze(mean(dat_tmp,2));
        end
    end
    
elseif strcmp(cfg.parameter, 'Phase')
    %average across freq if necessary
    if length(cfg.freq) ==1
        dat_tmp = squeeze(angle(data.fourierspctrm(:, id_ch,id_freq(1),id_time(1):id_time(2))));
    elseif length(cfg.freq) == 2
        dat_tmp = squeeze(angle(mean(data.fourierspctrm(:, id_ch,id_freq(1):id_freq(2),id_time(1):id_time(2))./abs(data.fourierspctrm(:, id_ch,id_freq(1):id_freq(2),id_time(1):id_time(2))),3)));
    end
end

%compute the order
if strcmp(cfg.sort, 'RT')
    [ord ord_idx] = sort(data.trialinfo(data.trialinfo(:,1) == cfg.event,2));
% elseif strcmp(cfg.sort, 'Amplitude') ...
%         || strcmp(cfg.sort, 'Power') ...
%         || strcmp(cfg.sort, 'Phase');
%         [ord ord_idx] = sort(dat_tmp(:,id_sort_tf));
%     % prepare the response plot
%     if isfield(cfg.response)
%         resp = data.trialinfo(data.trialinfo(:,1) == cfg.response,2);
%         resp = resp(ord_idx);
%     end
% else
%     %use external matrice, verify dim
%     if size(dat_tmp,1) == size(cfg.sort,1);
%     [ord ord_idx] = sort(cfg.sort);
%     else
%         error('the dimension of the input matrice used to sort the trials does not match the number of trials in the data');
%     end
end

%re-order the trials
dat_ord = dat_tmp(ord_idx,:);

%smoothing
switch isfield(cfg, 'smooth')
    case 1
%     dat_ord = filter(ones(1,cfg.smooth)/cfg.smooth,1,dat_ord);
     for a = 1:size(dat_ord,2);
         dat_ord(:,a) = smooth(dat_ord(:,a),cfg.smooth,'moving');
     end
    case 0
end

%check the min and max amplitude
if strcmp(cfg.parameter, 'Amplitude') || strcmp(cfg.parameter, 'Power')
        if strcmp(cfg.zlim, 'minmax')
            Ampl_minmax = [min(min(dat_ord)) max(max(dat_ord))];
        elseif isnumeric(cfg.zlim) && length(cfg.zlim) ==2
            Ampl_minmax = cfg.zlim;
        else
            error('amplitude wrongly stated');
        end
elseif strcmp(cfg.parameter, 'Phase')
    Ampl_minmax = [-pi/2 pi/2];
else
        Ampl_minmax = [min(min(dat_ord)) max(max(dat_ord))]*1.25;
end

%Figure
%surface(dat_ord,'EdgeColor','none','facecolor','interp'); hold on;
%set(gcf,'Renderer','zbuffer');

if isfield(cfg,'parent');
    fig = cfg.parent;
elseif isfield(cfg, 'size')
    fig=figure('Position',[50 50 cfg.size(1) cfg.size(2)]);
end

hold on
%define the display range (axis, labels and ticks)
xlabel(' [ms] ','fontsize',12,'fontweight','b');
ylabel(' sorted trials # ','fontsize',12);
axis([0 ((abs(cfg.time(1)) + abs(cfg.time(2))) * 1000 / t_step) 0 size(dat_ord,1)]);
Xtick       = [1:100/t_step:((abs(cfg.time(1)) + abs(cfg.time(2))) * 1000 / t_step)];
set(gca, 'XTick', Xtick);
Xlab        = [cfg.time(1)*1000:100:cfg.time(2)*1000];
Xlab        = mat2cell(Xlab,1,ones(1,length(Xlab)));
Xlab(1)     = {''};
Xlab(end)   = {''};
set(gca, 'XTickLabel', Xlab);

%plot
imagesc(dat_ord, [Ampl_minmax(1) Ampl_minmax(2)]);

set(gcf,'colormap',cfg.colormap);
plot([id_time0 id_time0], [0 size(dat_ord,1)], 'k','LineStyle', '--');
if strcmp(cfg.sort, 'RT') || isnumeric(cfg.sort)
    plot(ord'*1000/ t_step+id_time0, 1:size(dat_ord,1), 'k', 'LineWidth',2);
elseif isfield(cfg.response)
    plot(resp'*1000/ t_step+id_time0, 1:size(dat_ord,1), 'k', 'LineWidth',2);
end
title(cfg.channel);
end