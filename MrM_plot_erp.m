function fig = MrM_plot_erp(cfg, varargin)
% ERP single plot
%
% cfg.dataset    = dataset name within the data structure
% cfg.channel    = cell containing selected channel
% cfg.time       = time period to be plot
%                  [sec sec] or 'all'
% cfg.baseline   = 'all' or [sec sec]
% cfg.amplitude  = 'minmax' (default) or [mV mV]
% cfg.ste        = standard error
% cfg.erp_p      = statistical results
% cfg.color_erp  = {'k' 'r' 'b' 'g' 'y' 'm' 'c' 'w'}
% cfg.color_ste  = {'k' 'r' 'b' 'g' 'y' 'm' 'c' 'w'}
% cfg.color_p    = {'k' 'r' 'b' 'g' 'y' 'm' 'c' 'w'}
% cfg.linestyle  = {'-' '--' ':' '-.' 'none'}
% cfg.linesize   = []
% cfg.labsize    = label size;
% cfg.legend     = 'yes' or 'no';
% cfg.legend_erp = {'condition1' 'condition2' ...}
% cfg.legend_p   = {'comparison1' 'comparison2' ...}
% cfg.size       = width x height ex: [1024 768]
% cfg.parent     = panel{1...n};
% cfg.style      = 'trace' 'ribbon'

%check the concistency between input and dataset
if length(varargin) ~= size(cfg.dataset,1)
    error('inconsistencies between the input files and the input dataset');
else
end

%set defaults values
if ~isfield(cfg, 'linesize')
    cfg.linesize = 2;
end

%subset selection for each input dataset
dat = 0;
for i = 1:length(varargin)
    %channel index
    switch isfield(cfg, 'channel')
        case 1
            if iscellstr(cfg.channel)
                id_ch(i) =find(strcmp(cfg.channel, varargin{i}.label));
            else
                error('channel wrongly stated');
            end
        case 0
            error('no channel defined');
    end
    %time index
    switch isfield(cfg, 'time')
        case 1
            if length(cfg.time) == 2 && isnumeric(cfg.time)
                id_time(i,:) =[find(round(varargin{i}.time*1000) == round(cfg.time(1)*1000)) find(round(varargin{i}.time*1000) == round(cfg.time(2)*1000))];
            elseif strcmp(cfg.time, 'all');
                id_time(i,:) = [varargin{i}.time(1+1) varargin{i}.time(end)];
            else
                error('time period wrongly stated');
            end
        case 0
            error('no time of interest defined');
    end
    tmp.time = varargin{i}.time;
    %baseline
    switch isfield(cfg, 'baseline')
        case 1
            if isnumeric(cfg.baseline) && length(cfg.baseline) ==2
                id_Bl(i,:) = [find(round(tmp.time*1000) == round(cfg.baseline(1)*1000)) find(round(tmp.time*1000) == round(cfg.baseline(2)*1000))];
            elseif strcmp(cfg.time, 'all');
                id_Bl(i,:) = [tmp.time(1) tmp.time(end)];
            end
        case 0
    end
    
    %erp, ste and p selection
    for d = 1:size(cfg.dataset,2)
        if ~isempty(cfg.dataset{i,d})
            dat = dat + 1;
            tmp.dat(dat,:)     = varargin{i}.(cfg.dataset{i,d})(id_ch(i), :);
            if isfield(cfg, 'ste')  && ~isempty(cfg.ste{i,d})
                tmp.ste(dat,:)     = varargin{i}.(cfg.ste{i,d})(id_ch(i), :);
            end
            if isfield(cfg, 'baseline')
                bl = mean(tmp.dat(dat,id_Bl(i,1):id_Bl(i,2)),2);
                tmp.dat(dat,:) = tmp.dat(dat,:) - repmat(bl, [1 size(tmp.dat,2)] );
            end
        end
        %select the (1 - p) values
        tmp.erp_p(i,:) = 1-varargin{i}.(cfg.erp_p)(id_ch(i), :);
        tmp.erp_mask(i,:) = varargin{i}.(cfg.erp_mask)(id_ch(i), :);
    end
    %find the origin
    id_time0 = find((tmp.time*1000) == 0);
end

%check the min and max amplitude
ratio = 1.25;
switch isfield(cfg, 'amplitude')
    case 1
        if strcmp(cfg.amplitude, 'minmax')
            Ampl_minmax = [min(min(tmp.dat)) max(max(tmp.dat))];
        elseif isnumeric(cfg.amplitude) && length(cfg.amplitude) ==2
            Ampl_minmax = cfg.amplitude;
        else
            error('amplitude wrongly stated');
        end
    case 0
        Ampl_minmax = [min(min(tmp.dat)) max(max(tmp.dat))]*ratio;
end

%define the display range (axis, labels and ticks)
Xtick     = [(id_time(1,1)):50:(id_time(1,2))];
Xlab      = [cfg.time(1)*1000:50:cfg.time(2)*1000];
Xlab      = mat2cell(Xlab,1,ones(1,length(Xlab)));
Xlab(1)   = {''};
Xlab(end) = {''};

if isfield(cfg,'parent')
    fig = cfg.parent;
elseif isfield(cfg, 'size'),
    fig = figure('Position', [50 50 cfg.size(1) cfg.size(2)]);
end;

% Plot the ERPs with ste
g1 = subplot(3,2,1:4,'parent',fig); hold on;
axis([id_time(i, 1) , id_time(i, 2) , Ampl_minmax(1), Ampl_minmax(2)]);
if isfield(cfg,'title')
    title(cfg.title);
end
set(g1,'XTick',Xtick);
set(g1,'XTickLabel',Xlab, 'fontsize',cfg.labsize,'fontweight','b');
ylabel(' [µV] ','fontsize',cfg.labsize,'fontweight','b');
for j = 1:size(tmp.dat,1)
    plot(g1,tmp.dat(j,:), cfg.color_erp{j}, 'LineWidth',cfg.linesize,'LineStyle',cfg.linestyle{j});
    switch isfield(cfg, 'ste')
        case 1
            plot(g1,tmp.dat(j,id_time(i, 1):id_time(i, 2))+tmp.ste(j,id_time(i, 1):id_time(i, 2)), cfg.color_erp{j}, 'LineStyle', ':');
            plot(g1,tmp.dat(j,id_time(i, 1):id_time(i, 2))-tmp.ste(j,id_time(i, 1):id_time(i, 2)), cfg.color_erp{j}, 'LineStyle', ':');
        case 0
    end
    
end
plot(g1,id_time(1,:),[0 0],'k');
plot(g1,[id_time0 id_time0],Ampl_minmax,'k');

if isfield(cfg, 'legend') && strcmp(cfg.legend, 'yes')
    legend(g1, cfg.legend_erp,'Location', 'NorthEast');
    legend(g1, 'boxoff');
end

%Plot the 1-p values
g2 = subplot(3,2,5:6,'parent',fig); hold on;
if strcmp(cfg.style, 'trace');
    axis([id_time(i, 1) , id_time(i, 2) , 0.975, 1]);
    xlabel(' [ms] ','fontsize',cfg.labsize,'fontweight','b');
    set(g2,'XTick',Xtick);
    set(g2,'XTickLabel',Xlab,'fontsize',cfg.labsize,'fontweight','b');
    ylabel(' [1-p] ','fontsize',cfg.labsize,'fontweight','b');
    set(g2,'YTick',linspace(0.975,1,6));
    set(g2,'YTickLabel',{'0.975' '' '' '' '' '1'},'fontsize',cfg.labsize,'fontweight','b');
    for j = 1:size(tmp.erp_p)
        plot(g2,(tmp.erp_p(j,:) .* tmp.erp_mask(j,:)), cfg.color_p{j}, 'LineWidth',cfg.linesize);
    end
    plot(g2,[id_time0 id_time0],[0.975 1],'k');
    plot(g2,id_time(1,:),[0 0],'k');
    if strcmp(cfg.legend, 'yes')
        legend(g2, cfg.legend_p, 'Location', 'NorthEast');
        legend(g2, 'boxoff');
    end
elseif strcmp(cfg.style, 'ribbon');
    mask_tmp = double(tmp.erp_mask);
    for j = 1:size(mask_tmp,1);
        %tricks for the filled plot
        tmp_t = tmp.time(id_time(j,1):id_time(j,2));
        tmp_m = mask_tmp(j,id_time(j,1):id_time(j,2));
        e = 0;
        for  p = 1:length(tmp_t)-1
            if tmp_m(p) == tmp_m(p+1)
                new_tmp_t(p+e) = tmp_t(p);
                new_tmp_m(p+e) = tmp_m(p);
            elseif tmp_m(p) < tmp_m(p+1)
                new_tmp_t(p+e) = tmp_t(p);
                new_tmp_m(p+e) = tmp_m(p);
                new_tmp_t(p+e+1) = tmp_t(p+1);
                new_tmp_m(p+e+1) = tmp_m(p);
                e=e+1;        
            elseif tmp_m(p) > tmp_m(p+1)
                new_tmp_t(p+e) = tmp_t(p);
                new_tmp_m(p+e) = tmp_m(p);
                new_tmp_t(p+e+1) = tmp_t(p);
                new_tmp_m(p+e+1) = tmp_m(p+1);
                e=e+1;
            end
        end
        new_tmp_m(1) = 0;
        new_tmp_m(end) = 0;
        new_tmp_m = new_tmp_m + (j-1);
        fill(new_tmp_t*1000,new_tmp_m,cfg.color_p{j}, 'EdgeColor', 'none');
        clear new_tmp_t;
        clear new_tmp_m;
    end
    xlabel(' [ms] ','fontsize',cfg.labsize,'fontweight','b');
    axis([1000*tmp_t(1) , 1000*tmp_t(end) , 0,size(mask_tmp,1)]);    
    set(g2,'XTick',[1000*tmp_t(1):50:1000*tmp_t(end)]);
    set(g2,'XTickLabel',Xlab,'fontsize',cfg.labsize,'fontweight','b');
    set(g2, 'YTicklabel',{});
    if strcmp(cfg.legend, 'yes')
        set(gca,'YTick',linspace(0.5,size(tmp.erp_mask,1) -0.5, size(tmp.erp_mask,1)));
        set(gca,'YTickLabel', cfg.legend_p,'fontsize',cfg.labsize);
    end
else
    error('style for stats subplot not specified');
end
end