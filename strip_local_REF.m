
function [data_csd] = strip_local_REF(felec,ssize,data,bad_index)

%% FUNCTION TO RE-REFERENCE AN INTRACRANIAL GRID:
% This function re-references each channel to up to two channels on the
% axe (as few as a single channel). It does not use bad
% channels. Those channels are NaNed, as are channels surrounded by bad
% channels.
%
% Inputs:
% felec = the number of the first electrode on the strip
% ssize = the number of electrodes on the strip
% data_erp = a trial's worth of voltage data (electrodes,amplitude
%            include the entire trial, not just data for grid electrodes
% bad_index = the list of bad electrodes
%
% Output:
% data_erp_csd = the rereferenced grid...
%%

grid = reshape(felec:(ssize + felec - 1),ssize,1);
data_csd = zeros(ssize,length(data));
grid = flipud(grid);

for yb = 1:ssize
    if yb == 1
        cnt = 1;
        V1 = zeros(1,length(data));
        E1 = 0;
        V2 = data(grid(yb + 1,1),:);
        E2 = grid(yb + 1,1);
    elseif yb > 1 && yb < ssize
        cnt = 2;
        V1 = data(grid(yb - 1,1),:);
        E1 = grid(yb - 1,1);
        V2 = data(grid(yb + 1,1),:);
        E2 = grid(yb + 1,1);
    elseif yb == ssize
        cnt = 1;
        V1 = data(grid(yb - 1,1),:);
        E1 = grid(yb - 1,1);
        V2 = zeros(1,length(data));
        E2 = 0;
    else error('didn''t match anything');
    end;
    %check for bad Channels and replace their amplitudes by zeros
    if any(bad_index == E1)
        V1 = zeros(1,length(data));
        cnt = cnt - 1;
    end;
    if any(bad_index == E2)
        V2 = zeros(1,length(data));
        cnt = cnt - 1;
    end;
    if any(bad_index == (grid(yb,1)))
        data(grid(yb,1),:) = zeros(1,length(data));
        V = zeros(1,length(data));
    else
        V = data(grid(yb,1),:);
    end;
    
    if cnt == 0
        V = NaN(1,length(data));
        fprintf('channel %d set to NaNs because it has no neighbouring channels\n', grid(yb,1));
        % error('better do something about this');
    end;
    if isempty(find(V,1)) %if this channel was already determined to be bad (only NaNs) don't subract surrounding channels (and create nonzero values)
        fprintf('channel %d set to NaNs because it was set has BAD\n', grid(yb,1));
        data_csd(grid(yb,1),:) = NaN(1,length(data));
    else
        data_csd(grid(yb,1),:) = V - (1/cnt * (V1 + V2));
        fprintf('channel %d is re-referenced\n', grid(yb,1));
    end;
end;

data_csd = data_csd(felec:grid(1,1),:);