
function [data_csd] = grid_local_REF(felec,xsize,ysize,data,bad_index)

%% FUNCTION TO RE-REFERENCE AN INTRACRANIAL GRID:
% This function re-references each channel to up to four channels on the
% 'x' and 'y' axes (as few as a single channel).
% It does not use bad channels. Those channels are NaNed,
% as are channels surrounded by bad channels.
%
% Inputs:
% felec = the number of the first electrode on the grid
% xsize = the number of electrodes on the grids x-coordinate
% ysize = the number of electrodes on the grids y_coordinate
% data_erp = a trial's worth of voltage data (electrodes x time, amplitude
%            include the entire trial, not just data for grid electrodes
% bad_index = the list of bad electrodes
%
% Output:
% data_erp_csd = the rereferenced grid...
%%

grid = reshape(felec:((xsize*ysize)+ felec - 1),xsize,ysize);
data_csd = zeros((xsize*ysize),length(data));

for xb = 1:xsize
    for yb = 1:ysize
        if xb == 1 && yb == 1
            cnt = 2;
            V1 = zeros(1,length(data));
            E1 = 0;
            V2 = data(grid(xb,yb + 1),:);
            E2 = grid(xb,yb + 1);
            V3 = data(grid(xb + 1,yb),:);
            E3 = grid(xb + 1,yb);
            V4 = zeros(1,length(data));
            E4 = 0;
        elseif xb == 1 && yb == ysize
            cnt = 2;
            V1 = zeros(1,length(data));
            E1 = 0;
            V2 = zeros(1,length(data));
            E2 = 0;
            V3 = data(grid(xb,yb - 1),:);
            E3 = grid(xb,yb - 1);
            V4 = data(grid(xb + 1,yb),:);
            E4 = grid(xb + 1,yb);
        elseif xb == xsize && yb == 1
            cnt = 2;
            V1 = data(grid(xb - 1,yb),:);
            E1 = grid(xb - 1,yb);
            V2 = data(grid(xb,yb + 1),:);
            E2 = grid(xb,yb + 1);
            V3 = zeros(1,length(data));
            E3 = 0;
            V4 = zeros(1,length(data));
            E4 = 0;
        elseif xb == xsize && yb == ysize
            cnt = 2;
            V1 = data(grid(xb - 1,yb),:);
            E1 = grid(xb - 1,yb);
            V2 = zeros(1,length(data));
            E2 = 0;
            V3 = zeros(1,length(data));
            E3 = 0;
            V4 = data(grid(xb,yb - 1),:);
            E4 = grid(xb,yb - 1);
        elseif xb == 1 && yb > 1 && yb < ysize
            cnt = 3;
            V1 = zeros(1,length(data));
            E1 = 0;
            V2 = data(grid(xb,yb + 1),:);
            E2 = grid(xb,yb + 1);
            V3 = data(grid(xb + 1,yb),:);
            E3 = grid(xb + 1,yb);
            V4 = data(grid(xb,yb - 1),:);
            E4 = grid(xb,yb - 1);
        elseif xb == xsize && yb > 1 && yb < ysize
            cnt = 3;
            V1 = data(grid(xb - 1,yb),:);
            E1 = grid(xb - 1,yb);
            V2 = data(grid(xb,yb + 1),:);
            E2 = grid(xb,yb + 1);
            V3 = zeros(1,length(data));
            E3 = 0;
            V4 = data(grid(xb,yb - 1),:);
            E4 = grid(xb,yb - 1);
        elseif yb == 1 && xb > 1 && xb < xsize
            cnt = 3;
            V1 = data(grid(xb - 1,yb),:);
            E1 = grid(xb - 1,yb);
            V2 = data(grid(xb,yb + 1),:);
            E2 = grid(xb,yb + 1);
            V3 = data(grid(xb + 1,yb),:);
            E3 = grid(xb + 1,yb);
            V4 = zeros(1,length(data));
            E4 = 0;
        elseif yb == ysize && xb > 1 && xb < xsize
            cnt = 3;
            V1 = data(grid(xb - 1,yb),:);
            E1 = grid(xb - 1,yb);
            V2 = zeros(1,length(data));
            E2 = 0;
            V3 = data(grid(xb + 1,yb),:);
            E3 = grid(xb + 1,yb);
            V4 = data(grid(xb,yb - 1),:);
            E4 = grid(xb,yb - 1);
        elseif yb > 1 && yb < ysize && xb > 1 && xb < xsize
            cnt = 4;
            V1 = data(grid(xb - 1,yb),:);
            E1 = grid(xb - 1,yb);
            V2 = data(grid(xb,yb + 1),:);
            E2 = grid(xb,yb + 1);
            V3 = data(grid(xb + 1,yb),:);
            E3 = grid(xb + 1,yb);
            V4 = data(grid(xb,yb - 1),:);
            E4 = grid(xb,yb - 1);
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
        if any(bad_index == E3)
            V3 = zeros(1,length(data));
            cnt = cnt - 1;
        end;
        if any(bad_index == E4)
            V4 = zeros(1,length(data));
            cnt = cnt - 1;
        end;
        if any(bad_index == (grid(xb,yb)))
            data(grid(xb,yb),:) = zeros(1,length(data));
            V = zeros(1,length(data));
        else
            V = data(grid(xb,yb),:);
        end;
        if cnt == 0
            V = NaN(1,length(data));
            fprintf('channel %d set to NaNs because it has no neighbouring channels\n', grid(xb,yb));
            % error('better do something about this');
        end;
        if isempty(find(V,1)) %if this channel was already determined to be bad (only NaNs) don't subract surrounding channels (and create nonzero values)
            fprintf('channel %d set to zeros because it was set has BAD\n', grid(xb,yb));
            data_csd(grid(xb,yb),:) = NaN(1,length(data));
        else
            data_csd(grid(xb,yb),:) = V - (1/cnt * (V1 + V2 + V3 + V4));
            fprintf('channel %d is re-referenced\n', grid(xb,yb));
        end;
    end;
end;

data_csd = data_csd(felec:grid(xb,yb),:);