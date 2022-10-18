function [out, t] = ccg(Spike, bin_size, window_size, show_full, calculate_fr)
if nargin < 2 || isempty(bin_size)
    bin_size = 0.0005;
end
if nargin < 3 || isempty(window_size)
    window_size = 0.004;
end
if nargin < 4 || isempty(show_full)
    show_full = false;
end
if nargin < 5 || isempty(calculate_fr)
    calculate_fr = false;
end

sampRate = Spike.P.sample_rate;
binSize = bin_size * sampRate;
window = window_size * sampRate;

bin = 0 : binSize : window;
nBin = length(bin) - 1;
nClu = length(Spike.time);

sp = floor(cell2mat(Spike.time) * sampRate);
spidx = cell2mat(cellfun(@(x, y) y * ones(length(x), 1), Spike.time, num2cell((1:nClu)'), 'UniformOutput', false));
spike = sortrows([sp, spidx]);
nSpikeClu = cellfun(@length, Spike.time);
nSpike = size(spike, 1);

out = zeros(nClu, nClu, nBin);
for shift = 1:nSpike - 1
    dsp = spike(shift + 1:end, 1) - spike(1:end - shift, 1);

    binIdx = discretize(dsp, bin);
    inBin = ~isnan(binIdx);
    if sum(inBin) == 0; break; end

    i_clu = spike(1:end - shift, 2);
    j_clu = spike(shift + 1:end, 2);
    idx = sub2ind(size(out), i_clu(inBin), j_clu(inBin), binIdx(inBin));
    uidx = unique(idx);
    if length(uidx) == 1
        out(uidx) = out(uidx) + length(idx);
    else
        count_idx = histc(idx, uidx);
        out(uidx) = out(uidx) + count_idx;
    end
end

if nargout == 2
    if show_full
        t = -window_size + bin_size/2 : bin_size : window_size - bin_size/2;
    else
        t = bin_size / 2 : bin_size : window_size - bin_size/2;
    end
end

if show_full
    other_half = flip(permute(out, [2, 1, 3]), 3);
    out = cat(3, other_half, out);
end

if calculate_fr
    out = out ./ (nSpikeClu * bin_size);
end
