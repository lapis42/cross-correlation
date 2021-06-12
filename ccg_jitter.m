function out = ccg_jitter(Spike, bin_size, window_size, normalize, show_full)
    if nargin < 2 || isempty(bin_size)
        bin_size = 0.001;
    end
    if nargin < 3 || isempty(window_size)
        window_size = 0.02;
    end
    if nargin < 4 || isempty(normalize)
        normalize = false;
    end
    if nargin < 5 || isempty(show_full)
        show_full = true;
    end

    sampRate = Spike.P.sample_rate;
    binSize = bin_size * sampRate;
    window = window_size * sampRate;
    jitter = 0.005 * sampRate;
    nIter = 1000;

    bin = 0:binSize:window;
    nBin = length(bin) - 1;

    nClu = length(Spike.time);
    sp = floor(cell2mat(Spike.time) * sampRate);
    spidx = cell2mat(cellfun(@(x, y) y * ones(length(x), 1), Spike.time, num2cell((1:nClu)'), 'UniformOutput', false));
    nSp = length(sp);
    nSpClu = cellfun(@length, Spike.time);

    params = struct();
    params.nClu = nClu;
    params.nSp = nSp;
    params.nSpClu = nSpClu;
    params.bin = bin;
    params.binSize = bin_size;
    params.nBin = nBin;
    params.seed = now;
    params.normalize = normalize;
    params.show_full = show_full;
    rng(params.seed);

    if params.show_full
        ccgJitter = zeros(nIter, nClu, nClu, 2*nBin);
    else
        ccgJitter = zeros(nIter, nClu, nClu, nBin);
    end
    parfor iIter = 1:nIter
        spike = sortrows([sp + (rand(nSp, 1) * 2 - 1) * jitter, spidx]);
        temp = calc_ccg(spike, params);
        ccgJitter(iIter, :) = temp(:);
    end

    % calculate global range and pointwise range
    global_min = squeeze(quantile(min(ccgJitter, [], 4), 0.01));
    global_max = squeeze(quantile(max(ccgJitter, [], 4), 0.99));
    point = quantile(ccgJitter, [0.01, 0.99]);
    med = squeeze(median(ccgJitter, 1));
    m = squeeze(mean(ccgJitter, 1));
    s = squeeze(std(ccgJitter, 0, 1));

    spike = sortrows([sp, spidx]);
    ccg = calc_ccg(spike, params);

    z = (ccg - m) ./ s;

    % result
    out = struct();
    if params.show_full
        out.t = -window_size + bin_size/2 : bin_size : window_size - bin_size/2;
    else
        out.t = bin_size/2:bin_size:window_size - bin_size/2;
    end
    inT = out.t >= 0.001 & out.t <= 0.004;
    out.up = false(nClu, nClu);
    out.down = false(nClu, nClu);
    out.up(any(ccg(:, :, inT) > global_max, 3)) = true;
    out.down(any(ccg(:, :, inT) < global_min, 3)) = true;
    out.point = point;
    out.median = med;
    out.params = params;

    out.global_min = global_min;
    out.global_max = global_max;
    out.mean = m;
    out.std = s;
    out.zscore = z;
    out.median = med;

    out.ccg = ccg;
    %out.ccgJitter = ccgJitter;

    % Exclude where the peak at time 0 is too high
    [out_i, out_j] = find(squeeze(out.ccg(:, :, 1)) > 2 * out.global_max);
    out.up(out_i, out_j) = false;
    out.up(out_j, out_i) = false;
    out.down(out_i, out_j) = false;
    out.down(out_i, out_j) = false;

    % Remove autocorrelation
    out_idx = sub2ind([nClu, nClu], 1:nClu, 1:nClu);
    out.up(out_idx) = false;
    out.down(out_idx) = false;
    
end


function out = calc_ccg(spike, params)
    nClu = params.nClu;
    nBin = params.nBin;
    nSp = params.nSp;
    bin = params.bin;
    nSpClu = params.nSpClu;
    binSize = params.binSize;

    out = zeros(nClu, nClu, nBin);
    for shift = 1:nSp - 1
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

    if params.show_full
        other_half = flip(permute(out, [2, 1, 3]), 3);
        out = cat(3, other_half, out);
    end

    if params.normalize
        out = out ./ (nSpClu * binSize);
    end

end

