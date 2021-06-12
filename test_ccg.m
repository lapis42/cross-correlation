%% Data structure
Spike = struct();
Spike.nUnit = 8;
Spike.time = cell(Spike.nUnit, 1);
Spike.P.sample_rate = 30000;

%% generate data (8 units, 1 hour)
rng(1);
for iU = 1:Spike.nUnit
    fr = rand() * 20;
    iti = -log(rand(20*60*60, 1)) / fr;
    iti = iti + rand(20*60*60, 1) * 0.002 + 0.001; % add refractory period
    sps = cumsum(iti);
    sps(sps>60*60) = []; % remove spike after 60 min
    Spike.time{iU} = sps;
end

%% inject correlated spikes from unit 1 to unit 2 (20%, 5+-1 ms)
sp = Spike.time{1};
nSp = length(sp);
nrSp = floor(nSp / 20);
random_sp = sp(randi(nSp, nrSp, 1)) + 0.001 * randn(nrSp, 1) + 0.005;
Spike.time{2} = sort([Spike.time{2}; random_sp]);


%% main
tic;
[out, t] = ccg(Spike, 0.001, 0.020, true);
toc;


%% plot
fig = figure(123);
for iU = 1:Spike.nUnit
    for jU = 1:Spike.nUnit
        subplot(8, 8, Spike.nUnit * (iU - 1) + jU);
        bar(t, squeeze(out(iU, jU, :)), 1, 'LineStyle', 'none', 'FaceColor', 'k');
    end
end
