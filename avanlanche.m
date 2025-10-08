r1 = 0.001;
r2 = 0.64;
threshold = 0.17;

T= 20000;
% 25x25 grid test
n = 800;
node_number = n ;  

% Create grid connectivity (each node connected to 4 neighbors: up, down, left, right)
addpath('NCC_toolbox/');
load("weight_matrix.mat");
weight_matrix = weightMatrix;

connectivity = (weight_matrix >0);
net = Network(weight_matrix);
net.normalize_weights;

params = GreenbergHastingsModel.getDefaultParams();
params.r2 = r2;
params.r1=r1;
params.threshold=threshold;
params.initial_excited = 0;
ghModel = GreenbergHastingsModel(net, params);
ghModel.run(T)
excitedActivity = (ghModel.stateHistory==1);
sumColumn = excitedActivity' * ones(size(excitedActivity, 1),1);
aboveBaseline = (sumColumn>60);
excitedActivity = excitedActivity .* repmat(aboveBaseline',size(excitedActivity,1),1);
asdf2 = rastertoasdf2(excitedActivity,0.01, 'deneme', 'spikes', '11.08.2025');
% Compute all avalanche properties
Av = avprops(asdf2, 'ratio', 'fingerprint');

%% Plot histogram of branching ratios

minBR = min(Av.branchingRatio);
maxBR = max(Av.branchingRatio);

nEdges = 25;

edges = minBR:((maxBR-minBR)/(nEdges - 1)):maxBR;

BRhist = histc(Av.branchingRatio, edges);

figure;
plot(edges, BRhist);

title('Histogram of Avalanche Branching Ratios', 'fontsize', 14)
xlabel('Branching Ratio (s_{t+1} / s_t)', 'fontsize', 14)
ylabel('Frequency of Occurrence', 'fontsize', 14)

%% Compute power-law parameters using macro

% size distribution (SZ)
[tau, xminSZ, xmaxSZ, sigmaSZ, pSZ, pCritSZ, ksDR, DataSZ] =...
    avpropvals(Av.size, 'size', 'plot');

% size distribution (SZ) with cutoffs
UniqSizes = unique(Av.size);
Occurances = hist(Av.size,UniqSizes);
AllowedSizes = UniqSizes(Occurances >= 20);
AllowedSizes(AllowedSizes < 4) = [];
LimSize = Av.size(ismember(Av.size,AllowedSizes));
[tau, xminSZ, xmaxSZ, sigmaSZ, pSZ, pCritSZ, DataSZ] =...
    avpropvals(LimSize, 'size', 'plot');

% duration distribution (DR)
[alpha, xminDR, xmaxDR, sigmaDR, pDR, pCritDR, ksDR, DataDR] =...
    avpropvals(Av.duration, 'duration', 'plot');

% size given duration distribution (SD)
[sigmaNuZInvSD, waste, waste, sigmaSD] = avpropvals({Av.size, Av.duration},...
    'sizgivdur', 'durmin', xminDR{1}, 'durmax', xmaxDR{1}, 'plot');

%% Perform avalanche shape collapse for all shapes

% compute average temporal profiles
avgProfiles = avgshapes(Av.shape, Av.duration, 'cutoffs', 4, 20);

% plot all profiles
figure;
for iProfile = 1:length(avgProfiles)
    hold on
    plot(1:length(avgProfiles{iProfile}), avgProfiles{iProfile});
end
hold off

xlabel('Time Bin, t', 'fontsize', 14)
ylabel('Neurons Active, s(t)', 'fontsize', 14)
title('Mean Temporal Profiles', 'fontsize', 14)

% compute shape collapse statistics (SC) and plot
[sigmaNuZInvSC, secondDrv, range, errors] = avshapecollapse(avgProfiles, 'plot');

sigmaSC = avshapecollapsestd(avgProfiles);

title(['Avalanche Shape Collapse', char(10), '1/(sigma nu z) = ',...
    num2str(sigmaNuZInvSC), ' +/- ', num2str(sigmaSC)], 'fontsize', 14)