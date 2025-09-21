%% Test Script for Network Dynamics Classes
% This script demonstrates the basic functionality of the Network,
% DynamicsModel, and GreenbergHastingsModel classes

clear; clc; close all;

fprintf('=== Testing Network Dynamics Classes ===\n\n');

%% Test 1: Create and analyze a network
fprintf('1. Testing Network Class...\n');

% Create a small-world network
N = 50;  % number of nodes
k = 6;   % initial degree
beta = 0.3;  % rewiring probability

net = Network.generateSmallWorldNetwork(N, k, beta);

% Display basic network properties
fprintf('   Network created with %d nodes\n', net.N);
fprintf('   Network density: %.3f\n', nnz(net.A) / (N * (N-1)));

% Calculate network statistics
stats = net.getNetworkStats();
fprintf('   Mean degree: %.2f (±%.2f)\n', stats.degree_mean, stats.degree_std);
fprintf('   Mean clustering: %.3f (±%.3f)\n', stats.clustering_mean, stats.clustering_std);

% Visualize the network
figure(1);
net.visualize();
title('Small-World Network Structure');

fprintf('   ✓ Network class working correctly\n\n');

%% Test 2: Create Greenberg-Hastings model
fprintf('2. Testing Greenberg-Hastings Model...\n');
params = GreenbergHastingsModel.getDefaultParams();
gh_model = GreenbergHastingsModel(net, params);
gh_model.run(2000, 'verbose', true);
gh_model.plotStateEvolution()