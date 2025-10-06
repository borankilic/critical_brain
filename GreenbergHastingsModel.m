classdef GreenbergHastingsModel < DynamicsModel
    % Greenberg-Hastings cellular automaton model for excitable dynamics
    % States: 0 = refractory, 1 = quiescent, 2 = excited
    
    properties
        r2              % probability of transitioning from refractory to quiescent
        threshold       % threshold for quiescent nodes to become excited
        r1      % External excitation probability for quiscent nodes to go excited
    end
    
    properties (Access = private)
        numExcited      % track number of excited nodes over time
        numRefractory   % track number of refractory nodes over time
        numQuiescent    % track number of quiescent nodes over time
    end
    
    methods
        function obj = GreenbergHastingsModel(net, params)
            % Constructor for Greenberg-Hastings model
            % Required parameters:
            %   - r2: probability of refractory → quiescent transition
            %   - threshold: excitation threshold for weighted neighbor sum
            % Optional parameters:
            %   - seed: random seed
            %   - initial_excited: fraction or indices of initially excited nodes
            
            obj@DynamicsModel(net, params);
            
            % Validate required parameters
            requiredFields = {'r2', 'threshold'};
            obj.validateParameters(requiredFields);
            
            % Set model-specific properties
            obj.r2 = params.r2;
            obj.threshold = params.threshold;
            obj.r1 = params.r1;
            
            % Validate parameter values
            if obj.r2 < 0 || obj.r2 > 1
                error('r2 must be between 0 and 1');
            end
            if obj.threshold < 0
                error('threshold must be non-negative');
            end
        end
        
        function initialize(obj)
            % Initialize the Greenberg-Hastings model
            
            % Reset random seed for reproducibility
            if ~isempty(obj.rngSeed)
                rng(obj.rngSeed);
            end
            
            % Initialize all nodes to quiescent state (state 2)
            obj.state = 2*ones(obj.net.N, 1);
            
            % Set initial excitation pattern
            if isfield(obj.params, 'initial_excited')
                initial = obj.params.initial_excited;
                
                if isscalar(initial) && initial <= 1
                    % Treat as fraction of nodes to excite randomly
                    numToExcite = round(initial * obj.net.N);
                    excitedNodes = randperm(obj.net.N, numToExcite);
                    obj.state(excitedNodes) = 1; % Set to excited state
                    
                elseif isvector(initial) && all(initial == round(initial))
                    % Treat as specific node indices
                    if any(initial < 1) || any(initial > obj.net.N)
                        error('Initial excited node indices out of range');
                    end
                    obj.state(initial) = 1; % Set to excited state
                    
                else
                    error('initial_excited must be a fraction (0-1) or vector of node indices');
                end
            else
                % Default: excite a small random fraction
                numToExcite = max(1, round(0.01 * obj.net.N));
                excitedNodes = randperm(obj.net.N, numToExcite);
                obj.state(excitedNodes) = 1; % Set to excited state
            end
            
            obj.isInitialized = true;
            obj.currentTime = 0;
            
            % Initialize tracking arrays
            obj.numExcited = [];
            obj.numRefractory = [];
            obj.numQuiescent = [];
        end
        
        function step(obj, t)
            % Perform one time step of Greenberg-Hastings dynamics
            % Rules:
            % 1. Refractory (0) → Quiescent (2) with probability r2
            % 2. Excited (1) → Refractory (0) always
            % 3. Quiescent (2) → Excited (1) if weighted sum of excited neighbors > threshold
            
            currentState = obj.state;
            newState = currentState;
            
            % Get weighted adjacency matrix
            A = obj.net.A; % This should be weighted
            
            % Find excited neighbors for each node
            excitedNodes = (currentState == 1);
            
            % Calculate weighted sum of excited neighbors for each node
            weightedExcitedNeighbors = A * double(excitedNodes);
            
            % Apply GH rules:
            
            % 1. Refractory nodes (state 0) transition to quiescent with probability r2
            refractoryNodes = (currentState == 0);
            transitionToQuiescent = (rand(obj.net.N, 1) < obj.r2);
            newState(refractoryNodes & transitionToQuiescent) = 2; % → quiescent
            newState(refractoryNodes & ~transitionToQuiescent) = 0; % stay refractory
            
            % 2. Excited nodes (state 2) always become refractory (state 0)
            newState(currentState == 1) = 0;
            
            % 3. Quiescent nodes (state 2) become excited if weighted sum exceeds threshold
            quiescentNodes = (currentState == 2);
            shouldExcite = (weightedExcitedNeighbors > obj.threshold) |  (rand(obj.net.N, 1) < obj.r1);
            newState(quiescentNodes & shouldExcite) = 1; % → excited
            newState(quiescentNodes & ~shouldExcite) = 2; % stay quiescent
            
            % Update state
            obj.state = newState;
            
            % Record statistics
            obj.recordOutputs();
        end
        
        function recordOutputs(obj)
            % Record current state and update statistics
           
            
            % Update state-specific counters
            obj.numExcited(end+1) = sum(obj.state == 2);
            obj.numRefractory(end+1) = sum(obj.state == 0);
            obj.numQuiescent(end+1) = sum(obj.state == 1);
        end
        
        function activity = getPopulationActivity(obj)
            % Get fraction of excited nodes over time
            if isempty(obj.stateHistory)
                error('No simulation data available. Run simulation first.');
            end
            
            % Count excited nodes (state == 1) at each time step
            excitedHistory = (obj.stateHistory == 1);
            activity = sum(excitedHistory, 1) / obj.net.N;
        end
        
        function refractoryActivity = getRefractoryActivity(obj)
            % Get fraction of refractory nodes over time
            if isempty(obj.stateHistory)
                error('No simulation data available. Run simulation first.');
            end
            
            % Count refractory nodes (state == 0) at each time step
            refractoryHistory = (obj.stateHistory == 0);
            refractoryActivity = sum(refractoryHistory, 1) / obj.net.N;
        end
        
        function quiescentActivity = getQuiescentActivity(obj)
            % Get fraction of quiescent nodes over time
            if isempty(obj.stateHistory)
                error('No simulation data available. Run simulation first.');
            end
            
            % Count quiescent nodes (state == 2) at each time step
            quiescentHistory = (obj.stateHistory == 2);
            quiescentActivity = sum(quiescentHistory, 1) / obj.net.N;
        end
        
        function plotStateDistribution(obj)
            % Plot distribution of node states over time
            if isempty(obj.stateHistory)
                error('No simulation data available. Run simulation first.');
            end
            
            % Count nodes in each state at each time step
            maxState = max(obj.stateHistory(:));
            stateCounts = zeros(maxState + 1, size(obj.stateHistory, 2));
            
            for state = 0:maxState
                stateCounts(state + 1, :) = sum(obj.stateHistory == state, 1);
            end
            
            % Convert to fractions
            stateFractions = stateCounts / obj.net.N;
            
            figure;
            timeSteps = 1:size(obj.stateHistory, 2);
            
            % Create stacked area plot
            area(timeSteps, stateFractions');
            
            % Create legend
            legendLabels = cell(maxState + 1, 1);
            legendLabels{1} = 'Quiescent';
            legendLabels{2} = 'Excited';
            legendLabels{3} = sprintf('Refractory');

            
            legend(legendLabels, 'Location', 'best');
            xlabel('Time Step');
            ylabel('Fraction of Nodes');
            title('Distribution of Node States Over Time');
            ylim([0, 1]);
            grid on;
        end

        function plotStateEvolution(obj)
            % Plot evolution of all three states over time
            if isempty(obj.stateHistory)
                error('No simulation data available. Run simulation first.');
            end
            
            figure;
            timeSteps = 1:size(obj.stateHistory, 2);
            
            % Calculate fractions
            excitedFrac = obj.getPopulationActivity();
            refractoryFrac = obj.getRefractoryActivity();
            quiescentFrac = obj.getQuiescentActivity();
            
            % Plot all three states
            plot(timeSteps, excitedFrac, 'r-', 'LineWidth', 2, 'DisplayName', 'Excited');
            hold on;
            plot(timeSteps, refractoryFrac, 'b-', 'LineWidth', 2, 'DisplayName', 'Refractory');
            plot(timeSteps, quiescentFrac, 'g-', 'LineWidth', 2, 'DisplayName', 'Quiescent');
            
            xlabel('Time Step');
            ylabel('Fraction of Nodes');
            title('State Evolution in Greenberg-Hastings Model');
            legend('show');
            grid on;
            ylim([0, 1]);
        end
        
        function avalanches = detectAvalanches(obj, varargin)
            % Detect avalanche events in the dynamics
            % An avalanche is a sequence of continuous activity above baseline
            
            p = inputParser;
            addParameter(p, 'threshold', 0.01, @isnumeric);  % activity threshold
            addParameter(p, 'minDuration', 2, @isnumeric);   % minimum avalanche duration
            parse(p, varargin{:});
            
            threshold = p.Results.threshold;
            minDuration = p.Results.minDuration;
            
            if isempty(obj.stateHistory)
                error('No simulation data available. Run simulation first.');
            end
            
            activity = obj.getPopulationActivity();
            
            % Find periods above threshold
            aboveThreshold = activity > threshold;
            
            % Detect avalanche boundaries
            avalanches = [];
            inAvalanche = false;
            avalancheStart = 1;
            
            for t = 1:length(aboveThreshold)
                if aboveThreshold(t) && ~inAvalanche
                    % Start of avalanche
                    inAvalanche = true;
                    avalancheStart = t;
                elseif ~aboveThreshold(t) && inAvalanche
                    % End of avalanche
                    inAvalanche = false;
                    duration = t - avalancheStart;
                    
                    if duration >= minDuration
                        % Calculate avalanche properties
                        avalancheActivity = activity(avalancheStart:t-1);
                        avalanche = struct();
                        avalanche.start = avalancheStart;
                        avalanche.end = t - 1;
                        avalanche.duration = duration;
                        avalanche.size = sum(avalancheActivity);
                        avalanche.peakActivity = max(avalancheActivity);
                        
                        avalanches = [avalanches, avalanche]; %#ok<AGROW>
                    end
                end
            end
            
            % Handle case where simulation ends during avalanche
            if inAvalanche
                duration = length(aboveThreshold) - avalancheStart + 1;
                if duration >= minDuration
                    avalancheActivity = activity(avalancheStart:end);
                    avalanche = struct();
                    avalanche.start = avalancheStart;
                    avalanche.end = length(aboveThreshold);
                    avalanche.duration = duration;
                    avalanche.size = sum(avalancheActivity);
                    avalanche.peakActivity = max(avalancheActivity);
                    
                    avalanches = [avalanches, avalanche];
                end
            end
        end

        function summary = getSummaryStats(obj)
            % Get summary statistics specific to Greenberg-Hastings model
            
            % Get base summary from parent class
            summary = getSummaryStats@DynamicsModel(obj);
            
            % Add GH-specific statistics
            summary.model = 'Greenberg-Hastings';
            summary.parameters.r2 = obj.r2;
            summary.parameters.threshold = obj.threshold;
            
            if ~isempty(obj.stateHistory)
                % State-specific statistics
                excitedFraction = sum(obj.stateHistory(:) == 2) / numel(obj.stateHistory);
                refractoryFraction = sum(obj.stateHistory(:) == 0) / numel(obj.stateHistory);
                quiescentFraction = sum(obj.stateHistory(:) == 1) / numel(obj.stateHistory);
                
                summary.stateFractions.excited = excitedFraction;
                summary.stateFractions.refractory = refractoryFraction;
                summary.stateFractions.quiescent = quiescentFraction;
                
                % Calculate mean activity levels
                summary.meanActivity.excited = mean(obj.getPopulationActivity());
                summary.meanActivity.refractory = mean(obj.getRefractoryActivity());
                summary.meanActivity.quiescent = mean(obj.getQuiescentActivity());
            end
        end
        
        function displaySummary(obj)
            % Display summary specific to Greenberg-Hastings model
            
            summary = obj.getSummaryStats();
            
            fprintf('\n=== Greenberg-Hastings Model Summary ===\n');
            fprintf('Network size: %d nodes\n', summary.networkSize);
            fprintf('Simulation length: %d time steps\n', summary.simulationLength);
            statusOptions = {'Incomplete', 'Completed'};
            fprintf('Status: %s\n', statusOptions{summary.isCompleted + 1});
            
            fprintf('\nModel Parameters:\n');
            fprintf('  r2 (refractory→quiescent probability): %.4f\n', summary.parameters.r2);
            fprintf('  Threshold for excitation: %.4f\n', summary.parameters.threshold);
            fprintf('  Random seed: %d\n', summary.rngSeed);
            
            if isfield(summary, 'stateFractions')
                fprintf('\nState Occupancy Fractions:\n');
                fprintf('  Quiescent: %.3f\n', summary.stateFractions.quiescent);
                fprintf('  Excited: %.3f\n', summary.stateFractions.excited);
                fprintf('  Refractory: %.3f\n', summary.stateFractions.refractory);
                
                fprintf('\nMean Activity Levels:\n');
                fprintf('  Quiescent: %.3f\n', summary.meanActivity.quiescent);
                fprintf('  Excited: %.3f\n', summary.meanActivity.excited);
                fprintf('  Refractory: %.3f\n', summary.meanActivity.refractory);
            end
            
            fprintf('========================================\n\n');
        end
    end

   
    
    methods (Static)
        function params = getDefaultParams()
            % Get default parameters for Greenberg-Hastings model
            params = struct();
            params.r2 = 0.01;
            params.threshold = 1.5;
            params.r1 = 0.001;
            params.initial_excited = 0.1;  % 1% initially excited
            params.seed = 12345;            % reproducible results
        end
        
        function model = createWithRandomNetwork(N, connectivity, varargin)
            % Create GH model with a random network
            % N: number of nodes
            % connectivity: connection probability or mean degree
            
            p = inputParser;
            addParameter(p, 'params', GreenbergHastingsModel.getDefaultParams(), @isstruct);
            addParameter(p, 'networkType', 'erdos_renyi', @ischar);
            parse(p, varargin{:});
            
            params = p.Results.params;
            networkType = p.Results.networkType;
            
            % Create network
            switch lower(networkType)
                case 'erdos_renyi'
                    if connectivity < 1
                        % Treat as probability
                        net = Network.generateRandomNetwork(N, connectivity);
                    else
                        % Treat as mean degree
                        p_conn = connectivity / (N - 1);
                        net = Network.generateRandomNetwork(N, p_conn);
                    end
                case 'small_world'
                    % For small world, connectivity is mean degree
                    beta = 0.3;  % rewiring probability
                    net = Network.generateSmallWorldNetwork(N, round(connectivity), beta);
                otherwise
                    error('Unsupported network type: %s', networkType);
            end
            
            % Create model
            model = GreenbergHastingsModel(net, params);
        end
    end
end