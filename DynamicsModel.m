classdef (Abstract) DynamicsModel < handle
    % Abstract base class for network dynamics models
    % Follows composition pattern - takes Network instance as dependency
    
    properties
        net           % instance of Network class
        T             % number of time steps
        state         % current state vector (N x 1)
        stateHistory  % state history matrix (N x T)
        params        % struct containing model parameters
        rngSeed       % random number generator seed for reproducibility
        currentTime   % current simulation time step
        dt            % time step size (for continuous models)
    end
    
    properties (SetAccess = protected)
        isInitialized = false  % flag to track initialization
        isCompleted = false    % flag to track if simulation completed
    end
    
    methods
        function obj = DynamicsModel(net, params)
            % Constructor for dynamics model
            % net: Network instance
            % params: struct with model parameters
            
            if ~isa(net, 'Network')
                error('First argument must be a Network instance');
            end
            
            obj.net = net;
            obj.params = params;
            obj.currentTime = 0;
            
            % Set random seed if provided
            if isfield(params, 'seed')
                obj.rngSeed = params.seed;
                rng(obj.rngSeed);
            else
                obj.rngSeed = rng().Seed; % store current seed
            end
            
            % Set time step for continuous models
            if isfield(params, 'dt')
                obj.dt = params.dt;
            else
                obj.dt = 0.01; % default time step
            end
        end
        
        function run(obj, T, varargin)
            % Run the dynamics simulation
            % T: number of time steps to simulate
            % Optional parameters: 'recordEvery', 'verbose'
            
            p = inputParser;
            addParameter(p, 'recordEvery', 1, @isnumeric);
            addParameter(p, 'verbose', false, @islogical);
            parse(p, varargin{:});
            
            recordEvery = p.Results.recordEvery;
            verbose = p.Results.verbose;
            
            obj.T = T;
            
            % Initialize if not already done
            if ~obj.isInitialized
                obj.initialize();
            end
            
            % Pre-allocate history matrix (record only every recordEvery steps)
            recordSteps = floor(T / recordEvery);
            obj.stateHistory = zeros(obj.net.N, recordSteps);
            
            if verbose
                fprintf('Running dynamics simulation for %d time steps...\n', T);
                if T > 1000
                    fprintf('Progress: ');
                end
            end
            
            recordIndex = 1;
            
            % Main simulation loop
            for t = 1:T
                obj.currentTime = t;
                obj.step(t);
                
                % Record state if needed
                if mod(t, recordEvery) == 0
                    obj.stateHistory(:, recordIndex) = obj.getState();
                    recordIndex = recordIndex + 1;
                end
                
                % Progress reporting for long simulations
                if verbose && T > 1000 && mod(t, floor(T/10)) == 0
                    fprintf('%d%% ', round(100*t/T));
                end
            end
            
            obj.isCompleted = true;
            
            if verbose
                if T > 1000, fprintf('\n'); end
                fprintf('Simulation completed.\n');
            end
        end
        
        function reset(obj)
            % Reset the dynamics model to initial state
            obj.isInitialized = false;
            obj.isCompleted = false;
            obj.currentTime = 0;
            obj.state = [];
            obj.stateHistory = [];
            obj.T = [];
            
            % Reset random seed
            if ~isempty(obj.rngSeed)
                rng(obj.rngSeed);
            end
        end
        
        function state = getState(obj)
            % Get current state vector
            state = obj.state;
        end
        
        function setState(obj, newState)
            % Set current state vector
            if length(newState) ~= obj.net.N
                error('State vector length must match number of nodes (%d)', obj.net.N);
            end
            obj.state = newState(:); % ensure column vector
        end
        
        function history = getStateHistory(obj)
            % Get complete state history
            if obj.isCompleted
                history = obj.stateHistory;
            else
                warning('Simulation not completed. Returning partial history.');
                history = obj.stateHistory;
            end
        end
        
        function activity = getPopulationActivity(obj)
            % Get population-level activity over time
            if isempty(obj.stateHistory)
                error('No simulation data available. Run simulation first.');
            end
            
            % This is a generic implementation - subclasses may override
            activity = mean(obj.stateHistory, 1);
        end
        
        function plotRaster(obj, varargin)
            % Plot raster plot of network activity
            if isempty(obj.stateHistory)
                error('No simulation data available. Run simulation first.');
            end
            
            p = inputParser;
            addParameter(p, 'threshold', 0.5, @isnumeric);
            addParameter(p, 'maxNodes', 100, @isnumeric);
            parse(p, varargin{:});
            
            threshold = p.Results.threshold;
            maxNodes = p.Results.maxNodes;
            
            % Limit number of nodes for visualization
            if obj.net.N > maxNodes
                nodeIndices = sort(randperm(obj.net.N, maxNodes));
                data = obj.stateHistory(nodeIndices, :);
                fprintf('Showing %d randomly selected nodes out of %d\n', maxNodes, obj.net.N);
            else
                data = obj.stateHistory;
                nodeIndices = 1:obj.net.N;
            end
            
            % Create raster plot
            figure;
            [nodeIdx, timeIdx] = find(data > threshold);
            scatter(timeIdx, nodeIndices(nodeIdx), 1, 'k', 'filled');
            
            xlabel('Time Step');
            ylabel('Node Index');
            title('Network Activity Raster Plot');
            ylim([0.5, length(nodeIndices) + 0.5]);
            grid on;
        end
        
        function plotPopulationActivity(obj)
            % Plot population-level activity over time
            if isempty(obj.stateHistory)
                error('No simulation data available. Run simulation first.');
            end
            
            activity = obj.getPopulationActivity();
            
            figure;
            plot(1:length(activity), activity, 'LineWidth', 2);
            xlabel('Time Step');
            ylabel('Population Activity');
            title('Population Activity Over Time');
            grid on;
        end
        
        
        
      %% Functions to animate network graphics        
        
      function animateNetworkActivity(obj, varargin)
        % VISUALIZE_NETWORK_DYNAMICS Creates an interactive visualization of network node states
        %
        % Inputs:
        %   positions - N x 3 matrix of node positions [x, y, z]
        %   states    - N x M matrix of node states (0 or 1) over M time steps
        %   varargin  - optional parameter-value pairs:
        %               'FrameRate', value - frames per second (default: 5)
        %               'NodeSize', value  - size of nodes (default: 50)
        %               'AutoPlay', true/false - start animation automatically (default: true)
        %               'Loop', true/false - loop animation (default: true)
        %
        % Example:
        %   visualize_network_dynamics(pos, states, 'FrameRate', 10, 'NodeSize', 100)
        %
        % Controls:
        %   - Play/Pause button to control animation
        %   - Slider to manually navigate through time steps
        %   - Speed control slider
        %   - Reset button to return to first frame
        states = obj.stateHistory;
        positions = obj.net.coords;
        % Parse optional arguments
        p = inputParser;
        addParameter(p, 'FrameRate', 10, @isnumeric);
        addParameter(p, 'NodeSize', 100, @isnumeric);
        addParameter(p, 'AutoPlay', true, @islogical);
        addParameter(p, 'Loop', true, @islogical);
        parse(p, varargin{:});
        
        frame_rate = p.Results.FrameRate;
        node_size = p.Results.NodeSize;
        auto_play = p.Results.AutoPlay;
        loop_animation = p.Results.Loop;
        
        % Get dimensions
        [N, ~] = size(positions);
        [~, M] = size(states);
        

        % Create figure with controls
        fig = figure('Position', [100, 100, 900, 700], 'Name', 'Network Dynamics Animation', ...
                     'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'figure');
        set(fig, 'Color', 'white');
        
        % Create main plot area
        main_panel = uipanel('Parent', fig, 'Position', [0, 0.15, 1, 0.85], ...
                            'BackgroundColor', 'white', 'BorderType', 'none');
        
        % Create control panel
        control_panel = uipanel('Parent', fig, 'Position', [0, 0, 1, 0.15], ...
                               'BackgroundColor', [0.94, 0.94, 0.94], 'BorderType', 'line');
        
        % Set up axes in main panel
        ax = axes('Parent', main_panel, 'Position', [0.1, 0.1, 0.8, 0.8]);
        

        % 2D plot
        axis equal;
        grid on;
        xlabel('X Position');
        ylabel('Y Position');
        title('Network Node Dynamics (2D)');
        
        % Set axis limits with some padding
        x_range = [min(positions(:,1)), max(positions(:,1))];
        y_range = [min(positions(:,2)), max(positions(:,2))];
        x_padding = 0.1 * (x_range(2) - x_range(1));
        y_padding = 0.1 * (y_range(2) - y_range(1));
        if x_padding == 0, x_padding = 1; end
        if y_padding == 0, y_padding = 1; end
        xlim([x_range(1) - x_padding, x_range(2) + x_padding]);
        ylim([y_range(1) - y_padding, y_range(2) + y_padding]);

        
        % Initialize scatter plot handle

        scatter_handle = scatter(ax, positions(:,1), positions(:,2), node_size, ...
                                   'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
        
        % Add time step counter
        time_text = text(ax, 0.02, 0.95, '', 'Units', 'normalized', 'FontSize', 14, ...
                         'FontWeight', 'bold', 'BackgroundColor', 'white', ...
                         'EdgeColor', 'black', 'Margin', 5);
        
        % Create UI controls
        button_height = 0.4;
        button_y = 0.3;
        
        % Play/Pause button
        play_button = uicontrol('Parent', control_panel, 'Style', 'pushbutton', ...
                               'String', 'Pause', 'Position', [20, button_y*50, 80, button_height*50], ...
                               'FontSize', 10, 'FontWeight', 'bold');
        
        % Reset button
        reset_button = uicontrol('Parent', control_panel, 'Style', 'pushbutton', ...
                                'String', 'Reset', 'Position', [110, button_y*50, 60, button_height*50], ...
                                'FontSize', 10);
        
        % Time slider
        uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'Time Step:', ...
                  'Position', [180, button_y*50+20, 80, 20], 'FontSize', 10, ...
                  'BackgroundColor', [0.94, 0.94, 0.94]);
        
        time_slider = uicontrol('Parent', control_panel, 'Style', 'slider', ...
                               'Min', 1, 'Max', M, 'Value', 1, ...
                               'Position', [270, button_y*50, 300, button_height*50], ...
                               'SliderStep', [1/(M-1), 10/(M-1)]);
        
        % Speed control
        uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'Speed (fps):', ...
                  'Position', [580, button_y*50+20, 80, 20], 'FontSize', 10, ...
                  'BackgroundColor', [0.94, 0.94, 0.94]);
        
        speed_slider = uicontrol('Parent', control_panel, 'Style', 'slider', ...
                                'Min', 1, 'Max', 20, 'Value', frame_rate, ...
                                'Position', [670, button_y*50, 150, button_height*50], ...
                                'SliderStep', [1/19, 5/19]);
        
        % Speed display
        speed_text = uicontrol('Parent', control_panel, 'Style', 'text', ...
                              'String', sprintf('%.1f', frame_rate), ...
                              'Position', [830, button_y*50, 40, button_height*50], ...
                              'FontSize', 10, 'BackgroundColor', 'white');
        
        % Animation state variables
        current_frame = 1;
        is_playing = auto_play;
        animation_timer = [];
        
        % Function to update display
        function update_display(frame_num)
            if frame_num < 1, frame_num = 1; end
            if frame_num > M, frame_num = M; end
            
            current_frame = frame_num;
            
            % Update node colors based on states
            colors = zeros(N, 3);
            inactive_nodes = states(:, current_frame) == 0;
            quiscent_nodes = states(:, current_frame) == 1;
            excited_nodes = states(:, current_frame) == 2;
            colors(inactive_nodes, :) = repmat([1, 1, 1], sum(inactive_nodes), 1);
            colors(excited_nodes, :) = repmat([1, 0, 0], sum(excited_nodes), 1);
            
            % Update scatter plot colors
            set(scatter_handle, 'CData', colors);
            
            % Update time counter
            set(time_text, 'String', sprintf('Time Step: %d/%d', current_frame, M));
            
            % Update slider without triggering callback
            set(time_slider, 'Value', current_frame);
            
            drawnow;
        end
        
        % Function to start/stop animation
        function toggle_animation()
            if is_playing
                % Stop animation
                if ~isempty(animation_timer)
                    stop(animation_timer);
                    delete(animation_timer);
                    animation_timer = [];
                end
                is_playing = false;
                set(play_button, 'String', 'Play');
            else
                % Start animation
                current_speed = get(speed_slider, 'Value');
                animation_timer = timer('ExecutionMode', 'fixedRate', ...
                                       'Period', 1/current_speed, ...
                                       'TimerFcn', @animate_step);
                start(animation_timer);
                is_playing = true;
                set(play_button, 'String', 'Pause');
            end
        end
        
        % Animation step function
        function animate_step(~, ~)
            if current_frame >= M
                if loop_animation
                    current_frame = 1;
                else
                    toggle_animation();
                    return;
                end
            else
                current_frame = current_frame + 1;
            end
            update_display(current_frame);
        end
        
        % Callback functions
        set(play_button, 'Callback', @(~,~) toggle_animation());
        set(reset_button, 'Callback', @(~,~) update_display(1));
        set(time_slider, 'Callback', @(src,~) update_display(round(get(src, 'Value'))));
        set(speed_slider, 'Callback', @(src,~) update_speed_callback(src));
        
        function update_speed_callback(src)
            new_speed = get(src, 'Value');
            set(speed_text, 'String', sprintf('%.1f', new_speed));
            
            % If animation is playing, restart timer with new speed
            if is_playing
                if ~isempty(animation_timer)
                    stop(animation_timer);
                    delete(animation_timer);
                end
                animation_timer = timer('ExecutionMode', 'fixedRate', ...
                                       'Period', 1/new_speed, ...
                                       'TimerFcn', @animate_step);
                start(animation_timer);
            end
        end
        
        % Clean up timer when figure is closed
        set(fig, 'CloseRequestFcn', @cleanup_and_close);
        
        function cleanup_and_close(~, ~)
            if ~isempty(animation_timer)
                stop(animation_timer);
                delete(animation_timer);
            end
            delete(fig);
        end
        
        % Initialize display
        update_display(1);
        
        % Start animation if auto_play is enabled
        if auto_play
            toggle_animation();
        end
        
        fprintf('Network dynamics visualization created.\n');
        fprintf('Controls:\n');
        fprintf('  - Play/Pause: Start/stop animation\n');
        fprintf('  - Reset: Return to first frame\n');
        fprintf('  - Time Slider: Navigate manually through time steps\n');
        fprintf('  - Speed Slider: Adjust animation speed (1-20 fps)\n');
        end



        %% Result and SUmmary
        function saveResults(obj, filepath, varargin)
            % Save simulation results to file
            
            p = inputParser;
            addParameter(p, 'saveNetwork', true, @islogical);
            addParameter(p, 'saveParameters', true, @islogical);
            parse(p, varargin{:});
            
            % Prepare data structure
            results = struct();
            results.stateHistory = obj.stateHistory;
            results.T = obj.T;
            results.currentTime = obj.currentTime;
            results.rngSeed = obj.rngSeed;
            results.isCompleted = obj.isCompleted;
            results.simulationDate = datestr(now);
            
            if p.Results.saveParameters
                results.params = obj.params;
                results.dt = obj.dt;
            end
            
            if p.Results.saveNetwork
                results.network = struct();
                results.network.A = obj.net.A;
                results.network.N = obj.net.N;
                results.network.directed = obj.net.directed;
                results.network.weighted = obj.net.weighted;
            end
            
            save(filepath, 'results');
            fprintf('Results saved to %s\n', filepath);
        end
        
        function summary = getSummaryStats(obj)
            % Get summary statistics of the simulation
            if isempty(obj.stateHistory)
                error('No simulation data available. Run simulation first.');
            end
            
            summary = struct();
            summary.simulationLength = obj.T;
            summary.networkSize = obj.net.N;
            summary.meanActivity = mean(obj.stateHistory(:));
            summary.stdActivity = std(obj.stateHistory(:));
            summary.maxActivity = max(obj.stateHistory(:));
            summary.minActivity = min(obj.stateHistory(:));
            
            % Population activity statistics
            popActivity = obj.getPopulationActivity();
            summary.meanPopActivity = mean(popActivity);
            summary.stdPopActivity = std(popActivity);
            summary.maxPopActivity = max(popActivity);
            summary.minPopActivity = min(popActivity);
            
            summary.rngSeed = obj.rngSeed;
            summary.isCompleted = obj.isCompleted;
        end
        
        function displaySummary(obj)
            % Display summary of simulation results
            if isempty(obj.stateHistory)
                fprintf('No simulation data available.\n');
                return;
            end
            
            summary = obj.getSummaryStats();
            
            fprintf('\n=== Dynamics Simulation Summary ===\n');
            fprintf('Network size: %d nodes\n', summary.networkSize);
            fprintf('Simulation length: %d time steps\n', summary.simulationLength);
            fprintf('Status: %s\n', char(summary.isCompleted * "Completed" + ~summary.isCompleted * "Incomplete"));
            fprintf('Random seed: %d\n', summary.rngSeed);
            fprintf('\nActivity Statistics:\n');
            fprintf('  Mean: %.4f (±%.4f)\n', summary.meanActivity, summary.stdActivity);
            fprintf('  Range: [%.4f, %.4f]\n', summary.minActivity, summary.maxActivity);
            fprintf('\nPopulation Activity:\n');
            fprintf('  Mean: %.4f (±%.4f)\n', summary.meanPopActivity, summary.stdPopActivity);
            fprintf('  Range: [%.4f, %.4f]\n', summary.minPopActivity, summary.maxPopActivity);
            fprintf('=====================================\n\n');
        end
    end
    
    
    methods (Access = protected)
        function validateParameters(obj, requiredFields)
            % Helper method to validate that required parameters are present
            for i = 1:length(requiredFields)
                if ~isfield(obj.params, requiredFields{i})
                    error('Required parameter "%s" not found in params struct', requiredFields{i});
                end
            end
        end
        
        function recordOutputs(obj, t)
            % Optional method for recording additional outputs
            % Can be overridden by subclasses
            % Default implementation does nothing
        end
        
        function layout = generateNetworkLayout(obj)
            % Generate 2D layout for network visualization using force-directed algorithm
            % This is a simplified spring-embedded layout
            
            N = obj.net.N;
            A = obj.net.A;
            
            % Initialize random positions
            layout = randn(N, 2) * 0.5;
            
            % Force-directed layout parameters
            numIterations = 100;
            dt = 0.01;
            k_spring = 0.1;      % spring constant for connected nodes
            k_repulsion = 0.05;  % repulsion constant for all node pairs
            damping = 0.9;       % velocity damping
            
            velocities = zeros(N, 2);
            
            % Precompute distance matrix for connected nodes
            [connectedI, connectedJ] = find(A);
            
            for iter = 1:numIterations
                forces = zeros(N, 2);
                
                % Repulsive forces between all pairs
                for i = 1:N
                    for j = i+1:N
                        diff = layout(i, :) - layout(j, :);
                        dist = norm(diff) + 1e-6;  % avoid division by zero
                        
                        if dist < 2.0  % only apply repulsion for close nodes
                            force = k_repulsion * diff / (dist^2);
                            forces(i, :) = forces(i, :) + force;
                            forces(j, :) = forces(j, :) - force;
                        end
                    end
                end
                
                % Attractive forces for connected nodes
                for edgeIdx = 1:length(connectedI)
                    i = connectedI(edgeIdx);
                    j = connectedJ(edgeIdx);
                    
                    diff = layout(j, :) - layout(i, :);
                    dist = norm(diff) + 1e-6;
                    
                    idealDist = 1.0;  % desired edge length
                    force = k_spring * (dist - idealDist) * diff / dist;
                    
                    forces(i, :) = forces(i, :) + force;
                    forces(j, :) = forces(j, :) - force;
                end
                
                % Update velocities and positions
                velocities = damping * velocities + dt * forces;
                layout = layout + dt * velocities;
                
                % Cooling schedule (reduce forces over time)
                if mod(iter, 20) == 0
                    dt = dt * 0.95;
                end
            end
        end
        
        function drawNetworkEdges(obj, nodePos, edgeColor, lineWidth)
            % Helper function to draw network edges
            A = obj.net.A;
            
            % Find edge indices
            if obj.net.directed
                [sourceNodes, targetNodes] = find(A);
            else
                [sourceNodes, targetNodes] = find(triu(A, 1));  % upper triangular to avoid duplicates
            end
            
            % Draw edges
            for edgeIdx = 1:length(sourceNodes)
                i = sourceNodes(edgeIdx);
                j = targetNodes(edgeIdx);
                
                plot([nodePos(i, 1), nodePos(j, 1)], [nodePos(i, 2), nodePos(j, 2)], ...
                     'Color', edgeColor, 'LineWidth', lineWidth);
            end
        end
        
        function drawTimeProgressBar(obj, currentTime, totalTime, xRange, yPos)
            % Helper function to draw time progress bar
            progress = currentTime / totalTime;
            barWidth = diff(xRange);
            
            % Background bar
            rectangle('Position', [xRange(1), yPos, barWidth, 0.02], ...
                     'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'k', 'LineWidth', 0.5);
            
            % Progress bar
            if progress > 0
                rectangle('Position', [xRange(1), yPos, barWidth * progress, 0.02], ...
                         'FaceColor', [0.2 0.6 0.2], 'EdgeColor', 'none');
            end
            
            % Time text
            text(xRange(1) + barWidth/2, yPos + 0.01, sprintf('t = %d / %d', currentTime, totalTime), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 8, 'FontWeight', 'bold');
        end
        
        function stateName = getStateName(obj, stateValue)
            % Get descriptive name for state value (can be overridden by subclasses)
            stateName = sprintf('State %g', stateValue);
        end
    end
end