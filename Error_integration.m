% Updated parameters for the modified HH model
E_Na = 115;     % Sodium reversal potential (mV)
E_K = -12;      % Potassium reversal potential (mV)
V_leak = 10.613;% Leak reversal potential (mV)
C = 1.0;         % Membrane capacitance (uF/cm^2)
G_K = 36.0;      % Potassium maximum conductance (mS/cm^2)
G_Na = 120.0;    % Sodium maximum conductance (mS/cm^2)
G_leak = 0.3;    % Leak maximum conductance (mS/cm^2)

% % Rest of the code remains the same as in the previous response
% % Parameters for the HH model
V_baseline = -70;% Resting membrane potential (mV)
% deltaV = 10;     % Step input amplitude (mV)

% Time parameters
final_time = 50; % Simulation time (ms)
dt_ref = 0.001;  % Reference time step size (ms)

% Define the HH model equations with modified alpha and beta functions
hh_deriv = @(t, V, m, h, n, I) [
    (-G_Na * m^3 * h * (V - E_Na) - G_K * n^4 * (V - E_K) - G_leak * (V - E_leak) + I) / C;
    0.1 * ((V - 25) / (1 - exp(-(V - 25) / 10))) * (1 - m) - 4 * exp(-V / 18) * m;
    0.07 * exp(-V / 20) * (1 - h) - 1 / (1 + exp(-(V - 30) / 10)) * h;
    0.01 * ((V - 10) / (1 - exp(-(V - 10) / 10))) * (1 - n) - 0.125 * exp(-V / 80) * n
];

% Initialize arrays to store results
time_steps = [dt_ref, 0.005,0.01, 0.02, 0.03, 0.05,0.07]; % Time step sizes to analyze

% time_steps = [dt_ref, 0.05]; % Time step sizes to analyze


membrane_potentials = cell(length(time_steps), 1);

% Simulate the HH model for each time step size
for i = 1:length(time_steps)
    dt = time_steps(i);
    tspan = 0:dt:final_time;
    if i == 1
        t_plot = 0:dt:final_time;
    end


    num_steps = length(tspan);
    
    % Initialize variables
    V = zeros(num_steps, 1);
    m = zeros(num_steps, 1);
    h = zeros(num_steps, 1);
    n = zeros(num_steps, 1);
    I = zeros(num_steps, 1);
    

    % Apply step input current
    t_start = 0; % Start of the step input (ms)
    t_end = 50;   % End of the step input (ms)
    I(tspan >= t_start & tspan <= t_end) = 10; % Amplitude of the step input
    
    % Initial conditions
    % V(1) = V_baseline; % Initial membrane potential with step input
    V = 0;
    alpha_n=0.01*(10-V)/(exp((10-V)/10)-1);
    beta_n=0.125*exp(-V/80);
    alpha_m=0.1*(25-V)/(exp((25-V)/10)-1);
    beta_m=4*exp(-V/18);
    alpha_h=0.07*exp(-V/20);
    beta_h=1/(exp((30-V)/10)+1);

    n(1)=alpha_n/(alpha_n+beta_n);
    m(1)=alpha_m/(alpha_m+beta_m);
    h(1)=alpha_h/(alpha_h+beta_h);   
    
    % Time-stepping loop (Euler method)
    for j = 1:num_steps-1
        % Calculate alpha and beta values for gating variables at time j
        alpha_n = 0.01 * (10 - V(j)) / (exp((10 - V(j)) / 10) - 1);
        beta_n = 0.125 * exp(-V(j) / 80);
        alpha_m = 0.1 * (25 - V(j)) / (exp((25 - V(j)) / 10) - 1);
        beta_m = 4 * exp(-V(j) / 18);
        alpha_h = 0.07 * exp(-V(j) / 20);
        beta_h = 1 / (exp((30 - V(j)) / 10) + 1);

        % Update gating variables
        n(j+1) = n(j) + dt * (alpha_n * (1 - n(j)) - beta_n * n(j));
        m(j+1) = m(j) + dt * (alpha_m * (1 - m(j)) - beta_m * m(j));
        h(j+1) = h(j) + dt * (alpha_h * (1 - h(j)) - beta_h * h(j));

        % Calculate membrane potential using updated gating variables
        I_Na = G_Na * m(j+1)^3 * h(j+1) * (V(j) - E_Na);
        I_K = G_K * n(j+1)^4 * (V(j) - E_K);
        I_leak = G_leak * (V(j) - V_leak);
        I_ion = I(j) - I_Na - I_K - I_leak;

        V(j+1) = V(j) + dt * I_ion / C;
    end
    
    membrane_potentials{i} = V;
end

% Calculate error metrics (absolute error)
error_metrics = cell(length(time_steps), 1);
for i = 1:length(time_steps)

    V_ref = membrane_potentials{1};
    V_i = membrane_potentials{i};


    V_ref_dt = time_steps(1)/dt_ref;
    V_i_dt = time_steps(i)/dt_ref;

    error_metrics{i} = abs( V_i(1:V_ref_dt:end) - V_ref(1:V_i_dt:end) );
    
    % length(V_i(1:V_ref_dt:end) )
    % length(V_ref(1:V_i_dt:end) )

end

% Create a convergence analysis plot
% Create a figure with adjusted width and height
figure('Position', [100, 100, 800, 600]);
line_width = 2;

% subplot(1,2,1)
for i = 2:length(time_steps)
    V_i_dt = time_steps(i) / dt_ref;
    plot(t_plot(1:V_i_dt:end), error_metrics{i}, 'DisplayName', sprintf('dt = %.3f ms', time_steps(i)), 'LineWidth',line_width);
    hold on;
end
font_size = 18;
set(gca,'FontSize',font_size)
xlabel('Time (ms)', 'FontSize', font_size); % Enlarge the font size for the labels
ylabel('Absolute Error (mV)', 'FontSize', font_size); % Enlarge the font size for the labels
title('Convergence Analysis', 'FontSize', font_size); % Enlarge the font size for the title
legend('Location', 'Best', 'FontSize', font_size,'NumColumns', 2); % Enlarge the font size for the legend
grid on;

% To save the plot as an image file (optional)
% saveas(gcf, 'convergence_plot.png');


%% Verify the HH model is correct

% subplot(1,2,2)
% for i = 2:length(time_steps)
%     V_i_dt = time_steps(i) / dt_ref;
%     v_plot = membrane_potentials{i};
% 
%     % length(t_plot(1:V_i_dt:end))
%     % length(v_plot(1:V_i_dt:end))
%     plot(t_plot(1:V_i_dt:end), v_plot, 'DisplayName', sprintf('dt = %.3f ms', time_steps(i)), 'LineWidth',line_width);
%     hold on;
% end
% font_size = 18;
% set(gca,'FontSize',font_size)
% xlabel('Time (ms)', 'FontSize', font_size); % Enlarge the font size for the labels
% ylabel('V_m (mV)', 'FontSize', font_size); % Enlarge the font size for the labels
% % title('', 'FontSize', font_size); % Enlarge the font size for the title
% legend('Location', 'Best', 'FontSize', font_size); % Enlarge the font size for the legend
% hold off

