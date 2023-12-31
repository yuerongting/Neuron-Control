%%
%%
%% IDNI
clear all
% IDNI = 0;
IDNI = 1;

desired_time = 1;

ref = 15 / (1 - exp(-desired_time));

num_traject = 4;

marker_set = ["-", "--", ":", "-."];
% color_map = distinguishable_colors(20);
% color_map = color_map(3:end,:);

color_map = [[0.4940 0.1840 0.5560]; [0 0 0] ; [0.8500 0.3250 0.0980]; [0.4660 0.6740 0.1880]];

% noise 
mu = 0;
noise_std = 20/100;
noise=noise_std*randn(1,1)+mu;


lambda_1_range = [0,   3e3, 0, 1e3];
lambda_2_range = [0,    0, 4e7, 3e3];

% lambda_1_range = [8e2];
% lambda_2_range = [3e3];

% lambda_1_range = [0,   3e3, 0];
% lambda_2_range = [0,    0, 4e7];

lambda_1_range = [10,   1, 10];
lambda_2_range = [1,    10, 10];

lambda_1_range = [1,   1, 3];
lambda_2_range = [1,    3, 1];

lambda_1_range = [0.1, 1, 10];
lambda_2_range = [0.1, 1, 10];

% lambda_1_range = [5, 3,  5, 5];
% lambda_2_range = [5, 5,  3, 5];

lambda_1_range = [1e2, 0, 1, 0];
lambda_2_range = [1e1, 0, 0, 1];

lambda_1_range = [1e2, 1e1, 0];
lambda_2_range = [0, 0, 0];

lambda_1_range = [ 20,  0.1,      0 ,100 ];
lambda_2_range = [0,   0,      0,      0,0];


lambda_1_range = [ 5,  5];
lambda_2_range = [  0,  5];

lambda_1_range = [ 10, 1];
lambda_2_range = [  0, 0];

lambda_1_range = [ 10, 10, 10];
lambda_2_range = [  0, 0, 0];

lambda_1_range = [ 20, ];
lambda_2_range = [  1, ];

deltaT=0.05;


% delta_control = 0.05;
delta_control = 0.05;


% small_value = 1e-6;
% small_value = 1e-10;
% Q = small_value * eye(4);  % Small process noise covariance
% R = small_value * eye(1);  % Small measurement noise covariance

% Q = diag(cov_scale);  % Small process noise covariance
% R = diag(cov_scale(1));  % Small measurement noise covariance


font_size = 16;


color_voltage = jet(5*num_traject);
color_current = turbo(5*num_traject);

line_width = 2;

subplot(1,3,[1,2])
% subplot(2,1,1)

if(IDNI == 1)
    %% Time setting

    T_start = 2;
    Total_time=100;  %100 ms

    % Total_time=1000;  %1000 ms

    t_start=0:deltaT:T_start;
    V_start = zeros(1,numel(t_start));
    t=T_start:deltaT:(Total_time+T_start);


    V =0; % baseline voltage



    I=zeros(1,numel(t));

    % K_range = [0.2, 0.5,1];
    lambda_dynamic_inversion = 1;
    I_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));
    V_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));

    E_Na = 115; %mV
    E_K = -12; %mV
    V_leak = 10.613; %mV
    C=1;
    G_K=36;
    G_Na=120;
    G_leak = 0.3;

    %% Taylor approx 1st order
    syms m_d n_d h_d V_m b_n b_h b_m a_m a_n a_h u_e real
    x_state = [V_m n_d m_d h_d]';
    f_x = [1/C * (G_K * n_d ^4 * (E_K - V_m) + G_Na * m_d^3 * h_d * (E_Na - V_m) + G_leak * (V_leak - V_m));
        a_n * V_m * (1-n_d) - b_n * V_m * n_d;
        a_m * V_m * (1-m_d) - b_m * V_m * m_d;
        a_h * V_m * (1-h_d) - b_h * V_m * h_d] ;
    g_x = [1/C, 0,0,0]' * u_e;
    % h_x = [1,0,0,0];
    h_x = V_m;

    u_0 = 0;

    delta_u_record = zeros();
    error_record_SRDI = zeros();
    
    record_t = zeros(length(lambda_1_range), numel(t));

    time_first_15mv_IDNI = zeros(1,length(lambda_1_range));
    e_record = zeros();
    t_above_55mv_record = zeros();
    



    %% Traj setting
    total_duration = Total_time+T_start; 
    % interval_min = 8 ; 
    % action_duration = 5;
    
    T_period = 16; %% refractory period = 16 ms 

    time_vector = 0:deltaT:total_duration;
    periodic_wave = zeros(size(time_vector));
    track = periodic_wave;



    % for j = 1:num_traject    %% number of uncertain simulations
    for j = 1:length(1:1)

    % for j = 1:2  %% number of uncertain simulations
        rng(j) % random seed

        V_start = zeros(1,numel(t_start));

        %% nominal parameters
        E_Na = 115; %mV
        E_K = -12; %mV
        V_leak = 10.613; %mV
        C=1;
        G_K=36;
        G_Na=120;
        G_leak = 0.3;


        % T_start = 2;
        T_start = 0;

    %         Total_time=8;  %10 ms

        t_start=0:deltaT:T_start;

        t=T_start:deltaT:(Total_time+T_start);

        lambda_dynamic_inversion = 1;
        I_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));
        V_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));

        r_record = zeros();
        r_d_record = zeros();

        
        




        V =0; % baseline voltage (displacement)

        I=zeros(1,numel(t));
        %% Uncertain model
        percent = 0.0;

        % if(j>1)
        %     percent = uncertain_percent;
        % end

        lambda_1 = lambda_1_range(1);
        lambda_2 = lambda_2_range(1); % almost no contribution to error in [0,1]

        % Uncertainty
        E_Na = normrnd(E_Na, E_Na * percent);
        E_K = normrnd(E_K, abs(E_K) * percent);
        V_leak = normrnd(V_leak, V_leak * percent);
        C = normrnd(C, C * percent);
        G_K = normrnd(G_K, G_K * percent);
        G_Na = normrnd(G_Na, G_Na * percent);
        G_leak = normrnd(G_leak, G_leak * percent);


        alpha_n=0.01*(10-V)/(exp((10-V)/10)-1);
        beta_n=0.125*exp(-V/80);
        alpha_m=0.1*(25-V)/(exp((25-V)/10)-1);
        beta_m=4*exp(-V/18);
        alpha_h=0.07*exp(-V/20);
        beta_h=1/(exp((30-V)/10)+1);

        n(1)=alpha_n/(alpha_n+beta_n);
        m(1)=alpha_m/(alpha_m+beta_m);
        h(1)=alpha_h/(alpha_h+beta_h);   

        I(1) = 0;
        
        %% Extended Kalman Filter
        
        % % State transition function
        % stateTransition = @(x, I) stateTransitionFunction(x, I, deltaT, G_Na, G_K, G_leak, E_Na, E_K, V_leak, C);
        % 
        % % Measurement function (for simplicity, assuming we can directly measure the membrane potential)
        % measurementFunction = @(x) x(1);
        % 
        % % % Process noise covariance (you may need to adjust this based on your scenario)
        % % Q = eye(4) * 0.01;
        % % % Measurement noise covariance
        % % R = 0.1;
        % 
        % small_value = 1e-6;
        % cov_scale = [0.1, 0.01, 0.01, 0.01];
        % 
        % 
        % % Initial state vector
        % x_0 = [V(1); n(1); m(1); h(1)];
        % 
        % % Initialize the EKF
        % ekf = extendedKalmanFilter('StateTransitionFcn', stateTransition, 'MeasurementFcn', measurementFunction, 'ProcessNoise', Q, 'MeasurementNoise', R);
        % 
        % ekf.State = x_0; %% Initial state 
        % ekf.StateCovariance = eye(4); %% Initial covariance 
        % 
        % % Arrays to store results
        % estimated_states = zeros(4, numel(t)-1);
        
        
        %% Trajectory for each "j"
        for i=1:numel(t)-1
            time = i * deltaT;
            alpha_n(i) = .01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) );
            beta_n(i) = .125*exp(-V(i)/80);
            alpha_m(i) = .1*( (25-V(i)) / (exp((25-V(i))/10)-1) );
            beta_m(i) = 4*exp(-V(i)/18);
            alpha_h(i) = .07*exp(-V(i)/20);
            beta_h(i) = 1/(exp((30-V(i))/10)+1);

            % V_th = -55;

            %% System Linearized at x_0, u_0

            x_state = [V_m n_d m_d h_d]';
            
            A_0 = jacobian((f_x+g_x),x_state);
            B_0 = jacobian(g_x,u_e);

            A_0 = subs(A_0, {V_m, n_d, m_d, h_d,      a_n, a_m, a_h,     b_n, b_m, b_h,     u_e}, ...
                {V(i), n(i), m(i), h(i),      alpha_n(i), alpha_m(i), alpha_h(i),   beta_n(i), beta_m(i),   beta_h(i)    u_0});

            B_0 = subs(B_0, {V_m, n_d, m_d, h_d,      a_n, a_m, a_h,     b_n, b_m, b_h,     u_e}, ...
                {V(i), n(i), m(i), h(i),      alpha_n(i), alpha_m(i), alpha_h(i),   beta_n(i), beta_m(i),   beta_h(i)    u_0});

            f_x_0 = subs(f_x,{V_m, n_d, m_d, h_d,      a_n, a_m, a_h,     b_n, b_m, b_h,     u_e}, ...
                {V(i), n(i), m(i), h(i),      alpha_n(i), alpha_m(i), alpha_h(i),   beta_n(i), beta_m(i),   beta_h(i)    u_0 });

            u_0 = I(i);
            x_0_dot = f_x_0 + subs(g_x,{u_e},{u_0});

            
            %% Ref trajectory design

            
            
            % T_period = interval_min + action_duration;  % Duration of one period (10 seconds)

            % r_d = ref * exp( -1 * (i-1) * deltaT);
            % r_d_record(i) = r_d;
            % 
            % track(i) = ref * (1 - exp( -1 *(i-1) * deltaT));
            
            index_traj = i - 1;
            % index_traj = i ;
            % r_d = ref * exp( -1 * (index_traj) * deltaT);
            % r_d_record(i) = r_d;

            % track(i) = ref * (1 - exp( -1 *(index_traj) * deltaT));

            % time_in_current_period = mod(  time_vector(index_traj)  ,  T_period )  ;
            time_in_current_period = mod(  time  ,  T_period + desired_time )  ;

            % Add spikes within each 10-second interval
            if time_in_current_period <= desired_time % Change 1.0 to adjust spike width
                spike_value = ref * (1 - exp( -1 * time_in_current_period)); % Adjust spike width as needed
                
                % r_d = ref * exp( -1 * (time_in_current_period) * deltaT);
                % r_d_record(i) = r_d;

                % periodic_wave(index_traj) = spike_value ;
                % display(time_in_current_period)
                r_d = ref * exp( -1 * (time_in_current_period) * deltaT);
                
                track(i) = spike_value;
                current_flag = 1;
            else 
                r_d = 0;
                current_flag = 0;
            end

            r_d_record(i) = r_d;

            r = track(i);
            r_record(i) = r;

            % if (i > desired_time/ deltaT)
            %    track(i) = 15;
            %    % record_t(i) = i;
            % end

            interval = 70 - 55;

            if (V(i) > interval)
               record_t(j,i) = i;
            end            
            first_nonzero_index = find(record_t(j,1:end) ~= 0, 1, 'first');%% Time when reaching the threshold


            x_state = [V(i), n(i), m(i), h(i)]';

            g_x_0 = B_0;

            x_dot_0 = x_0_dot;





            %%% measurement also uncertain
            hh = [1,0,0,0] * x_state;  
            

            %%% Error
            error = r - hh;

            e_record(j,i) = error * current_flag;


            %% Control sampling time
            if( mod(i+1, delta_control/deltaT) == 0)
            % if( (mod(i+1, delta_control/deltaT) == 1) && (r_d ~= 0) )
            %     delta_u = inv([1,0,0,0] * subs(B_0,{u_e},{u_0})  ) * ( r_d - [1,0,0,0]*(  x_dot_0 + (A_0 + g_x_0 * u_0) * (x_state - )) - lambda_1*(r - hh ) - lambda_2 * (r-hh)^3);
                % delta_u = inv([1,0,0,0] * subs(B_0,{u_e},{u_0})  ) * ( r_d - [1,0,0,0]*(  x_dot_0 ) - lambda_1*(r - hh ) - lambda_2 * (r-hh)^3);
                
                % delta_u = inv([1,0,0,0] * subs(B_0,{u_e},{u_0})  ) * ( r_d - [1,0,0,0]*(  x_dot_0 ) -  (lambda_1*(r - hh )) - (lambda_2 * (r-hh)^3)   );
                delta_u = current_flag * inv([1,0,0,0] * subs(B_0,{u_e},{u_0}) )  * ( r_d  + ( lambda_1*(r - hh ) ) + (lambda_2 * (r-hh)^3)  - [1,0,0,0]*(  x_dot_0 )  )  ;
                
                % double(delta_u * deltaT)
                
                display(strcat('i  :', num2str(double(i))))

                display(strcat('r_d :', num2str(double(r_d))))

                % r_d
                display(strcat('Error  :', num2str(  double(error)  )))

                display(strcat('Delta_u  :', num2str(  double(delta_u * deltaT)  )))

                display(strcat('x_dot_0  :', num2str(  double(x_dot_0(1)  ))))
                display(' ')

    % delta_u = inv([1,0,0,0] * subs(B_0,{u_e},{u_0})  ) * ( r_d - [1,0,0,0]*(  x_dot_0 ) - lambda_1*(r - hh ) - lambda_2 * (r-hh)^3);
    
                % disp(lambda_1)
                % disp(lambda_2)
                % delta_u_record(i) = delta_u;
    
                % Incremental input

                % I(i+1) = I(i) + delta_u * deltaT; %% Already contain
                % deltaT in "V(i+1) = V(i) + deltaT*I_ion/C;"

                I(i+1) = I(i) + delta_u;
            else
                delta_u = 0;
                I(i+1) = I(i);
            end
            

            delta_u_record(j,i) = delta_u;


            %% Noise term

            % noise=std*randn(1,1)+mu;
            % 
            % if I(i) == 0
            %     noise = 0;
            % end
            % 
            % noise = 0;  %%%%% without noise

            % I(i) = I(i) + noise;


            % Constrained input
            

            % if( (V(i)> interval) || V(i)<0 || time > 3 ) %% Positive current
            if (~isempty(first_nonzero_index))

                % if( (V(i)> interval) || V(i)<0 || i > first_nonzero_index) %% Positive current
                if( (V(i)> interval) ) %% Positive current
            % %     if((V(i)> interval) || V(i)<0)  %% Negative current
                   I(i+1) = 0;
                end
            end

            if( I(i+1) < 0 ) %% Positive current
        % %     if((V(i)> interval) || V(i)<0)  %% Negative current
               I(i+1) = 0;
            end

            % if( time_in_current_period >  desired_time )
            %     I(i+1) =0;
            % end

            % I(i+1) = 100;

            if I(i+1) == 0
                noise = 0;
            end
            % noise = 0;
%             noise = 0;  %%%%% test

            I(i+1) = I(i+1) + noise;


            if (V(i)< interval && time < 3) % record error
                error_record_SRDI(j,i+1) = error;
            end
            


            %% Step input after spike
            % I(i) = 20 * j;

            %% Model Dynamics 

            I_Na = double((m(i)^3)) *G_Na * double(h(i)) * (V(i)-E_Na); %Equations 3 and 14
            I_K = double((n(i)^4)) * G_K * (V(i)-E_K); %Equations 4 and 6
            I_leak = G_leak * (V(i)-V_leak);
            I_ion = I(i+1) - I_K - I_Na - I_leak;
            

            V(i+1) = V(i) + deltaT*I_ion/C;
            n(i+1) = n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i));
            m(i+1) = m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i));
            h(i+1) = h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i));
            

            
            %% Extended Kalman Filter

            % % Step 1: Update membrane potential and gating variables using Hodgkin-Huxley model
            % % [V, n, m, h] = hodgkinHuxleyUpdate(V, n, m, h, ...);  % Replace ... with necessary parameters
            % % V(i+1)
            % 
            % % Step 2: EKF Prediction Step
            % [x_kplus1, ~] = predict(ekf, deltaT);
            % estimated_states(:, i) = x_kplus1;
            % 
            % % Step 3: Obtain measurement (simulated or real)
            % V_measurement = V(i+1) + randn() * sqrt(R);  % Noisy measurement
            % 
            % % Step 4: EKF Correction Step
            % [x_kplus1_corrected, ~] = correct(ekf, V_measurement);
            % 
            % % Optionally, update the state for the next iteration
            % ekf.State = x_kplus1_corrected;
            % 
            % 
            % if i > first_nonzero_index
            %     V(i+1) = V(i) + deltaT*I_ion/C;
            %     n(i+1) = n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i));
            %     m(i+1) = m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i));
            %     h(i+1) = h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i));
            % else
            % 
            % % Update state for the next iteration
            % % estimated_states(:, i) = x_kplus1; % Store the estimated states
            %     V(i+1) = x_kplus1_corrected(1);
            %     n(i+1) = x_kplus1_corrected(2);
            %     m(i+1) = x_kplus1_corrected(3);
            %     h(i+1) = x_kplus1_corrected(4);
            % end

        end

            display(j);


            V_rest = -70;
            % V_rest = 0;

            V = V+V_rest; %Set resting potential to -70mv

            I_record(1, :) = I;
            V_record(1, :) = V;

            I_start = V_start;
            V_start = V_start + V_rest;

            % t_above_55mv = find(V_record>=-55);
            % t_above_55mv = (t_above_55mv(1) -1) *deltaT + T_start;  %%  start from index = 1 so minus 1
            % t_above_55mv_record(j) = t_above_55mv;

            % time_first_15mv_IDNI(1,j) =  find( V > -55, 1 ) *  deltaT ;
%             display(time_first_15mv_IDNI(1,j));
%             display(find( V > -55, 1 ) *  deltaT );
            
            
            %% Nominal model

            if j==1
                V_record_nominal = V_record;
                I_record_nominal = I_record;
            end

            yyaxis left
            if j>1
                yyaxis left
                % plot(cat(2,t_start,t),cat(2,V_start,V_record(1,:)), marker_set(1), 'LineWidth',line_width, 'color', color_map(j,:),'DisplayName',strcat( ' Membrane Potential Voltage ', '\lambda_1=',num2str(-lambda_1_range(j)), ' ,\lambda_2=',num2str(-lambda_2_range(j))));
                plot(t, V_record(1,:), marker_set(1), 'LineWidth',line_width, 'color', color_map(j,:),'DisplayName',strcat( ' Membrane Potential Voltage ', '\lambda_1=',num2str(lambda_1_range(j)), ' ,\lambda_2=',num2str(lambda_2_range(j))));
                yyaxis right
                plot(t, I_record(1,:), marker_set(3), 'LineWidth',line_width, 'color',color_current(j,:),'DisplayName',strcat(' Current Injection ', '\lambda_1=',num2str(lambda_1_range(j)), ' ,\lambda_2=',num2str(lambda_2_range(j))));hold on
            else
                yyaxis left
                plot(t, V_record(1,:),marker_set(1),'LineWidth',line_width, 'color', color_map(j,:), 'DisplayName',strcat(' Membrane Potential Voltage ','\lambda_1=',num2str(lambda_1_range(j)), ' ,\lambda_2=',num2str(lambda_2_range(j))));hold on
                yyaxis right
                plot(t, I_record(1,:),marker_set(3),'LineWidth',line_width, 'color',color_current(j,:),'DisplayName',strcat(' Current Injection ','\lambda_1=',num2str(lambda_1_range(j)), ' ,\lambda_2=',num2str(lambda_2_range(j))))
                % plot(cat(2,t_start,t),cat(2,I_start,I_record(1,:)),'LineWidth',line_width, 'color',distinguishable_colors(j, :),'DisplayName',strcat(' Current Injection ','\lambda_1=',num2str(-lambda_1_range(j)), ' ,\lambda_2=',num2str(-lambda_2_range(j))))
                
            end
            hold on

    end

    
    %% Plot
    
    yyaxis left
    yline(-55,'--b', 'LineWidth',3, 'DisplayName', ' Threshold', 'color', [0 0.4470 0.7410]) % Threshold

    track_plot = track(1:end-1) + V_rest;

    set(gca,'FontSize',font_size)
    color_simu = [0 0.4470 0.7410];



    % i = 1;
    hold on
    yyaxis left
    % plot(cat(2,t_start,t),cat(2,V_start,V_record_nominal(i,:)),'--','LineWidth',1, 'Color', '#00008B') % uncertain model label

    % add label for uncertain model
    % plot(cat(2,t_start,t),cat(2,V_start,V_record(1,:)), '--','LineWidth',1, 'color', '#00008B','HandleVisibility','off');
    

    
    plot( t(1:end-2) , track_plot(1:length(t(1:end-2))) , 'color', "#FF00FF",'LineWidth', 1, 'Marker', '.', 'DisplayName', ' Reference Trajectory')  %% ref track
    hold on
    ylabel('Voltage (mV)','FontSize',font_size)
    xlabel('Time (ms)','FontSize',font_size)



    yyaxis right
    % plot(cat(2,t_start,t),cat(2,I_start,I_record_nominal(i,:)),'LineWidth',1, 'color','r')

    hold on
    ylabel('Current (nA)' ,'FontSize',font_size)

    grid on
    % end 

    % legend('show')
    legend('Location','northeastoutside'); 

    % title('SRDI','FontSize',font_size)
    % legend( 'Reference Trajectory', 'Threshold (-55mV)' , 'Voltage (Uncertain Model)', 'Current Injection','Input Contraint (15 nA)', 'FontSize', font_size)

    % h_legend = legend( 'Voltage (Nominal Model)', 'Threshold (-55mV)' , 'Voltage (Uncertain Model)', 'Reference Trajectory', 'Current (Nominal Model)','Current (Uncertain Model)', 'FontSize', font_size);
    % set(h_legend, 'location', 'northeastoutside')
    
    % xlim([1.5, 6])

    % 'Threshold (-55mV)' ,
    hold on
end
% hold off

subplot(1,3,3)

% for j = 1:num_traject
for j = 1:length(lambda_1_range)
    plot(t(2:end),abs(e_record(j,:)).', 'LineWidth',line_width, 'color',color_map(j,:),'DisplayName',strcat(' Current Injection ','\lambda_1=',num2str(lambda_1_range(j)), ' ,\lambda_2=',num2str(lambda_2_range(j))))
    hold on
end

set(gca,'FontSize',font_size)
% xlim([T_start, T_start + T_period-1])
% xlim([T_start, T_start + 5])
ylabel('Error (mV)','FontSize',font_size)
xlabel('Time (ms)','FontSize',font_size)
% legend([strcat(' Current Injection ','\lambda_1=',num2str(-lambda_1_range(j)), ' ,\lambda_2=',num2str(-lambda_2_range(j)))])
legend()
% hold off