close all;clear all;clc;

num_traject = 30; 
% num_traject = 5; 
desired_time = 0.8;

ref = 15 / (1 - exp(-desired_time));

uncertain_percent = 0.1;
font_size = 23;

deltaT=0.05;
T_start = 0;
Total_time=2;  %2 ms

% noise 
mu = 0;
noise_std = 20/100;
noise=noise_std*randn(1,1)+mu;

% noise = 0; %% test


%(a) dynamic inversion
Dynamic_inversion = 1; 
lambda_dynamic_inversion = 30; % Dynamic inversion


%(b) IDNI
IDNI = 1; 
lambda_1 = 20;
lambda_2 = 1.0; % almost no contribution to error in [0,1]

%(C) MPC
MPC = 1; 
predict_horizon = 20; %% Predict horizon
control_horizon = 20; %% Control horizon

time_first_15mv_DI = zeros(1,num_traject);

t=T_start:deltaT:(Total_time+T_start);
for i=1:numel(t)
    track(i) = ref * (1 - exp( -1 *(i-1)*deltaT));
    if (i > desired_time/ deltaT)
       track(i) = 15;
    end
end



if(Dynamic_inversion == 1)
    subplot(3,3,[1,2])
    E_Na = 115; %mV
    E_K = -12; %mV
    V_leak = 10.613; %mV
    C=1;
    G_K=36;
    G_Na=120;
    G_leak = 0.3;

    % % deltaT=0.01;
    % deltaT=0.05;
    % 
    % % T_start = 2;
    % Total_time=2;  %10 ms

    t_start=0:deltaT:T_start;
    V_start = zeros(1,numel(t_start));
    t=T_start:deltaT:(Total_time+T_start);

    V =0; % baseline voltage
    alpha_n=0.01*(10-V)/(exp((10-V)/10)-1);
    beta_n=0.125*exp(-V/80);
    alpha_m=0.1*(25-V)/(exp((25-V)/10)-1);
    beta_m=4*exp(-V/18);
    alpha_h=0.07*exp(-V/20);
    beta_h=1/(exp((30-V)/10)+1);
    n(1)=alpha_n/(alpha_n+beta_n);
    m(1)=alpha_m/(alpha_m+beta_m);
    h(1)=alpha_h/(alpha_h+beta_h);
    I=zeros(1,numel(t));

    I_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));
    V_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));
    error_record_CDI = zeros();
    error_time = zeros();
    

    
%     track(52:end) = track(51);
    
    for j = 1:num_traject

        rng(j) % random seed

        t_start=0:deltaT:T_start;
        V_start = zeros(1,numel(t_start));
        t=T_start:deltaT:(Total_time+T_start);

        for k = 1:numel(lambda_dynamic_inversion)
            K = lambda_dynamic_inversion(k);
            
            
            
            V =0; % baseline voltage
            alpha_n=0.01*(10-V)/(exp((10-V)/10)-1);
            beta_n=0.125*exp(-V/80);
            alpha_m=0.1*(25-V)/(exp((25-V)/10)-1);
            beta_m=4*exp(-V/18);
            alpha_h=0.07*exp(-V/20);
            beta_h=1/(exp((30-V)/10)+1);
            n(1)=alpha_n/(alpha_n+beta_n);
            m(1)=alpha_m/(alpha_m+beta_m);
            h(1)=alpha_h/(alpha_h+beta_h);
            I=zeros(1,numel(t));


            E_Na = 115; %mV
            E_K = -12; %mV
            V_leak = 10.613; %mV
            C=1;
            G_K=36;
            G_Na=120;
            G_leak = 0.3;

            %% Uncertain model
            if(j==1)
                percent = 0;
            end

            if(j>1)
                percent = uncertain_percent;
            end

            % Uncertainty
            E_Na = normrnd(E_Na, E_Na * percent);
            E_K = normrnd(E_K, abs(E_K) * percent);
            V_leak = normrnd(V_leak, V_leak * percent);
            C = normrnd(C, C * percent);
            G_K = normrnd(G_K, G_K * percent);
            G_Na = normrnd(G_Na, G_Na * percent);
            G_leak = normrnd(G_leak, G_leak * percent);

            for i=1:numel(t)-1

                alpha_n(i) = .01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) );
                beta_n(i) = .125*exp(-V(i)/80);
                alpha_m(i) = .1*( (25-V(i)) / (exp((25-V(i))/10)-1) );
                beta_m(i) = 4*exp(-V(i)/18);
                alpha_h(i) = .07*exp(-V(i)/20);
                beta_h(i) = 1/(exp((30-V(i))/10)+1);

                %% Dynamic inversion with Constain input current
    %             V_th = -55; % track

%                 track(i) = ref * (1 - exp( -1 *(i-1)*deltaT));
                V_th = track(i);

                I(i) = K * C* V_th - G_K*(n(i)^4)*(E_K-V(i)) - G_Na*(m(i)^3)*h(i)*(E_Na-V(i)) - G_leak*(V_leak-V(i)) - K * V(i);
    %             I(i) = -I(i);

              
                noise=noise_std*randn(1,1)+mu;

%                 noise = 0; %% test

                I(i) = I(i) + noise;

                time = i * deltaT;

                interval = 70 - 55;
                
                if (V(i)< interval && time < 3) % record error
                    error_record_CDI(j,i) = V(i) - track(i);
%                     error_time(j) = i;

                end
                
                if( (V(i)> interval) || V(i)<0 || time > 3 || I(i) < 0) %% Positive current
            % %     if((V(i)> interval) || V(i)<0)  %% Negative current
                   I(i) = 0;
                end
                
                

                %% Model Dynamics
                I_Na = (m(i)^3) *G_Na * h(i) * (V(i)-E_Na); %Equations 3 and 14
                I_K = (n(i)^4) * G_K * (V(i)-E_K); %Equations 4 and 6
                I_leak = G_leak * (V(i)-V_leak);
                I_ion = I(i) - I_K - I_Na - I_leak;

                V(i+1) = V(i) + deltaT*I_ion/C;
                n(i+1) = n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i));
                m(i+1) = m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i));
                h(i+1) = h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i));

            end

            
            V_rest = -70;
            % V_rest = 0;

            V = V+V_rest; %Set resting potential to -70mv



            I_record(j, :) = I;
            V_record(j, :) = V;
            if(~isempty(find( V > -55, 1 ) ))
                time_first_15mv_DI(1,j) =  find( V > -55, 1 ) *  deltaT ;
            end


            I_start = V_start;
            V_start = V_start + V_rest;


            %% Nominal model

            if j==1
                V_record_nominal = V_record;
                I_record_nominal = I_record;
            end
            
            colororder({'r','b'})
            % yyaxis right
            % if j>1
            %     yyaxis right
            %     plot(cat(2,t_start,t),cat(2,V_start,V_record(j,:)), '--','LineWidth',1, 'color', '#00008B','HandleVisibility','off');
            %     yyaxis left
            %     plot(cat(2,t_start,t),cat(2,I_start,I_record(j,:)), '--','LineWidth',1, 'color','r','HandleVisibility','off');
            % else
            %     yyaxis right
            %     plot(cat(2,t_start,t),cat(2,V_start,V_record(j,:)),'LineWidth',3, 'color', [0 0.4470 0.7410]);hold on
            %     yyaxis left
            %     plot(cat(2,t_start,t),cat(2,I_start,I_record(j,:)),'LineWidth',3, 'color','r')
            % end
            if j==1
                yyaxis right
                plot(cat(2,t_start,t),cat(2,V_start,V_record(j,:)),'LineWidth',3, 'color', [0 0.4470 0.7410]);hold on
                yyaxis left
                plot(cat(2,t_start,t),cat(2,I_start,I_record(j,:)),'LineWidth',3, 'color','r')
            end
            % if j==1
            %     yyaxis right
            %     plot(cat(2,t_start,t),cat(2,V_start,V_record(j,:)),'LineWidth',3, 'color', [0 0.4470 0.7410]);hold on
            %     yyaxis left
            %     plot(cat(2,t_start,t),cat(2,I_start,I_record(j,:)),'LineWidth',3, 'color','r')
            % end
            hold on
        end
    end

    
    yyaxis right
    yline(-55,'--b', 'LineWidth',3) % Threshold

    track_plot = track + V_rest;

    set(gca,'FontSize',font_size)
    color_simu = [0 0.4470 0.7410];
    


    i = 1;
    hold on
    yyaxis right
    % plot(cat(2,t_start,t),cat(2,V_start,V_record_nominal(i,:)),'--','LineWidth',1, 'Color', '#00008B') % uncertain model label

    % add label for uncertain model
    % plot(cat(2,t_start,t),cat(2,V_start,V_record(1,:)), '--','LineWidth',1, 'color', '#00008B','HandleVisibility','off');
    
    % ref traj
    plot(t(1:end-1), track_plot(1,1:end-1), 'color', "#FF00FF",'LineWidth', 2, 'Marker', '_')  %% ref track
    
    %% plot confidence interval for V
    min_values = min(V_record,[],1);
    max_values = max(V_record,[],1);
    
    % Calculate the mean values for each trajectory (optional, for reference)
    mean_values = mean(V_record, 1);
    
    % Calculate the standard error of the mean for each trajectory
    n = size(V_record, 2);  % Number of samples per trajectory
    standard_error = std(V_record, 0, 1) / sqrt(n);
    confidence_interval = 1.96 * standard_error;  % 1.96 is the critical value for a 95% confidence interval
    
    % Create vectors for the ascending and descending x values
    x = 1:size(V_record, 2);  % x-axis represents trajectories
    x = x-1;
    x_ascend = x;
    x_descend = flip(x);
    % Create a loop to fill the confidence intervals for each trajectory
    min_ascend = mean_values - confidence_interval;
    max_descend = flip(mean_values + confidence_interval);
    
    fill([x_ascend, x_descend] * deltaT, [min_ascend, max_descend], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    fill([x_ascend, x_descend] * deltaT, [min_values, flip(max_values)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    %% Continue plotting
    hold on;

    % ylabel('Current (nA)','FontSize',font_size)
    ylabel('Voltage (mV)','FontSize',font_size)
    % xlabel('Time (ms)','FontSize',font_size)

    % set(gca,'color','r');
    % gca.YAxis(1).Color = 'b';
    % gca.YAxis(2).Color = 'b';
    

    yyaxis left
    % plot(cat(2,t_start,t),cat(2,I_start,I_record_nominal(i,:)),'LineWidth',1, 'color','r')
    hold on

    %% plot confidence interval for I
    min_values = min(I_record,[],1);
    max_values = max(I_record,[],1);

    % Calculate the mean values for each trajectory (optional, for reference)
    mean_values = mean(I_record, 1);

    % Calculate the standard error of the mean for each trajectory
    n = size(I_record, 2);  % Number of samples per trajectory
    standard_error = std(I_record, 0, 1) / sqrt(n);
    confidence_interval = 1.96 * standard_error;  % 1.96 is the critical value for a 95% confidence interval
    
    % Create vectors for the ascending and descending x values
    x = 1:size(I_record, 2);  % x-axis represents trajectories
    x = x-1;
    x_ascend = x;
    x_descend = flip(x);
    % Create a loop to fill the confidence intervals for each trajectory
    min_ascend = mean_values - confidence_interval;
    max_descend = flip(mean_values + confidence_interval);

    %% 95 confidence region 
    fill([x_ascend, x_descend] * deltaT, [min_ascend, max_descend], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    %% All possible region
    fill([x_ascend, x_descend] * deltaT, [min_values, flip(max_values)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    %% Continue plotting
    hold on
    % ylabel('Voltage (mV)' ,'FontSize',font_size)
    ylabel('Current (nA)' ,'FontSize',font_size)
    
    % leftLabel = get(gca, 'ylabel');  % Get the label of the left yyaxis
    % rightLabel = get(gca, 'ylabel', 'right');  % Get the label of the right yyaxis
    % 
    % set(gca, 'ylabel', rightLabel);  % Set the label of the left yyaxis to the label of the right yyaxis
    % set(gca, 'ylabel', leftLabel, 'YAxisLocation', 'right');  % Set the label of the right yy

    grid on
    % end 

    title('(a) CDI','FontSize',font_size)
    % legend( 'Reference Trajectory', 'Threshold (-55mV)' , 'Voltage (Uncertain Model)', 'Current Injection','Input Contraint (15 nA)', 'FontSize', font_size)

    % h_legend = legend( 'Voltage (Nominal Model)', 'Threshold (-55mV)' , 'Reference Trajectory', 'Voltage (95% Confidence Region)', 'Voltage (All Possible Region)','Current (Nonimal Model)','Current (95% Confidence Region)', 'Current (All Possible Region)', 'FontSize', font_size);
    % set(h_legend, 'location', 'northeastoutside')
    
    % 'Threshold (-55mV)' ,
    hold on


    
end












%%
%%
%% IDNI

% IDNI = 0;

if(IDNI == 1)
    %% Time setting
    if (MPC==1)
        subplot(3,3,[4,5])
    end

    % deltaT=0.01;

    % T_start = 2;
    % Total_time=4;  %10 ms

    t_start=0:deltaT:T_start;
    V_start = zeros(1,numel(t_start));
    t=T_start:deltaT:(Total_time+T_start);


    V =0; % baseline voltage


    alpha_n=0.01*(10-V)/(exp((10-V)/10)-1);
    beta_n=0.125*exp(-V/80);
    alpha_m=0.1*(25-V)/(exp((25-V)/10)-1);
    beta_m=4*exp(-V/18);
    alpha_h=0.07*exp(-V/20);
    beta_h=1/(exp((30-V)/10)+1);
    n(1)=alpha_n/(alpha_n+beta_n);
    m(1)=alpha_m/(alpha_m+beta_m);
    h(1)=alpha_h/(alpha_h+beta_h);
    I=zeros(1,numel(t));

    % K_range = [0.2, 0.5,1];
    lambda_dynamic_inversion = 1;
    % I_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));
    % V_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));

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

    A_0 = jacobian((f_x+g_x),x_state);
    B_0 = jacobian(g_x,u_e);

    V =0;
    i=1;
    alpha_n(i) = .01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) );
    beta_n(i) = .125*exp(-V(i)/80);
    alpha_m(i) = .1*( (25-V(i)) / (exp((25-V(i))/10)-1) );
    beta_m(i) = 4*exp(-V(i)/18);
    alpha_h(i) = .07*exp(-V(i)/20);
    beta_h(i) = 1/(exp((30-V(i))/10)+1);

    n(1)=alpha_n/(alpha_n+beta_n);
    m(1)=alpha_m/(alpha_m+beta_m);
    h(1)=alpha_h/(alpha_h+beta_h);
    alpha_n=0.01*(10-V)/(exp((10-V)/10)-1);
    beta_n=0.125*exp(-V/80);
    alpha_m=0.1*(25-V)/(exp((25-V)/10)-1);
    beta_m=4*exp(-V/18);
    alpha_h=0.07*exp(-V/20);
    beta_h=1/(exp((30-V)/10)+1);

    % if ( i == 1) 
    x_0 = [V(i), n(i), m(i), h(i)]';
    u_0 = 0;


    i=1;

    A_0 = subs(A_0, {V_m, n_d, m_d, h_d,      a_n, a_m, a_h,     b_n, b_m, b_h,     u_e}, ...
        {V(i), n(i), m(i), h(i),      alpha_n(i), alpha_m(i), alpha_h(i),   beta_n(i), beta_m(i),   beta_h(i)    u_0});

    % 
    B_0 = subs(B_0, {V_m, n_d, m_d, h_d,      a_n, a_m, a_h,     b_n, b_m, b_h,     u_e}, ...
        {V(i), n(i), m(i), h(i),      alpha_n(i), alpha_m(i), alpha_h(i),   beta_n(i), beta_m(i),   beta_h(i)    u_0});

    f_x_0 = subs(f_x,{V_m, n_d, m_d, h_d,      a_n, a_m, a_h,     b_n, b_m, b_h,     u_e}, ...
        {V(i), n(i), m(i), h(i),      alpha_n(i), alpha_m(i), alpha_h(i),   beta_n(i), beta_m(i),   beta_h(i)    u_0 });


    % deltaT=0.01;


    % T_start = 2;
    % Total_time=4;  %10 ms

    t_start=0:deltaT:T_start;

    r_record = zeros();
    r_d_record = zeros();

    delta_u_record = zeros();
    error_record_SRDI = zeros();

    time_first_15mv_IDNI = zeros(1,num_traject);
    
    I_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));
    V_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));


    for j = 1:num_traject  %% number of uncertain simulations

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
    %         Total_time=8;  %10 ms

        t_start=0:deltaT:T_start;

        t=T_start:deltaT:(Total_time+T_start);

        lambda_dynamic_inversion = 1;
        

        V =0; % baseline voltage

        I=zeros(1,numel(t));
        %% Uncertain model
        if(j==1)
            percent = 0;
        end

        if(j>1)
            percent = uncertain_percent;
        end

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

        for i=1:numel(t)-1
            time = i * deltaT;

            alpha_n(i) = .01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) );
            beta_n(i) = .125*exp(-V(i)/80);
            alpha_m(i) = .1*( (25-V(i)) / (exp((25-V(i))/10)-1) );
            beta_m(i) = 4*exp(-V(i)/18);
            alpha_h(i) = .07*exp(-V(i)/20);
            beta_h(i) = 1/(exp((30-V(i))/10)+1);

            V_th = -55;
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

            r_d = ref * exp( -1 * (i-1)*deltaT);
            r_d_record(i) = r_d;

        %     track(i) = y_des_dot - ref;

            % track(i) = ref * (1 - exp( -1 *(i-1)*deltaT));
            % 
            % if (i > desired_time/ deltaT)
            %    track(j,i) = 15;
            %    % record_t(i) = i;
            % end

            x_state = [V(i), n(i), m(i), h(i)]';

            g_x_0 = B_0;

            x_dot_0 = x_0_dot;

            r = track(i);
            r_record(i) = r;

            hh = [1,0,0,0] * x_state;  %%% measurement also uncertain

        %     delta_u = inv([1,0,0,0] * subs(B_0,{u_e},{u_0})  ) * ( r_d - [1,0,0,0]*(  x_dot_0 + (A_0 + g_x_0 * u_0) * (x_state - )) - lambda_1*(r - hh ) - lambda_2 * (r-hh)^3);
            % delta_u = inv([1,0,0,0] * subs(B_0,{u_e},{u_0})  ) * ( r_d - [1,0,0,0]*(  x_dot_0 ) - lambda_1*(r - hh ) - lambda_2 * (r-hh)^3);
            delta_u = inv([1,0,0,0] * subs(B_0,{u_e},{u_0}) )  * ( r_d  + ( lambda_1*(r - hh ) ) + (lambda_2 * (r-hh)^3)  - [1,0,0,0]*(  x_dot_0 )  )  ;
                
            delta_u_record(i) = delta_u;

            % Incremental input
            I(i+1) = I(i) + delta_u;

            %% Noise term

            noise=noise_std*randn(1,1)+mu;

            if I(i+1) == 0
                noise = 0;
            end
            % noise = 0;
%             noise = 0;  %%%%% test

            I(i+1) = I(i+1) + noise;

            if( I(i+1) < 0 ) %% Positive current
        % %     if((V(i)> interval) || V(i)<0)  %% Negative current
               I(i+1) = 0;
            end
            % Constrained input
            interval = 70 - 55;
        %     if((i>1)&& (V(i)<V(i-1)) ||(V(i)> interval) || V(i)<0) %% Positive current

            if( (V(i)> interval) || V(i)<0 || time > 3 ) %% Positive current
        % %     if((V(i)> interval) || V(i)<0)  %% Negative current
               I(i+1) = 0;
            end

            if (V(i)< interval && time < 3) % record error
                error_record_SRDI(j,i) = r - hh;
            end
            


            %% Model Dynamics

            I_Na = double((m(i)^3)) *G_Na * double(h(i)) * (V(i)-E_Na); %Equations 3 and 14
            I_K = double((n(i)^4)) * G_K * (V(i)-E_K); %Equations 4 and 6
            I_leak = G_leak * (V(i)-V_leak);
            I_ion = I(i+1) - I_K - I_Na - I_leak;


            V(i+1) = V(i) + deltaT*I_ion/C;
            n(i+1) = n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i));
            m(i+1) = m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i));
            h(i+1) = h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i));


        end
            display(j);


            V_rest = -70;
            % V_rest = 0;

            V = V+V_rest; %Set resting potential to -70mv

            I_record(j, :) = I;
            V_record(j, :) = V;

            I_start = V_start;
            V_start = V_start + V_rest;

            


            time_first_15mv_IDNI(1,j) =  find( V > -55, 1 ) *  deltaT ;
%             display(time_first_15mv_IDNI(1,j));
%             display(find( V > -55, 1 ) *  deltaT );
            
            
            %% Nominal model

            if j==1
                V_record_nominal = V_record;
                I_record_nominal = I_record;
            end

            % yyaxis right
            % if j>1
            %     yyaxis right
            %     plot(cat(2,t_start,t),cat(2,V_start,V_record(j,:)), '--','LineWidth',1, 'color', '#00008B','HandleVisibility','off');
            %     yyaxis left
            %     plot(cat(2,t_start,t),cat(2,I_start,I_record(j,:)), '--','LineWidth',1, 'color','r','HandleVisibility','off')
            % else
            %     yyaxis right
            %     plot(cat(2,t_start,t),cat(2,V_start,V_record(j,:)),'LineWidth',3, 'color', [0 0.4470 0.7410]);hold on
            %     yyaxis left
            %     plot(cat(2,t_start,t),cat(2,I_start,I_record(j,:)),'LineWidth',3, 'color','r')
            % end
            % hold on
            

            
            
            % Plot the mean values
            % plot(x, mean_values, 'bo-', 'LineWidth', 1.5);
            
            if j==1
                yyaxis right
                plot(cat(2,t_start,t),cat(2,V_start,V_record(j,:)),'LineWidth',3, 'color', [0 0.4470 0.7410]);hold on
                yyaxis left
                plot(cat(2,t_start,t),cat(2,I_start,I_record(j,:)),'LineWidth',3, 'color','r')
            end

    end


    %% Plot
    yyaxis right
    yline(-55,'--b', 'LineWidth',3) % Threshold

    track_plot = track + V_rest;

    set(gca,'FontSize',font_size)
    color_simu = [0 0.4470 0.7410];



    i = 1;
    hold on
    yyaxis right
    % plot(cat(2,t_start,t),cat(2,V_start,V_record_nominal(i,:)),'--','LineWidth',1, 'Color', '#00008B') % uncertain model label

    % add label for uncertain model
    % plot(cat(2,t_start,t),cat(2,V_start,V_record(1,:)), '--','LineWidth',1, 'color', '#00008B','HandleVisibility','off');
    
    % ref traj
    plot(t(1:end-1), track_plot(1,1:end-1), 'color', "#FF00FF",'LineWidth', 2, 'Marker', '_')  %% ref track
    
    %% plot confidence interval for V
    min_values = min(V_record,[],1);
    max_values = max(V_record,[],1);
    
    % Calculate the mean values for each trajectory (optional, for reference)
    mean_values = mean(V_record, 1);
    
    % Calculate the standard error of the mean for each trajectory
    n = size(V_record, 2);  % Number of samples per trajectory
    standard_error = std(V_record, 0, 1) / sqrt(n);
    confidence_interval = 1.96 * standard_error;  % 1.96 is the critical value for a 95% confidence interval
    
    % Create vectors for the ascending and descending x values
    x = 1:size(V_record, 2);  % x-axis represents trajectories
    x = x-1;
    x_ascend = x;
    x_descend = flip(x);
    % Create a loop to fill the confidence intervals for each trajectory
    min_ascend = mean_values - confidence_interval;
    max_descend = flip(mean_values + confidence_interval);
    
    fill([x_ascend, x_descend] * deltaT, [min_ascend, max_descend], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    fill([x_ascend, x_descend] * deltaT, [min_values, flip(max_values)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    %% Continue plotting
    hold on;

    ylabel('Voltage (mV)' ,'FontSize',font_size)
    % ylabel('Current (nA)' ,'FontSize',font_size)
    % xlabel('Time (ms)','FontSize',font_size)



    yyaxis left
    % plot(cat(2,t_start,t),cat(2,I_start,I_record_nominal(i,:)),'LineWidth',1, 'color','r')
    hold on

    %% plot confidence interval for I
    min_values = min(I_record,[],1);
    max_values = max(I_record,[],1);

    % Calculate the mean values for each trajectory (optional, for reference)
    mean_values = mean(I_record, 1);

    % Calculate the standard error of the mean for each trajectory
    n = size(I_record, 2);  % Number of samples per trajectory
    standard_error = std(I_record, 0, 1) / sqrt(n);
    confidence_interval = 1.96 * standard_error;  % 1.96 is the critical value for a 95% confidence interval
    
    % Create vectors for the ascending and descending x values
    x = 1:size(I_record, 2);  % x-axis represents trajectories
    x = x-1;
    x_ascend = x;
    x_descend = flip(x);
    % Create a loop to fill the confidence intervals for each trajectory
    min_ascend = mean_values - confidence_interval;
    max_descend = flip(mean_values + confidence_interval);

    fill([x_ascend, x_descend] * deltaT, [min_ascend, max_descend], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    fill([x_ascend, x_descend] * deltaT, [min_values, flip(max_values)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    %% Continue plotting
    hold on
    % ylabel('Voltage (mV)' ,'FontSize',font_size)
    ylabel('Current (nA)' ,'FontSize',font_size)

    grid on
    % end 

    title('(b) SRDI','FontSize',font_size)
    % legend( 'Reference Trajectory', 'Threshold (-55mV)' , 'Voltage (Uncertain Model)', 'Current Injection','Input Contraint (15 nA)', 'FontSize', font_size)

    h_legend = legend( 'Current (Nominal Model)', 'Current (95% Confidence Region)', 'Current (All Possible Region)', 'Voltage (Nonimal Model)', 'Threshold (-55mV)' , 'Reference Trajectory', 'Voltage (95% Confidence Region)', 'Voltage (All Possible Region)', 'FontSize', font_size);
    set(h_legend, 'location', 'northeastoutside')
    
    % 'Threshold (-55mV)' ,
    hold on





end


% for j = 1:num_traject
%     plot( cat(2,linspace(1, length(t_start), length(t_start)), linspace(1+ length(t_start), length(error_record(j,:)) + length(t_start), length(error_record(j,:)) )    )  , cat(2,zeros(length(t_start),1)', error_record(j,:) )  )
%     hold on
% end
% xlim([180 300])
% xlabel('Time (0.1 nS)')
% ylabel('Current (nA)')
% title('Dynamic Inversion Error','FontSize',font_size)




%%
%%
%% MPC

if(MPC == 1)
    % num_traject = 2;
    subplot(3,3,[7,8])
    Ts = 0.01;

    ref = 15 / (1 - exp(-0.5));

    % font_size = 18;
    % uncertain_percent = 0.2;


    E_Na = 115; %mV
    E_K = -12; %mV
    V_leak = 10.613; %mV
    C=1;
    G_K=36;
    G_Na=120;
    G_leak = 0.3;

    % deltaT=0.01;

    % T_start = 2;
    % Total_time=4;  %10 ms

    t_start=0:deltaT:T_start;
    V_start = zeros(1,numel(t_start));
    t=T_start:deltaT:(Total_time+T_start);

    V =0; % baseline voltage
    alpha_n=0.01*(10-V)/(exp((10-V)/10)-1);
    beta_n=0.125*exp(-V/80);
    alpha_m=0.1*(25-V)/(exp((25-V)/10)-1);
    beta_m=4*exp(-V/18);
    alpha_h=0.07*exp(-V/20);
    beta_h=1/(exp((30-V)/10)+1);
    n(1)=alpha_n/(alpha_n+beta_n);
    m(1)=alpha_m/(alpha_m+beta_m);
    h(1)=alpha_h/(alpha_h+beta_h);
    I=zeros(1,numel(t));

    % I_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));
    % V_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));
    error_record_MPC = zeros();
    error_time = zeros();

    % for i=1:numel(t)
    %     track(i) = ref * (1 - exp( -1 *(i-1)*deltaT));
    % end
    % track(52:end) = track(51);

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

    J_A = jacobian(f_x,x_state);
    J_B = jacobian(g_x,u_e);

    k = 1;

    I_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));
    V_record = zeros(numel(lambda_dynamic_inversion)-1, numel(t));

    %% Uncertain Model
    Model_0 = 1;
    if(Model_0 ==1)
        for j = 1:num_traject
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
            t_start=0:deltaT:T_start;

            t=T_start:deltaT:(Total_time+T_start);

            lambda_dynamic_inversion = 1;
            

            V =1e-10; % baseline voltage

            I=zeros(1,numel(t));
            %% Uncertain model
            if(j==1)
                percent = 0;
            end

            if(j>1)
                percent = uncertain_percent;
            end

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


            for i=1:numel(t)-1

                alpha_n(i) = .01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) );
                beta_n(i) = .125*exp(-V(i)/80);
                alpha_m(i) = .1*( (25-V(i)) / (exp((25-V(i))/10)-1) );
                beta_m(i) = 4*exp(-V(i)/18);
                alpha_h(i) = .07*exp(-V(i)/20);
                beta_h(i) = 1/(exp((30-V(i))/10)+1);

                %% System Linearized at x_0, u_0

                x_0 = [V(i), n(i), m(i), h(i)]';
                u_0 = I(i);


                A_0 = subs(J_A, {V_m, n_d, m_d, h_d,      a_n, a_m, a_h,                 b_n, b_m, b_h}, ...
                    {V(i), n(i), m(i), h(i),      alpha_n(i), alpha_m(i), alpha_h(i),    beta_n(i), beta_m(i),   beta_h(i)});
                B_0 = subs(J_B, {V_m, n_d, m_d, h_d,      a_n, a_m, a_h,     b_n, b_m, b_h,     u_e}, ...
                    {V(i), n(i), m(i), h(i),      alpha_n(i), alpha_m(i), alpha_h(i),   beta_n(i), beta_m(i),   beta_h(i),    u_0});

                f_x_0 = subs(f_x,{V_m, n_d, m_d, h_d,      a_n, a_m, a_h,     b_n, b_m, b_h,     u_e}, ...
                    {V(i), n(i), m(i), h(i),      alpha_n(i), alpha_m(i), alpha_h(i),   beta_n(i), beta_m(i),   beta_h(i),    u_0});

                C_0 = [1,0,0,0];

                sys = ss(eval(A_0), eval(B_0), C_0, 0);
                D_0 = sys.d;

                %% MPC

                p=predict_horizon;

                plant = sys;
%                 plant = c2d(sys,Ts);

                % Q = [1,0,0,0];
                Q = 10 * eye(1);

                R = 1;
                [K,Qp] = lqry(plant,Q,R);
                newPlant = plant;

                d = eig(Qp);
                isposdef = all(d > 1e-4);
                if(isposdef && issymmetric(Qp))
                    L = chol(Qp);
                    newPlant = plant;
                    set(newPlant,'C',[C_0;L],'D',[D_0;zeros(4,1)] );


                    newPlant = setmpcsignals(newPlant,'MO',[1],'UO',[2 3 4 5]);

                    setpoint = cat(2,0,track(1:end));
                    if((i+p < length(track))  )
                        xref = zeros(length(setpoint(i:i+p-1)),5);
                        xref(:,1) = setpoint(i:i+p-1); %setpoint(i+1:i+p)
                        xref(:,2) = zeros(length(setpoint(i:i+p-1)),1); %zeros(length(setpoint(i:i+p-1)),1);
                        xref(:,3) = zeros(length(setpoint(i:i+p-1)),1);
                        xref(:,4) = zeros(length(setpoint(i:i+p-1)),1);
                        xref(:,5) = zeros(length(setpoint(i:i+p-1)),1);
                    end

                    mpcobj = mpc(newPlant,Ts);
                    mpcobj.Weights.MV = 1;
                    mpcobj.Weights.ManipulatedVariables =1;
                    mpcobj.Weights.MVRate = 0.0001;
                    mpcobj.Weights.ManipulatedVariablesRate = 0.0001;  
                    mpcobj.PredictionHorizon = p;
                    mpcobj.PredictionHorizon = p;
                    mpcobj.ControlHorizon = control_horizon;
                    mpcobj.Weights.OutputVariables = [5 0 0 0 0];

                    setoutdist(mpcobj,'model',ss(zeros(5,1))); % disturbance
                    simOpt = mpcsimopt(mpcobj);

                    x0 = [V(i), n(i), m(i), h(i)]';
                    simOpt.PlantInitialState = x0;
            %         Tstop = 25;
                    r = xref;

%                     xc = mpcstate(mpcobj);

                    y = V(i);

                    if((V(i) < 15) && (i*deltaT < 3))
%                         u = mpcmove(mpcobj,xc,y,r);
                        [y,ttt,u] = sim(mpcobj,p,r,simOpt);

                        I(i) = u(2);
                    else
%                         u=0;
                        I(i) = 0;
                    end
                else
                    I(i) = 0;
                end

%                     %% Noise term
%                     mu = 0;
%                     std = 80/100;
                noise=noise_std*randn(1,1)+mu;
% 
%                     noise = 0; %% test

                I(i) = I(i) + noise;

                % Constrained input
                interval = 70 - 55;
                

            %     if((i>1)&& (V(i)<V(i-1)) ||(V(i)> interval) || V(i)<0) %% Positive current
                if( (V(i)> interval) || (i * deltaT) > 3 || V(i)<0 || I(i)<0 ) %% Positive current
            %     if((V(i)> interval) || V(i)<0)  %% Negative current
%                         I(i+1) = 0;
                    I(i) = 0;
                end
                
                if (V(i)< interval && (i * deltaT) < 3) % record error
                    error_record_MPC(j,i) = V(i) - setpoint(i);
                end
                
%                     if (V(i)< interval && time < 3) % record error
%                         error_record_CDI(j,i) = V_th - V(i);
% 
%                     end
                %% Model Dynamics

                I_Na = double((m(i)^3)) *G_Na * double(h(i)) * (V(i)-E_Na); %Equations 3 and 14
                I_K = double((n(i)^4)) * G_K * (V(i)-E_K); %Equations 4 and 6
                I_leak = G_leak * (V(i)-V_leak);
                I_ion = I(i) - I_K - I_Na - I_leak;


                V(i+1) = V(i) + deltaT*I_ion/C;
                n(i+1) = n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i));
                m(i+1) = m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i));
                h(i+1) = h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i));

            end


    %% Record and plot
            V_rest = -70;
            % V_rest = 0;

            V = V+V_rest; %Set resting potential to -70mv
            I_record(j, :) = I;
            V_record(j, :) = V;
            I_start = V_start;
            V_start = V_start + V_rest;

            % time_first_15mv_MPC(1,j) =  find( V > -55, 1 ) *  deltaT ;

            %% Nominal model
            if j==1
                V_record_nominal = V_record;
                I_record_nominal = I_record;
            end

            % yyaxis right
            % if j>1
            %     yyaxis right
            %     plot(cat(2,t_start,t),cat(2,V_start,V_record(j,:)), '--','LineWidth',1, 'color', '#00008B','HandleVisibility','off');
            %     yyaxis left
            %     plot(cat(2,t_start,t),cat(2,I_start,I_record(j,:)), '--','LineWidth',1, 'color','r','HandleVisibility','off')
            % else
            %     yyaxis right
            %     plot(cat(2,t_start,t),cat(2,V_start,V_record(j,:)),'LineWidth',3, 'color', [0 0.4470 0.7410]);hold on
            %     yyaxis left
            %     plot(cat(2,t_start,t),cat(2,I_start,I_record(j,:)),'LineWidth',3, 'color','r')
            % end
            if j==1
                yyaxis right
                plot(cat(2,t_start,t),cat(2,V_start,V_record(j,:)),'LineWidth',3, 'color', [0 0.4470 0.7410]);hold on
                yyaxis left
                plot(cat(2,t_start,t),cat(2,I_start,I_record(j,:)),'LineWidth',3, 'color','r')
            end
            hold on
        end
    end
    yyaxis right
    yline(-55,'--b', 'LineWidth',3) % Threshold

    track_plot = track + V_rest;

    set(gca,'FontSize',font_size)
    color_simu = [0 0.4470 0.7410];



    i = 1;
    hold on
    yyaxis right
    % plot(cat(2,t_start,t),cat(2,V_start,V_record_nominal(i,:)),'--','LineWidth',1, 'Color', '#00008B') % uncertain model label

    % add label for uncertain model
    % plot(cat(2,t_start,t),cat(2,V_start,V_record(1,:)), '--','LineWidth',1, 'color', '#00008B','HandleVisibility','off');
    
    % ref traj
    plot(t(1:end-1), track_plot(1,1:end-1), 'color', "#FF00FF",'LineWidth', 2, 'Marker', '_')  %% ref track
    
    %% plot confidence interval for V
    min_values = min(V_record,[],1);
    max_values = max(V_record,[],1);
    
    % Calculate the mean values for each trajectory (optional, for reference)
    mean_values = mean(V_record, 1);
    
    % Calculate the standard error of the mean for each trajectory
    n = size(V_record, 2);  % Number of samples per trajectory
    standard_error = std(V_record, 0, 1) / sqrt(n);
    confidence_interval = 1.96 * standard_error;  % 1.96 is the critical value for a 95% confidence interval
    
    % Create vectors for the ascending and descending x values
    x = 1:size(V_record, 2);  % x-axis represents trajectories
    x = x-1;
    x_ascend = x;
    x_descend = flip(x);
    % Create a loop to fill the confidence intervals for each trajectory
    min_ascend = mean_values - confidence_interval;
    max_descend = flip(mean_values + confidence_interval);
    
    fill([x_ascend, x_descend] * deltaT, [min_ascend, max_descend], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    fill([x_ascend, x_descend] * deltaT, [min_values, flip(max_values)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    %% Continue plotting
    hold on;

    ylabel('Voltage (mV)' ,'FontSize',font_size)
    % ylabel('Current (nA)' ,'FontSize',font_size)
    xlabel('Time (ms)','FontSize',font_size)



    yyaxis left
    % plot(cat(2,t_start,t),cat(2,I_start,I_record_nominal(i,:)),'LineWidth',1, 'color','r')
    hold on

    %% plot confidence interval for I
    min_values = min(I_record,[],1);
    max_values = max(I_record,[],1);

    % Calculate the mean values for each trajectory (optional, for reference)
    mean_values = mean(I_record, 1);

    % Calculate the standard error of the mean for each trajectory
    n = size(I_record, 2);  % Number of samples per trajectory
    standard_error = std(I_record, 0, 1) / sqrt(n);
    confidence_interval = 1.96 * standard_error;  % 1.96 is the critical value for a 95% confidence interval
    
    % Create vectors for the ascending and descending x values
    x = 1:size(I_record, 2);  % x-axis represents trajectories
    x = x-1;
    x_ascend = x;
    x_descend = flip(x);
    % Create a loop to fill the confidence intervals for each trajectory
    min_ascend = mean_values - confidence_interval;
    max_descend = flip(mean_values + confidence_interval);

    fill([x_ascend, x_descend] * deltaT, [min_ascend, max_descend], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    fill([x_ascend, x_descend] * deltaT, [min_values, flip(max_values)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    %% Continue plotting
    hold on
    % ylabel('Voltage (mV)' ,'FontSize',font_size)
    ylabel('Current (nA)' ,'FontSize',font_size)

    grid on
    % end 

    title('(c) MPC','FontSize',font_size)
    % legend( 'Reference Trajectory', 'Threshold (-55mV)' , 'Voltage (Uncertain Model)', 'Current Injection','Input Contraint (15 nA)', 'FontSize', font_size)

    % h_legend = legend( 'Voltage (Nominal Model)', 'Threshold (-55mV)' , 'Reference Trajectory', 'Voltage (95% Confidence Region)', 'Voltage (All Possible Region)','Current (Nonimal Model)','Current (95% Confidence Region)', 'Current (All Possible Region)', 'FontSize', font_size);
    % set(h_legend, 'location', 'northeastoutside')
    
    % 'Threshold (-55mV)' ,
    hold on

    % %% Reference track
    % set(gca,'FontSize',font_size)
    % color_simu = [0 0.4470 0.7410];
    % 
    % i = 1;
    % hold on
    % yyaxis right
    % plot(t(1:end), track_plot, 'color', "#FF00FF",'LineWidth', 1, 'Marker', '_')  %% ref track
    % 
    % ylabel('Current (nA)','FontSize',font_size)
    % xlabel('Time (ms)','FontSize',font_size)
    % 
    % yyaxis left
    % plot(cat(2,t_start,t),cat(2,I_start,I_record_nominal(i,:)),'LineWidth',1, 'color','r')
    % hold on
    % ylabel('Voltage (mV)' ,'FontSize',font_size)
    % line(nan, nan, 'color', '#00008B', 'linestyle', '--', 'linewidth', 1);
    % line(nan, nan, 'color','r',      'linestyle', '--', 'linewidth', 1);
    % grid on
    % 
    % 
    % title('(c) MPC','FontSize',font_size)
    % h_legend = legend( 'Voltage (Nominal Model)', 'Threshold (-55mV)' , 'Reference Trajectory', 'Current (Nominal Model)','Current (Uncertain Model)', 'Voltage (Uncertain Model)', 'FontSize', font_size);
    % set(h_legend, 'location', 'northeastoutside')
    % set(gca,'FontSize',font_size);
    % hold off
end







%% Plot
Erro_plot = 0;
if(Erro_plot == 1)
    clear noise_std
    
    time_first_15mv_DI = nonzeros(time_first_15mv_DI);
    
    y = {time_first_15mv_DI; time_first_15mv_IDNI; time_first_15mv_MPC};
    
    mean_velocity = zeros();
    for i = 1:3
        mean_velocity(i) = mean(y{i} + 2); % mean velocity
    end
    std_velocity = zeros();
    for i = 1:3
        std_velocity(i) = sqrt(var(y{i}));
    end
    
    colorstring = 'bgry';
    
    
    
    figure
    xData = 1:3;
    bar(xData, mean_velocity,'FaceAlpha',0.5 )
    yline(2.5,'--r','2.5 ms','FontSize',font_size, LineWidth=2);
    hold on 
    errlow = 2 * std_velocity;
    errhigh = 2 * std_velocity;
    er = errorbar(xData,mean_velocity',errlow,errhigh , '.r', "CapSize",30, 'MarkerSize',1, 'MarkerEdgeColor', 'g', 'LineWidth',2,'Markersize',3);    
    % er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    ylim([2.2 2.9])
    % set(gca,'xtick',['CDI', 'SRDI', 'MPC'])
    ylabel("Time (ms)")
    set(gca,'xticklabel',{'CDI';'SRDI'; 'MPC'})
    set(gca,'FontSize',font_size);
    hold off
    
    
    %% Error dyanamic
    font_size = 23;
    if(Dynamic_inversion == 1)
        subplot(1,3,1)
        for j = 1:num_traject
            plot( cat(2,linspace(1, length(t_start), length(t_start)), linspace(1+ length(t_start), length(error_record_CDI(j,:)) + length(t_start), length(error_record_CDI(j,:)) )    )  , cat(2,zeros(length(t_start),1)', abs(error_record_CDI(j,:)) ), 'color', '#00008B'  )
            hold on
        end
        xlim([180 300])
        xlabel('Time (0.05 ms)')
        ylabel('Current (nA)')
        title('Error','FontSize',font_size)
        set(gca,'FontSize',font_size);
        hold on
    end
    
    if(IDNI == 1)
        subplot(1,3,2)
        for j = 1:num_traject
            plot( cat(2,linspace(1, length(t_start), length(t_start)), linspace(1+ length(t_start), length(error_record_SRDI(j,:)) + length(t_start), length(error_record_SRDI(j,:)) )    )  , cat(2,zeros(length(t_start),1)', abs(error_record_SRDI(j,:)) ), 'color', '#00008B' )
            hold on
        end
        xlim([180 260])
        xlabel('Time (0.05 ms)')
        ylabel('Current (nA)')
        title('Error','FontSize',font_size)
    %     xlim()
    set(gca,'FontSize',font_size);
        hold on
    end
    
    if(MPC == 1)
        subplot(1,3,3)
        for j = 1:num_traject
            plot( cat(2,linspace(1, length(t_start), length(t_start)), linspace(1+ length(t_start), length(error_record_MPC(j,:)) + length(t_start), length(error_record_MPC(j,:)) )    )  , cat(2,zeros(length(t_start),1)', abs(error_record_MPC(j,:)) ), 'color', '#00008B' )
            hold on
        end
        xlim([180 300])
        xlabel('Time (0.05 ms)')
        ylabel('Current (nA)')
        title('Error','FontSize',font_size)
        set(gca,'FontSize',font_size);
        hold off
    end
end