% Biopotential simulation of an excitable cell membrane according to Hodgkin-Huxley
% membrane model as described by J Kenyon
clear variables;

%% Specify fix parameters that define the HH simulation
delta_t = 1E-3; % Time step size [ms]
T_tot = 20; % Total simulation time [ms]
noT = T_tot/delta_t; % No of time steps in simulation

CM = 1; % membrane capacitance [uF/ cm^2]
V_Na = -115; % Nernst potential for Na [mV]
V_K = 12; % [mV]
V_L = -10.613; % [mV]
g_Na = 120; % scales Na conductance; [m mho/ cm^2] mho = 1/Ohm
g_K = 20; % scales K conductance; [m mho/ cm^2]
g_L = 0.3; % scales leakage conductance; [m mho/ cm^2]

I_t0 = 1; % Current stimuli timepoint [ms]
I_dt = 0.1; % Current stimuli duration [ms]
I_amp = 500; % Current stimuli amplitud [?A/cm^2]
I_second = 8;

%% Allocate time dependent variables.
%    Most of these variables are stored for post-ploting purposes.
%    alfa*, beta*, n, m and h are HH specific parameter for calculating the
%    conductance G* (see Kenyon for details).

t = (1:noT) * delta_t;  % time vector
Vm = zeros(1, noT);     % Membrane voltage/potential centered around the average resting potential (see Kenyon for details).
I_stim = zeros(1, noT); % Current stimuli
alfa_n = zeros(1, noT);
beta_n = zeros(1, noT);
alfa_m = zeros(1, noT);
beta_m = zeros(1, noT);
alfa_h = zeros(1, noT);
beta_h = zeros(1, noT);
delta_V = zeros(1, noT);    % Change in membrane voltage
n = zeros(1, noT);
m = zeros(1, noT);
h = zeros(1, noT);
G_Na = zeros(1, noT);   % Na conductance
G_K = zeros(1, noT);    % K conductance
G_L = zeros(1, noT);    % Leakage conductance
I_Na = zeros(1, noT);   % Na current
I_K = zeros(1, noT);    % K current
I_L = zeros(1, noT);    % Leakage current
I_C = zeros(1, noT);    % Capacitive current

%% Run simulation using Eulers method
for i = 1 : noT-1
    % Define current stimuli  -  Use this for task 2, 6 and 7.
    if ((t(i)>=I_t0 && t(i)<I_t0+I_dt) || (t(i)>=I_t0+I_second && t(i)<I_t0+I_dt+I_second))
        I_stim(i) = I_amp;
    end
    %% Add a voltage stimuli  -  Use this for task 3, 4 and 5.
%     if t(i)==I_t0
%         Vm(i) = -20;
%     end
    
	alfa_n(i) = 0.01 * (Vm(i) + 10)/ (exp(1 + 0.1 * Vm(i)) - 1);
	beta_n(i) = 0.125 * exp(Vm(i)/ 80);
	alfa_m(i) = 0.1 * (Vm(i) + 25)/(exp(0.1 * Vm(i) + 2.5) - 1);
	beta_m(i) = 4 * exp(Vm(i) / 18);
	alfa_h(i) = 0.07 * exp(Vm(i) / 20);
    beta_h(i) = 1 / (exp(3 + 0.1 * Vm(i)) + 1);
    
    if i==1 
        % Calculate initiate state parameters
        n(1) = alfa_n(1) / (alfa_n(1) + beta_n(1));
        m(1) = alfa_m(1) / (alfa_m(1) + beta_m(1));
        h(1) = alfa_h(1) / (alfa_h(1) + beta_h(1));
    end
    
    G_Na(i) = g_Na * m(i)^3 * h(i);
    I_Na(i) = G_Na(i) * (Vm(i) - V_Na);
    
    G_K(i) = g_K * n(i)^4;
    I_K(i) = G_K(i) * (Vm(i) - V_K);
 
    G_L(i) = g_L;
    I_L(i) = G_L(i) * (Vm(i) - V_L);
    
    %% Add the calculation of the derivative of the four state parameters
    delta_m = alfa_m(i) * (1 - m(i)) - beta_m(i) * m(i);
    delta_n = alfa_n(i) * (1 - n(i)) - beta_n(i) * n(i);
    delta_h = alfa_h(i) * (1 - h(i)) - beta_h(i) * h(i);
    delta_V(i) = -(1/CM)*(I_Na(i) + I_K(i) + I_L(i) + I_stim(i));
    
    %% Add the itterative update of the four state parameters using Eulers method
    m(i + 1) =  m(i) + delta_m * delta_t;
    n(i + 1) =  n(i) + delta_n * delta_t;
    h(i + 1) =  h(i) + delta_h * delta_t;
    Vm(i + 1) =  Vm(i) + delta_V(i) * delta_t;
    
    I_C(i) = -CM * delta_V(i); % Capacitive current
end

figure(101);
% Plot the membrain potential and the stimuli current assuming a resting membrane potential of -60mV.
P = plot(t,-Vm-60,'b',t,I_stim,'r');
set(P,'linewidth',2);
set(gca, 'fontsize', 14);
L1 = legend('V_m','I_{stim}');
xlabel('[ms]');
ylabel('[mV] or [uA/cm^2]');
title(strcat('Time step:',num2str(delta_t)));
% hold on


