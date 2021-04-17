%% Main function simulating and plotting ODE equations
function HCV_main

% define all variables that are used by other scripts as global (for example, here we call Virology_HCV_ode.m) 
global nt;
global bt;
global T_0;
global del;
global ep;
global p;
global c;

% define parameter values
nt=0.5; % represents reduction in infectivity due to antiviral treatment
bt=7*10^-9 ; % represents viral infectivity
T_0=10^7; % represents the constant number/concentration of target cells including at time t=0
del=0.2; % represents the death rate of infected cells 
ep=0.998; % represents efficacy of the drug in inhibiting the viral production
c=23; % representing viral clearance

% Since the system is assumed to be in steady state at t=0 (the time of start of antiviral treatment)
% at which we also are often aware of viral loads (V_0), we additionally have the following constraints 
V_0=10^7 ; % assumed (usually known) copies/mL level in chrnoic HCV patients
I_0=bt*V_0*T_0/del; % from dI/dt equation
p=c*V_0/I_0; % from dV/dt equation


%% Simulation for defined parameters
xinit=[I_0;V_0]; % define initial values of dependent variables in the system of ODEs in the same order as they are described 
tstop=30; % specify how long the simulation needs to be run (here we assume 30 days of treatment)
tspan=[0:tstop]; 

[t,y]=ode23s(@Virology_HCV_ode,tspan,xinit); % calling ode solver

%% Plotting simulation with y-axis on log scale using semilogy
figure(1)
semilogy(t,y(:,1),'Color','black','LineStyle','-','LineWidth',1) % plot infected cell population
hold on
semilogy(t,y(:,2),'Color','red','LineStyle','-','LineWidth',1)  % plot viral loads 
legend('Infected cells','Viral loads')
ylim([10^-2 10^8])
xlim([0 tstop])
ylabel('Concentration of')
xlabel('time (days)')


