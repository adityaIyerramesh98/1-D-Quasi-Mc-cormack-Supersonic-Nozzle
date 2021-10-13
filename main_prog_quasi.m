%Main program for solving the Quasi 1-D Supersonic Nozzle Flow using
%Mac-cormack's method
%Author: Aditya Iyer Ramesh

close all
clc

%Properties wto be applied as Input
nt = 1400;   %Number of timesteps

n = 31; %Number of grid points

gamma = 1.4;    %sp. heat capacity  ratio
x = linspace(0,3,n);   %range for initial points
dx = x(2) - x(1);      %numerical derivative

C = 0.5;    %Defining Courant Number


%Declaration for Function of Non-Conservative form
tic;
[mass_flow_rate_non_cons, pressure_non_cons, mach_number_non_cons, rho_non_cons, V_non_cons, T_non_cons, rho_throat_non_cons, V_throat_non_cons, T_throat_non_cons, mass_flow_rate_throat_non_cons, pressure_throat_non_cons, mach_number_throat_non_cons] = non_conservative_form(x,dx,n,nt,gamma,C);
Elapsed_Time_non_cons = toc;
fprintf("Elapsed Time for Non-Conservative Form: 0.4%f", Elapsed_Time_non_cons);
hold on

%Graph Plots for Various Parameters
%1) Variation of Timestep properties at the area w.r.t "Throat" in
%Non-conservative form
figure(1)

%TEMPERATURE
subplot(4,1,1)
plot(linspace(1,nt,nt),T_throat_non_cons,"color","r")
ylabel("TEMPERATURE")
legend("Temp at Throat");
axis([0 1400 0.6 1])
grid minor;
title("Variation of Timestep at Area w.r.t Throat - Non Conservative")

%DENSITY
subplot(4,1,2)
plot(linspace(1,nt,nt),rho_throat_non_cons,"color","b")
ylabel("DENSITY")
legend("Dens at Throat");
axis([0 1400 0 1.3])
grid minor;

%MACH NUMBER
subplot(4,1,4)
plot(linspace(1,nt,nt),mach_number_throat_non_cons,"color","g")
ylabel("MACH NUMBER")
legend("Mach No. at Throat");
axis([0 1400 0.6 1.4])
grid minor;

%PRESSURE
subplot(4,1,3)
plot(linspace(1,nt,nt),pressure_throat_non_cons,"color","m")
ylabel("PRESSURE")
legend("Press at Throat");
axis([0 1400 0.4 1.1])
grid minor;

%2) Simulation for SS values of prime significance in Non-conservative Form
figure(2)

%PRESSURE
subplot(4,1,1)
plot(x,pressure_non_cons,"color","g")
ylabel("PRESSURE")
legend("Pressure");
axis([0 3 0 1])
grid minor;
title("Flow-rate distribution variation - Non conservative form")

%MACH NUMBER
subplot(4,1,2)
plot(x,mach_number_non_cons,"color","c")
ylabel("MACH NUMBER")
legend("Mach Number");
axis([0 3 0 4])
grid minor;

%DENSITY
subplot(4,1,3)
plot(x,rho_non_cons,"color","y")
ylabel("DENSITY")
legend("Density");
axis([0 3 0 1])
grid minor;

%TEMPERATURE
subplot(4,1,4)
plot(x,T_non_cons,"color","m")
ylabel("TEMPERATURE")
legend("Temperature");
axis([0 3 0 1])
grid minor;

%Declaration for Function of Conservative form
tic;
[mass_flow_rate_cons, pressure_cons, mach_number_cons, rho_cons, V_cons, T_cons, rho_throat_cons, V_throat_cons, T_throat_cons, mass_flow_rate_throat_cons, pressure_throat_cons, mach_number_throat_cons] = conservative_form(x,dx,n,nt,gamma,C);
Elapsed_Time_cons = toc;
fprintf("Elapsed Time for Conservative Form: 0.4%f", Elapsed_Time_cons);
hold on

%Graph Plots for Various Parameters
%3) Variation of Timestep properties at the area w.r.t "Throat" in
%Conservative form
figure(3)

%TEMPERATURE
subplot(4,1,1)
plot(T_throat_cons,"color","r")
ylabel("TEMPERATURE")
legend("Temp at Throat");
axis([0 1400 0.6 1])
grid minor;
title("Variation of Timestep at Area w.r.t Throat - Conservative")

%DENSITY
subplot(4,1,2)
plot(rho_throat_cons,"color","b")
ylabel("DENSITY")
legend("Dens at Throat");
axis([0 1400 0 1])
grid minor;

%MACH NUMBER
subplot(4,1,4)
plot(mach_number_throat_cons,"color","g")
ylabel("MACH NUMBER")
legend("Mach No. at Throat");
axis([0 1400 0.6 1.4])
grid minor;

%PRESSURE
subplot(4,1,3)
plot(pressure_throat_non_cons,"color","m")
ylabel("PRESSURE")
legend("Press at Throat");
axis([0 1400 0.2 0.8])
grid minor;

%4) Simulation for SS values of prime significance in Conservative Form
figure(4)

%PRESSURE
subplot(4,1,1)
plot(x,pressure_cons,"color","g")
ylabel("PRESSURE")
legend("Pressure");
axis([0 3 0 1])
grid minor;
title("Flow-rate distribution variation - Conservative form")

%MACH NUMBER
subplot(4,1,2)
plot(x,mach_number_cons,"color","c")
ylabel("MACH NUMBER")
legend("Mach Number");
axis([0 3 0 4])
grid minor;

%DENSITY
subplot(4,1,3)
plot(x,rho_cons,"color","y")
ylabel("DENSITY")
legend("Density");
axis([0 3 0 1])
grid minor;

%TEMPERATURE
subplot(4,1,4)
plot(x,T_cons,"color","m")
ylabel("TEMPERATURE")
legend("Temperature");
axis([0 3 0 1])
grid minor;

%Plot of Comparison for MASS FLOW RATES (Normalized) for both forms
figure(5)

%MASS FLOW RATE
hold on
plot(x,mass_flow_rate_non_cons,"color","b",'linewidth','1.5')
hold on
plot(x,mass_flow_rate_cons,"color","r","linewidth","1.5")
hold on
grid on

legend("Non-conservative VS Conservative")
xlabel("domain-length")
ylabel("mass-flow-rate")
title("differentiating norm. mass flow rates btw both forms") 

