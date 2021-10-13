%Function for Non-conservative form of quasi 1-D supersonic flow
%Author: Aditya Iyer Ramesh

function [mass_flow_rate_non_cons, pressure_non_cons, mach_number_non_cons, rho_non_cons, V_non_cons, T_non_cons, rho_throat_non_cons, V_throat_non_cons, T_throat_non_cons, mass_flow_rate_throat_non_cons, pressure_throat_non_cons, mach_number_throat_non_cons] = non_conserv(x,dx,n,nt,gamma,C)

%Initial BCs
rho_non_cons = 1 - 0.31468*x;
T_non_cons = 1 - 0.2341*x;
V_non_cons = (0.1 + 1.09*x).*T_non_cons.^0.5;

%Profiles for Area and Throat Region
a = 1 + 2.2*(x - 1.5).^2;
throat = find(a==1);

%Time-loop (outer)
for k = 1:nt
    
    %control for time-steps
    dt = min(C.*dx./(sqrt(T_non-cons) + V_non_cons));
    
    %A copy of previous values being saved
    rho_old = rho_non_cons;
    T_old = T_non_cons;
    V_old = V_non_cons;
    
    %Predictor Procedure
    for j = 2:n-1
        
        %Simplification of spatial terms
        dV_dx = (V_non_cons(j+1) - V_non_cons(j))/dx;
        drho_dx = (rho_non_cons(j+1) - rho_non_cons(j))/dx;
        dT_dx = (T_non_cons(j+1) - T_non_cons(j))/dx;
        dlog_a_dx = (log(a(j+1)) - log(a(j)))/dx;
        
        %Equation of continuity
        drho_dt_P(j) = -rho(j)*dv_dx - rho_non_cons(j)*V_non_cons(j)*dlog_a_dx - V_non_cons(j)*drho_dx;
        
        %Equation for Momentum
        dV_dt_P(j) = -V_non_cons(j)*dV_dx - (1/gamma)*(dT_dx + (T_non_cons(j)*drho_dx/rho_non_cons(j)));
        
        %Equation for Energy
        dT_dt_P(j) = -V_non_cons(j)*dT_dx - (gamma -1)*T_non_cons(j)*(dV_dx + V-non_cons(j)*dlog_a_dx);
        
        %Updation of solution after every loop
        rho_non_cons(j) = rho_non_cons(j) + drho_dt_P(j)*dt;
        V_non_cons(j) = V_non_cons(j) + dV_dt_P(j)*dt;
        T_non_cons(j) = T_non_cons(j) + dT_dt_P(j)*dt;
        
    end
    
    %Corrector Procedure
    for j = 2:n-1
        
        %Simplification of spatial terms
        dV_dx = (V_non_cons(j-1) - V_non_cons(j))/dx;
        drho_dx = (rho_non_cons(j-1) - rho_non_cons(j))/dx;
        dT_dx = (T_non_cons(j-1) - T_non_cons(j))/dx;
        dlog_a_dx = (log(a(j-1)) - log(a(j)))/dx;
        
        %Equation of Continuity
        drho_dt_C(j) = -rho(j)*dv_dx - rho_non_cons(j)*V_non_cons(j)*dlog_a_dx - V_non_cons(j)*drho_dx;
        
        %Equation for Momentum
        dV_dt_C(j) = -V_non_cons(j)*dV_dx - (1/gamma)*(dT_dx + (T_non_cons(j)*drho_dx/rho_non_cons(j)));
        
        %Equation for Energy
        dT_dt_C(j) = -V_non_cons(j)*dT_dx - (gamma -1)*T_non_cons(j)*(dV_dx + V-non_cons(j)*dlog_a_dx);
        
    end
    
    %Computation for Time-derivative (Average)
    drho_dt = 0.5.*(drho_dt_P + drho_dt_C);
    dV_dt = 0.5.*(dV_dt_P + dV_dt_C);
    dT_dt = 0.5.*(dT_dt_P + dT_dt_C);
    
    %Updating the Final Solution
    for i = 2:n-1
        rho_non_cons(i) = rho_old(i) + drho_dt(i)*dt;
        V_non_cons(i) = V_old(i) + dV_dt(i)*dt;
        T_non_cons(i) = T_old(i) + dT_dt(i)*dt;
        
    end
    
    %Inlet BCs
    rho_non_cons(1) = 1;
    V_non_cons(1) = 2*V_non_cons(2) - V_non_cons(3);
    T_non_cons(1) = 1;
    
    %Outlet BCs
    rho_non_cons(n) = 2*rho_non_cons(n-1) - rho_non_cons(n-2);
    V_non_cons(n) = 2*V_non_cons(n-1) - V_non_cons(n-2);
    T_non_cons(n) = 2*T_non_cons(n-1) - T_non_cons(n-2);
    
    %Equations for Mass Flow Rate, Mach Number, Pressure (Definition)
    mass_flow_rate_non_cons = rho_non_cons.*a.*V_non_cons;
    pressure_non_cons = rho_non_cons.*T_non_cons;
    mach_number_non_cons = (V_non_cons./sqrt(T_non_cons));
    
    %Throat section - Calculation
    rho_throat_non_cons(k) = rho_non_cons(throat);
    T_throat_non_cons(k) = T_non_cons(throat);
    V_throat_non_cons(k) = V_non_cons(throat);
    mass_flow_rate_throat_non_cons(k) = mass_flow_rate_non_cons(throat);
    mach_number_throat_non_cons(k) = mach_number_non_cons(throat);
    pressure_throat_non_cons(k) = pressure_non_cons(throat);
    
    %Graph Plots for juxtaposing non-dimensional mass flow rates at
    %different time-steps
    figure(6);
    if k == 50
        plot(x,mass_flow_rate_non_cons,'y','linewidth',1.25);
        hold on
    elseif k == 100
        plot(x,mass_flow_rate_non_cons,'g','linewidth',1.25)
        hold on
    elseif k == 200
        plot(x,mass_flow_rate_non_cons,'m','linewidth',1.25)
        hold on
    elseif k == 400
        plot(x,mass_flow_rate_non_cons,'c','linewidth',1.25)
        hold on
    elseif k == 800
        plot(x,mass_flow_rate_non_cons,'k','linewidth',1.25)
        hold on
    end
    
    %Inscribing titles
    title("Mass Flow Rates at Diff. Timesteps - Non Conservative Form")
    xlabel("Dist.")
    ylabel("Mass Flow Rate")
    legend(["50^t^h Timestep","100^t^h Timestep","200^t^h Timestep","400^t^h Timestep","800^t^h Timestep"])
    axis([0 3 0 2])
    grid on
end

end 