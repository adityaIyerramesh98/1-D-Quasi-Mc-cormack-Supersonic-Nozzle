%Function for Conservative form of quasi 1-D supersonic flow
%Author: Aditya Iyer Ramesh

function [mass_flow_rate_cons, pressure_cons, mach_number_cons, rho_cons, V_cons, T_cons, rho_throat_cons, V_throat_cons, T_throat_cons, mass_flow_rate_throat_cons, pressure_throat_cons, mach_number_throat_cons] = conserv_form(x,dx,n,nt,gamma,C);

%Profile for Area and Throat Region
a = 1 + 2.2*(x - 1.5).^2;
throat = find(a==1);

%Profiles at the start for Density and Temperature
for i = 1:n
    
    if (x(i) >= 0 && x(i) <= 0.5)
        rho_cons(i) = 1;
        T_cons(i) = 1;
    elseif (x(i) >= 0.5 && x(i) <= 1.5)
        rho_cons(i) = 1 - 0.366*(x(i) - 0.5);
        T_cons(i) = 1 - 0.167*(x(i) - 0.5);
    elseif (x(i) >= 1.5 && x(i) <= 3)
        rho_cons = 0.634 - 0.3879*(x(i) - 1.5);
        T_cons = 0.833 - 0.3507*(x(i) - 1.5);
    end
end

%Profiles at start for Velocity and Pressure
V_cons = 0.59./(rho_cons.*a);
pressure_cons = (rho_cons.*T_cons);

%Determining Solution vectors for inclusion of constants
U1 = rho_cons.*a;
U2 = rho_cons.*a.*V_cons;
U3 = rho_cons.*a.*(T_cons./(gamma-1)) + ((gamma/2).*V_cons.^2);

%Time-loop (1): Outer Loop
for k = 1 : nt

    %Controlling of time step by inserting CFL equation
    dt = min(C.*dx./(sqrt(T_cons) + V_cons));
    
    %Storing old values of defined solution vectors
    U1_old = U1;
    U2_old = U2;
    U3_old = U3;
    
    %For Flux Vectors, naming them as F1, F2 and F3
    F1 = U2;
    F2 = ((U2.^2)./U1) + ((gamma - 1)/gamma)*(U3 - (gamma/2)*((U2.^2)./U1));
    F3 = ((gamma*U2.*U3)./U1) - (gamma*(gamma - 1)*((U2.^3)./(2*U1.^2)));
    
    %Procedure for Predictor method as done in non-conservative form
    for j = 2:n-1
        
        %Determining the term J (for predictor)
        J2(j) = (1/gamma)*rho_cons(j)*T_cons(j)*((a(j+1) - a(j))/dx);
        
        %Eqn of Continuity
        dU1_dt_P(j) = -((F1(j+1) - F1(j))/dx);
        
        %Momentum Eqn
        dU2_dt_P(j) = -((F2(j+1) - F2(j))/dx) + J2(j);
        
        %Energy Eqn
        dU3_dt_P(j) = -((F3(j+1) - F3(j))/dx);
        
        %Storing New Values
        U1(j) = U1(j) + dU1_dt_P(j)*dt;
        U2(j) = U2(j) + dU2_dt_P(j)*dt;
        U3(j) = U3(j) + dU3_dt_P(j)*dt;
        
        %Storing New Flux Vectors
        F1 = U2;
        F2 = ((U2.^2)./U1) + ((gamma - 1)/gamma)*(U3 - ((gamma*0.5*(U2.^2)./U1)));
        F3 = ((gamma*U2.*U3)./U1) - ((gamma*0.5*(gamma - 1)*(U2./3))./(U1.^2));
    end
    
    %Procedure for Corrector Method as done in Non-conservative form
        for j = 2:n-1
            
            %Determining the term J (for corrector)
            J2(j) = (1/gamma)*rho_cons(j)*T_cons(j)*((a(j) - a(j-1))/dx);
            
            %Classifying J's spatial terms in detail
            dU1_dt_C(j) = -((F1(j) - F1(j-1))/dx);
            dU2_dt_C(j) = -((F2(j) - F2(j-1))/dx) + J2(j);
            dU3_dt_C(j) = -((F3(j) - F3(j-1))/dx);
            
        end
        
        %Calculating the Avg Time Derivative
        dU1_dt = 0.5*(dU1_dt_P + dU1_dt_C);
        dU2_dt = 0.5*(dU2_dt_P + dU2_dt_C);
        dU3_dt = 0.5*(dU3_dt_P + dU3_dt_C);
        
        %Updating our Final Solution
        for i = 2:n-1
            U1(i) = U1_old(i) + dU1_dt(i)*dt;
            U2(i) = U2_old(i) + dU2_dt(i)*dt;
            U3(i) = U3_old(i) + dU3_dt(i)*dt;
        end
        
        %Declaring BCs
        %Inlet BCs
        U1(1) = rho_cons(1)*a(1);
        U2(1) = 2*U2(2) - U2(3);
        U3(1) = U1(1)*((T_cons(1)/(gamma - 1)) + ((gamma/2)*(V_cons(1)).^2));
        V_cons = U2(1)./U1(1);
        
        %Outlet BCs
        U1(n) = 2*U1(n-1) - U1(n-2);
        U2(n) = 2*U2(n-1) - U2(n-2);
        U3(n) = 2*U3(n-1) - U3(n-2);
        
        %Updating first hand parameters
        rho_cons = U1./a;
        V_cons = U2./U1;
        T_cons = (gamma - 1)*((U3./U1 - (gamma/2)*(V_cons).^2));
        
        %Determining Pressure, Mass Flow Rate and Mach Number
        mass_flow_rate_cons = rho_cons.*a.*V_cons;
        pressure_cons = rho_cons.*T_cons;
        mach_number_cons = (V_cons./sqrt(T_cons));
        
        %Finally, determining variables at throat
        rho_throat_cons(k) = rho_cons(throat);
        T_throat_cons(k) = T_cons(throat);
        V_throat_cons(k) = V_cons(throat);
        mass_flow_rate_throat_cons(k) = mass_flow_rate_cons(throat);
        pressure_throat_cons(k) = pressure_cons(throat);
        mach_number_throat_cons(k) = mach_number_cons(throat);
        
        %Graph Plots for Comparison of Non-Dimensional Mass Flow Rates at
        %Different Timesteps
        figure(7)
        if k == 50
            plot(x, mass_flow_rate_cons,"r","linewidth",1.5);
            hold on
        elseif k == 100
            plot(x, mass_flow_rate_cons,"m","linewidth",1.5);
            hold on
        elseif k == 200
            plot(x, mass_flow_rate_cons,"c","linewidth",1.5);
            hold on
        elseif k == 400
            plot(x, mass_flow_rate_cons,"y","linewidth",1.5);
            hold on
        elseif k == 800
            plot(x, mass_flow_rate_cons,"k","linewidth",1.5);
            hold on
        end

        %Allocating appropriate labels for the above plots
        title("Mass Flow Rates for Diff Timesteps - Conservative form")
        xlabel("Domain Length")
        ylabel("Mass Flow Rates")
        legend({'50^t^h Timestep','100^t^h Timestep','200^t^h Timestep','400^t^h Timestep','800^t^h Timestep'});
        axis([0 3 0.4 1])
        grid on
end

end
