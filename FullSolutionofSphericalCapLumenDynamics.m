%%% Full solution of spherical cap lumen dynamics



%% Assign values for model parameters
Lambda=3e-9;  % Reference experimental value: 3e−9~3e−6 m/s
Lambdaw=1e-13;  % Reference experimental value: 1e−14~1e−11 m^3/(N∙s)
Jp=6e16;  % Reference experimental value: 1e17~1e18 m^-2∙s^-1
Jw=5e-10;  % Reference experimental value: 3e-10~20e-10 m∙s^-1
jp=1e6;
jpw=1.5e-19;
Sigmar=-5e-11;
Gamma=0.5e-4;  % Reference experimental value: 10^−4 N/m
Gammaj=0.4e-4;  % Reference experimental value: 10^−4 N/m
KBT=4.1e-21;
L=1.2e-5;  % cell diameter
C_cell=7.3e25;  % cell ion concentration



%% Solve DAEs
tlength=1000;
t_solve=linspace(0,120,tlength);
t_step=t_solve(2)-t_solve(1);

R_solve=zeros(tlength,1);
theta_solve=zeros(tlength,1);
TD_solve=zeros(tlength,1);
LD_solve=zeros(tlength,1);
deltaC_solve=zeros(tlength,1);
N_solve=zeros(tlength,1);
DeltaP_solve=zeros(tlength,1);
lumen_volume_solve=zeros(tlength,1);

step=5000;  % grid size
% initial condition assignment
R_solve(t_start*tlength/1000-1,i_solve)=2e-6;
theta_solve(t_start*tlength/1000-1,i_solve)=1;
N_solve(t_start*tlength/1000-1,i_solve)=1e10;
DeltaP_solve(t_start*tlength/1000-1,i_solve)=100;
lumen_volume_solve(t_start*tlength/1000-1,i_solve)=4/3*pi*(2e-6)^3*(1-3/2*cos(1)+1/2*cos(1)^3);

for i=t_start*tlength/1000:tlength
    R_min=1e-8;R_max=4e-5;R_step=(R_max-R_min)/step;

    error_best=10000;
    for R_search=R_min:R_step:R_max
        RT=R_solve(i-1)*sin(theta_solve(i-1));
        deltaC_solve(i-1)=N_solve(i-1)/lumen_volume_solve(i-1)-C_cell;
        dNdt=4*pi*R_solve(i-1)^2*(1-cos(theta_solve(i-1)))*(-Lambda*deltaC_solve(i-1)+Jp)-2*pi*RT*jp/(L-R_solve(i-1)*sin(theta_solve(i-1)));
        N_solve(i)=N_solve(i-1)+dNdt*t_step;
        deltaPi_next=2*KBT*deltaC_solve(i-1);
        DeltaP_solve(i-1)=Gamma*R_solve(i-1)/R_interp(i-1)^2;
        dvdt=4*pi*R_solve(i-1)^2*(1-cos(theta_solve(i-1)))*(Lambdaw*(deltaPi_next-DeltaP_solve(i-1))+Jw)-2*pi*RT*jpw/(L-R_solve(i-1)*sin(theta_solve(i-1)));
        lumen_volume_solve(i)=lumen_volume_solve(i-1)+dvdt*t_step;

        theta_search=acos((Gammaj-Sigmar)/(DeltaP_solve(i-1)*R_search));
        error=abs(lumen_volume_solve(i)-(2/3*pi*R_search^3*(1-cos(theta_search))^2*(2+cos(theta_search))));

        if error<error_best
            R_solve(i)=R_search;
            theta_solve(i)=theta_search;
            error_best=error;
        end
    end

    %%% calculate TD and LD
    TD_solve(i)=2*R_solve(i)*sin(theta_solve(i));
    LD_solve(i)=2*R_solve(i)*(1-cos(theta_solve(i)));
end

