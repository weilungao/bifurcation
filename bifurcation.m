%%%%bifurcation model%%%%
%model described in Salter et al. 2017 [JGR-ES] (SPV)
%it builds on the model of Bolla Pitta et al. 2003 [WRR] (BRT)
%this model runs under the settings of SPV by default
%but can run with the settings of BRT

makemovie=0; %set to 1 to generate a movie (see supplemental of SPV)

%% set parameter

%%%set governing dimensionless parameters
taus_a=1; %dimenionless shields stress of branch a
asprat=100; %aspect ratio B_a/h_a
S_a=.001; %slope of branch a
Lstar_b=450; %L_b/B_a, branch length excludes upstream cell
Lstar_c=450; 
alpha=3;%parameter for how far upstream morphology adjusts
r=.5;%parameter in slope term of cross stream sed flux feedback
tilt=0; %sets tilt rate, specified as fraction of aggradation rate A, e.g. .5 = 50% of A,0 = no tilt (see SPV)

%%%switches - activate different model settings 
reset_topo=1;%reset current topography? 1=yes, 0=no, 2=yes, start from equilibrium profile
transport_form=2;%sed xport formula, 1=MPM,2=EH,3=custom (BRT uses 1)
down_bc=2 ;%set downstream b.c. 1=equal water surface elevation (BRT), 2=sed flux (SPV)
variable_chezy=0;%0 (SPV) is fixed chezy, 1 to have chezy coefficent which varies with branch entrance depth (BRT)
sealevel_mode=1;%1 to have time varying sea level H_sea, else set to 0
Habc_set=1;%if 0, H_a is used instead of Habc (SPV), set to 1 for Habc (BRT)

%%%set other parameters
B_a=100;%width upstream channel (m)
B_b=B_a*.5;%width of downstream channels, B_i=B_a*.5 (SPV)
B_c=B_a*.5;
lamp=.4;%bed porosity
g=9.81;%gravity 
ps=2650;%particle density kg/m^3
p=1000;%density of water
Rr = (ps-p)/p;%submerged specific gravity of sediment
Cz_a=15;%upstream dimensionless chezy, value is overwritten if using variable Cz (BRT)
timeyear = 31557600;%seconds in a year
L_b=B_a*Lstar_b; %total length of reach for branch b excluding first cell 
L_c=B_a*Lstar_c; 

%%%set numerical parameters
dt_years=.00005;%timestep in years
dt = dt_years * timeyear;%dt in seconds
dx_up=alpha*B_a;%length of first cell - fundamental model parameter
dx=alpha*B_a*.5;%length of downstream cells, SPV uses dx=alpha*B_a*.5 (half of dx_up)
N_b=1+L_b/dx;%branch b # of cells incl. first cell, must be integer
N_c=1+L_c/dx;
timesteps=2*10^6;%number of timesteps
finaltime=timesteps*dt_years;%year model ends
perturb_amp=1*B_a*10^-8;% initial "seed" transverse bed elevation perturbation, (m)
timestep_perturb=10000;%timestep at which perturbation is applied
watlev_accuracy=1*10^-12;%how accurately you iteratively solve for the constancy of water surface elevation (in meters)
recovery=.999;%for recovering from complete avulsion within while loop, see SPV for explanation

%%%set time vector
if reset_topo==0 %initial start time equals end time of previous simulation
    t=linspace(t(j),t(j)+finaltime,timesteps); %time vector (years)
else %initial start time is year 0
    t=linspace(0,finaltime,timesteps);
end


%% set switch-specific parameters

%%%sediment transport formula
    if transport_form==1 %MPM
        n=8;%coefficient in transport formula
        nb=n*ones(1,N_b+1);%specifies n as a vector, for MPM n is spatially and temporally constant
        nc=n*ones(1,N_b+1);
        m=1.5;%exponent in transport formula
        tausc=.047;%critical stress in transport formula
    elseif transport_form==2 %EH
        n=.05*Cz_a^2;
        nb=n*ones(1,N_b+1);%initialize vector of n values, will change within time loop as Cz_b changes. Not spatially constant, because Cz_b can vary with x.
        nc=n*ones(1,N_b+1);
        m=2.5;
        tausc=0;
    elseif transport_form==3 %custom 
        n=.05*Cz_a^2;
        nb=n*ones(1,N_b+1); %here, n is spatially and temporally constant
        nc=n*ones(1,N_b+1);
        m=5;
        tausc=.1;
    end 
    
%%%downstream boundary condition
    if down_bc==2 %B.C. used by SPV
            F=1; %sets bypass fraction, 0 = no flux B.C., 1 = pure bypass (fixed downstream bed elevations)
    elseif down_bc==1
        eta_downstream_average=0;%final cells of branches b and c average to this for entire simulation (ensures bypass)
    end
    
%%%sea level
if sealevel_mode==0 %sea level option toggled off
    sealevel_b=-10000*ones(1,timesteps);%set below minimum bed elevation, if sealevel is < bed elevation it has no effect
    sealevel_c=-10000*ones(1,timesteps);
elseif sealevel_mode==1 %sea level option toggled on
    sea_rate_b=0;%RSL rise rate in m/year (branch b)
    sea_rate_c=0;
    sea_initial_b=5;%initial sea level in meters (branch b)
    sea_initial_c=5;
    if reset_topo==0 %sets sea level to what it was at end of previous run
        sea_initial_b=sealevel_b(j);
        sea_initial_c=sealevel_c(j);
    end

    sealevel_b=linspace(sea_initial_b,sea_rate_b*finaltime+sea_initial_b,timesteps);
    sealevel_c=linspace(sea_initial_c,sea_rate_c*finaltime+sea_initial_c,timesteps);
end 

if variable_chezy==1;
    minchezy=.01; %minimum allowed value of the chezy coefficient.
end


    
%% model initiation

    %%%calculated values
    H_a=B_a/asprat; %upstream channel depth (m)
    ds=H_a*S_a/(Rr*taus_a); %grain size (m)
    
    if variable_chezy==1
        Cz_a=6+2.5*log(H_a/(2.5*ds));% from log law of the wall
        if reset_topo==0 %allows for lengthening/shortening of the branches without resetting Cz vector
            if length(Cz_b)<N_b+1
                Cz_b=horzcat(Cz_b,Cz_b(length(Cz_b)-1)*ones(1,N_b+1-length(Cz_b)));
            elseif length(Cz_b)>N_b+1
                Cz_b=Cz_b(1:N_b);
            end
            if length(Cz_c)<N_c+1
                Cz_c=horzcat(Cz_c,Cz_c(length(Cz_c)-1)*ones(1,N_c+1-length(Cz_c)));
            elseif length(Cz_c)>N_c+1
                Cz_c=Cz_c(1:N_c);
            end
        end
    end
    if reset_topo~=0
    Cz_b=Cz_a*ones(1,N_b+1);
    Cz_c=Cz_a*ones(1,N_c+1); 
    end
    
    Q_a=sqrt(((H_a^2)*((B_a^2)*(Cz_a^2)*taus_a*(ps-p)*g*ds))/p); %from normal flow assumption
    Qs_a=(n*(taus_a-tausc).^m)*((Rr * g * ds) ^ 0.5) * ds*B_a; %using sed xport formula set above
    Q_b=(Q_a*B_b/(B_b+B_c)); %water discharge introduced uniformly between two upstream cells
    Q_c=(Q_a*B_c/(B_b+B_c));
    if reset_topo==0 %sets starting water discharge to previous value if restarting simulation
        Q_b=Qb_save(length(j-1));
        Q_c=Qc_save(length(j-1));
    end
    

%%%create vectors of x values
xb=zeros(1,N_b);%x vector for branch b, values correspond to position of cell center (meters)
xb(1)=dx_up/2;
xb(2)=xb(1)+(dx+dx_up)/2;
for i=3:length(xb)
    xb(i)=xb(i-1)+dx;
end
xc=zeros(1,N_c);
xc(1)=dx_up/2;
xc(2)=xc(1)+(dx+dx_up)/2;
for i=3:length(xc)
    xc(i)=xc(i-1)+dx;
end
xb_face=zeros(1,N_b+1);%set values of x in branch b on cell faces
xc_face=zeros(1,N_c+1);
xb_face(1)=0;
xc_face(1)=0;
xb_face(2)=dx_up;
xc_face(2)=dx_up;
for i=3:length(xb_face)
    xb_face(i)=xb_face(i-1)+dx;
end
for i=3:length(xc_face)
    xc_face(i)=xc_face(i-1)+dx;
end
    
%%%initialize topography
Abran=B_b*L_b+B_c*L_c;%total area of each branch
Atot=Abran+B_b*dx_up+B_c*dx_up;%total area including upstream cells
Qsi0=(Qs_a*Abran/Atot)/(2*(1-F)+2*(Abran/Atot)*F);%equilibrium sed flux at branch entrances    

if reset_topo==1       
    %set initial condition to be a triangular wedge on a flat basement
    S_initialb=.001;%initial slope of triangular wedge
    bedheightb=4;%maximum bed height of triangular wedge (upstream end)
    minetab=0; %basement elevation
    S_initialc=.001;
    bedheightc=4;
    minetac=0; 

    eta_b=max(bedheightb-S_initialb*xb(1:N_b),minetab);
    eta_c=max(bedheightc-S_initialc*xc(1:N_c),minetac);

elseif reset_topo==2 
%reset_topo==2: start from an equilibrium symmetric bifurcation
%requires symmetric bifurcation, i.e. no tilt or unequal B or L, also does not account for Cz variation
    if down_bc==2

            Arun=zeros(1,N_b+1);%upstream area as a function of x 
            Arun(2)=Arun(1)+2*B_b*dx_up;
            for i=3:N_b+1
            Arun(i)=Arun(i-1)+2*B_b*dx;
            end
            Qsi_eq=(F*Qsi0-(Qs_a/2)).*(Arun)./(Atot)+(Qs_a/2); %equilibrium sediment flux vector as a function of x

    elseif down_bc==1 %bypass so sed flux profile is constant (BRT)
        Qsi_eq=Qs_a/2*ones(1,length(xb_face));
    end
    eta_b=zeros(1,N_b);%initialize vector of bed elevations
    eta_c=zeros(1,N_c);

    %compute the slope between cells using the computed(3:fc+1)
    %equilibrium sed flux profile and chosen sed xport formula
    Sb_profile=(((Qsi_eq./(n*B_b*sqrt(Rr*g*(ds^3)))).^(1/m)+tausc).*(Rr*(B_b^(2/3)).*(Cz_b(2)^(2/3))*(g^(1/3))*ds/((Q_a/2)^(2/3)))).^(3/2);
    Sc_profile=(((Qsi_eq./(n*B_c*sqrt(Rr*g*(ds^3)))).^(1/m)+tausc).*(Rr*(B_c^(2/3)).*(Cz_c(2)^(2/3))*(g^(1/3))*ds/((Q_a/2)^(2/3)))).^(3/2);

    %from slope values computed above, calculate the equilibrium 
    %bed elevation profile
    eta_b(2)=eta_b(1)-Sb_profile(2)*(dx_up+dx)/2;
    eta_c(2)=eta_c(1)-Sc_profile(2)*(dx_up+dx)/2;
    for i=3:N_b
    eta_b(i)=eta_b(i-1)-Sb_profile(i)*dx;
    eta_c(i)=eta_c(i-1)-Sc_profile(i)*dx;
    end
    eta_b=eta_b-eta_b(N_b);
    eta_c=eta_c-eta_c(N_c);

elseif reset_topo==0 %use bed elevation profile from end of a previous run
    
    %this section shortens or lengthens profile if L (and therefore N) is 
    %changed from previous run. Works only if dx isn't changed
    if length(eta_b)~=N_b %if eta from previous run has length other than N
        if length(eta_b)<N_b %if new branch length is longer than old
            %linearly extrapolate bed elevations from upstream to downstream 
            %using the old profile, i.e. constant slope in new downstream cells
            Nprime_b=N_b-length(eta_b);
            Sprime_b=(eta_b(N_b-Nprime_b-1)-eta_b(N_b-Nprime_b))/dx;
            append_b=0;%reset previous append2
            append_b(1:Nprime_b)=eta_b(N_b-Nprime_b)-dx*Sprime_b*(1:Nprime_b);
            eta_b=horzcat(eta_b,append_b);
        elseif length(eta_b)>N_b% if new branch length is shorter than old
            eta_b=eta_b(1:N_b); %remove excess downstream cells
        end
    end
    
    %likewise for the other branch
    if length(eta_c)~=N_c
        if length(eta_c)<N_c
            Nprime_c=N_c-length(eta_c);
            Sprime_c=(eta_c(N_c-Nprime_c-1)-eta_c(N_c-Nprime_c))/dx;
            append_c=0;
            append_c(1:Nprime_c)=eta_c(N_c-Nprime_c)-dx*Sprime_c*(1:Nprime_c);
            eta_c=horzcat(eta_c,append_c);
        elseif length(eta_c)>N_c

        eta_c=eta_c(1:N_c);
        end

    end
end

%%%initialize model vectors
H_b=H_a*ones(1,N_b+1); %depth in branch b
H_c=H_a*ones(1,N_c+1);
taus_b=taus_a*ones(1,N_b+1); %shields stress branch b
taus_c=taus_a*ones(1,N_c+1);
Qs_b=(Qs_a*B_b/(B_b+B_c))*ones(1,N_b+1); %sed flux branch b
Qs_c=(Qs_a*B_c/(B_b+B_c))*ones(1,N_c+1); %sed flux branch b

%%%initialize data vectors

%data collected every timestep
Qb_save=zeros(1,timesteps);%water discharge
Qc_save=zeros(1,timesteps);
Qsb_save=zeros(1,timesteps);%sediment flux at branch inlet
Qsc_save=zeros(1,timesteps); 
Qsb_down_save=zeros(1,timesteps); %sediment flux at downstream boundary
Qsc_down_save=zeros(1,timesteps); 
Qsy_save=zeros(1,timesteps); %transeverse sediment flux from cell b1 to cell c1 (immediately upstream of bifurcation)
etab_save=zeros(1,timesteps); %saves bed elevation is cell b1 
etac_save=zeros(1,timesteps);
etab_branchav_save=zeros(1,timesteps);%saves average bed elevation from cells 2 to N
etac_branchav_save=zeros(1,timesteps);
tausb_save=zeros(1,timesteps);%saves shields stress at branch entrance
tausc_save=zeros(1,timesteps);
Hb_save=zeros(1,timesteps);%water depth at branch entrance
Hc_save=zeros(1,timesteps);
Sb_save=zeros(1,timesteps);%slope at branch entrance
Sc_save=zeros(1,timesteps);
frontb=zeros(1,timesteps);%front position (center position of final subaerial cell)
frontc=zeros(1,timesteps);
if variable_chezy==1
Czb_save=zeros(1,timesteps);%value of chezy coefficient in each branch
Czc_save=zeros(1,timesteps);
end

%profiles - collected periodically
Nprofiles=10; %during run, how many total profiles to collect?
etab_profile=zeros(Nprofiles,N_b);
etac_profile=zeros(Nprofiles,N_c);
Qsb_profile=zeros(Nprofiles,N_b+1);
Qsc_profile=zeros(Nprofiles,N_c+1);
t_save=zeros(1,Nprofiles); 

if makemovie==1;
counter_movie=0;%used for movie making
vid_interval=1000;%every x timesteps create frame for movie
end

counter_profile=0;%used for saving profiles

%mass conservation check
mass_initial=(eta_b(1)*B_b*dx_up+eta_c(1)*B_c*dx_up+sum(eta_b(2:N_b).*B_b)*dx+sum(eta_c(2:N_c).*B_c)*dx)*(1-lamp);
aggradationrate=(2*(.5*Qs_a-F*Qsi0)/Qs_a)*Qs_a/((dx_up*B_b+dx_up*B_c+L_b*B_b+L_c*B_c)*(1-lamp))*timeyear;

t_star=t*timeyear*Qs_a/(B_a*alpha*B_a*H_a*(1-lamp)); %dimensionless time 

%% time evolution

for j=1:timesteps
%%%apply perturbation to kick away from symmetric solution
    if j==timestep_perturb     
            eta_b=eta_b+perturb_amp; %bed elevation in branch b is perturbed upwards
            eta_c=eta_c-perturb_amp; %bed elevation in branch c is perturbed downwards
    end

%%%ensure that BRT boundary condition is still applied following perturbation
    if down_bc==1
        eta_b(N_b)=eta_downstream_average+(eta_b(1)-eta_c(1))/2;
        eta_c(N_c)=eta_downstream_average-eta_b(N_b);
    end
    
%%%find front position 
%finds center position of final subaerial cell center
%reference frame is x=0 at the upstream boundary

    f2=find(eta_b<=sealevel_b(j));%find all inundated cells  
    if isempty(f2)==0%at least one cell is inundated
        if f2(1)==1%all cells inundated
            frontb(j)=0;
            fb=0;
        else
        fb=f2(1)-1;%index of final subaerial cell, used later to determine where to compute depths
        frontb(j)=xb(fb); %front position (m)
        end   
    else %no cells inundated
        fb=N_b-1;%subtract 1 from index, because when used later, we set the downstream boundary separately
        frontb(j)=xb(fb+1);
    end
    
    f3=find(eta_c<=sealevel_c(j));%find all inundated cells  
    if isempty(f3)==0%at least one cell is inundated
        if f3(1)==1%all cells inundated
            frontc(j)=0;
            fc=0;
        else
        fc=f3(1)-1;%index of final subaerial cell 
        frontc(j)=xc(fc);
        end   
    else %no cells inundated
        fc=N_c-1;
        frontc(j)=xc(fc+1);
    end
        
%%%iteration to find xi_b=xi_c (constant water surface at bifurcation point)
    
    S_b=(eta_b(1)-eta_b(2))/(.5*dx_up+.5*dx);%slopes defined at cell boundaries
    S_c=(eta_c(1)-eta_c(2))/(.5*dx_up+.5*dx);
    if S_b<=0 %due to normal flow assumption, discharge is set to 0 if slope is negative     
            Q_b=0;
            Q_c=Q_a;
            Hb_temp=0;
            Hc_temp=(Q_c^(2/3))/((B_c^(2/3))*(Cz_c(2)^(2/3))*(S_c^(1/3))*(g^(1/3)));           
    elseif S_c<=0    
            Q_b=Q_a;
            Q_c=0;
            Hc_temp=0;
            Hb_temp=(Q_b^(2/3))/((B_b^(2/3))*(Cz_b(2)^(2/3))*(S_b^(1/3))*(g^(1/3)));     
    else %both slopes are positive
        Hb_temp=(Q_b^(2/3))/((B_b^(2/3))*(Cz_b(2)^(2/3))*(S_b^(1/3))*(g^(1/3)));%normal flow
        Hc_temp=(Q_c^(2/3))/((B_c^(2/3))*(Cz_c(2)^(2/3))*(S_c^(1/3))*(g^(1/3)));
        
        if variable_chezy==1 
            Cz_b(2)=max(minchezy,6+2.5*log(Hb_temp/(2.5*ds)));%update branch chezy based on flow depth 
            Cz_c(2)=max(minchezy,6+2.5*log(Hc_temp/(2.5*ds)));
        end

    xi_b=Hb_temp+eta_b(1);%water surface elevation at cell boundary
    xi_c=Hc_temp+eta_c(1);
                    
    while abs(xi_b-xi_c)>watlev_accuracy %iteration to find fluxes Q2 and Q3 such that xi_b=xi_c (to specified accuracy)

                if Q_b==0%allows for recovery from complete avulsion (otherwise, iteration fails)
                    Q_b=Q_a-recovery*Q_a;
                elseif Q_c==0
                    Q_c=Q_a-recovery*Q_a;
                    Q_b=Q_a-Q_c;
                end
                
                %improve estimate of Q2 and Q3
                %with newton raphson iteration method
                Q_b=Q_b-(3/2)*(xi_b-xi_c)/((Q_b^(-1/3))/((B_b^(2/3))*(Cz_b(2)^(2/3))*(S_b^(1/3))*(g^(1/3)))+((Q_a-Q_b)^(-1/3))/((B_c^(2/3))*(Cz_c(2)^(2/3))*(S_c^(1/3))*(g^(1/3))));%newton raphson iteration
                Q_c=Q_a-Q_b;
                %if discharge is below a threshold, complete avulsion
                %occurs and exit while loop
                if Q_b<0
                    Q_b=0;
                    Q_c=Q_a;
                    Hb_temp=0;
                    Hc_temp=(Q_c^(2/3))/((B_c^(2/3))*(Cz_c(2)^(2/3))*(S_c^(1/3))*(g^(1/3)));

                    break
                end
                if Q_c<0
                    Q_c=0;
                    Q_b=Q_a;
                    Hc_temp=0;
                    Hb_temp=(Q_b^(2/3))/((B_b^(2/3))*(Cz_b(2)^(2/3))*(S_b^(1/3))*(g^(1/3)));

                    break
                end
                %after finding Q2 and Q3
                Hb_temp=(Q_b^(2/3))/((B_b^(2/3))*(Cz_b(2)^(2/3))*(S_b^(1/3))*(g^(1/3)));%normal flow
                Hc_temp=(Q_c^(2/3))/((B_c^(2/3))*(Cz_c(2)^(2/3))*(S_c^(1/3))*(g^(1/3)));
   
                if variable_chezy==1 %update chezy coefficient 
                Cz_b(2)=max(minchezy,6+2.5*log(Hb_temp/(2.5*ds)));
                Cz_c(2)=max(minchezy,6+2.5*log(Hc_temp/(2.5*ds)));
                end
         
                xi_b=Hb_temp+eta_b(1);%recalculate water surface elevations
                xi_c=Hc_temp+eta_c(1);
                                  
    end 
    end 
    
%%%cross-stream flux computation
    Qy=(B_b/(B_c+B_b))*Q_a-Q_c;%water discharge from cell c1 to b1 (from continuity)
    H_b(2)=(Q_b^(2/3))./((B_b^(2/3))*(Cz_b(2)^(2/3))*((S_b).^(1/3))*(g^(1/3)));%water depth at bifurcation point
    H_c(2)=(Q_c^(2/3))./((B_c^(2/3))*(Cz_c(2)^(2/3))*((S_c).^(1/3))*(g^(1/3)));
    
    Sy=(eta_c(1)-eta_b(1))/((B_b+B_c)/2); %transverse bed slope
    
    if Habc_set==1 %average water level Habc (used in transverse sed xport formula)
        Habc=.5*(H_b(2)*B_b/(B_b+B_c)+H_c(2)*B_c/(B_b+B_c)+H_a);
    else
        Habc=H_a;
    end

    Qsy=Qs_a*(Qy*H_a/(Q_a*Habc)+r*alpha*Sy/sqrt(taus_a)); %transverse sed flux from cell c1 to b1

    %compute the rest of the depth vector  

    H_b(3:fb+1)=(Q_b^(2/3))./((B_b^(2/3)).*(Cz_b(3:fb+1).^(2/3)).*(((eta_b(2:fb)-max(sealevel_b(j),eta_b(3:fb+1)))/dx).^(1/3))*(g^(1/3)));%downstream reach
    H_b(fb+2:N_b)=sealevel_b(j)-eta_b(fb+2:N_b);%depth beyond shoreline
    H_c(3:fc+1)=(Q_c^(2/3))./((B_c^(2/3)).*(Cz_c(3:fc+1).^(2/3)).*(((eta_c(2:fc)-max(sealevel_c(j),eta_c(3:fc+1)))/dx).^(1/3))*(g^(1/3)));%downstream reach
    H_c(fc+2:N_c)=sealevel_c(j)-eta_c(fc+2:N_c);
    if variable_chezy==1
    %the error due to not iterating for the values of Cz and H is very
    %small because Cz changes slowly.
        Cz_b=max(minchezy,6+2.5.*log(H_b/(2.5*ds)));
        Cz_c=max(minchezy,6+2.5.*log(H_c/(2.5*ds)));
        if transport_form==2; %change n to reflect change in Cz in EH formula
            nb=.05*Cz_b.^2;
            nc=.05*Cz_c.^2;
        end
    end
    
    %compute shields stress vector
    taus_b(2)=p*((Q_b).^2)./((ps-p)*g*ds*(Cz_b(2)^2)*(B_b^2)*(H_b(2).^2));
    taus_c(2)=p*((Q_c).^2)./((ps-p)*g*ds*(Cz_c(2)^2)*(B_c^2)*(H_c(2).^2));
    taus_b(3:fb+1)=p*((Q_b).^2)./((ps-p)*g*ds*(Cz_b(3:fb+1).^2)*(B_b^2).*(H_b(3:fb+1).^2));
    taus_b(fb+2:N_b)=0;%beyond shoreline taus=0,meaning no flux out of first underwater cell
    taus_c(3:fc+1)=p*((Q_c).^2)./((ps-p)*g*ds*(Cz_c(3:fc+1).^2)*(B_c^2).*(H_c(3:fc+1).^2));
    taus_c(fc+2:N_c)=0;
    
    %sed flux vector
    for i=2:N_b %for loop to compute sediment flux in downstream branches using MPM
       if taus_b(i)>tausc
           Qs_b(i) = (nb(i)*(taus_b(i)-tausc).^m)*((Rr * g * ds) ^ 0.5) * ds*B_b;
       else
           Qs_b(i)=0;
       end
    end
    for i=2:N_c
       if taus_c(i)>tausc
           Qs_c(i) = (nc(i)*(taus_c(i)-tausc).^m)*((Rr * g * ds) ^ 0.5) * ds*B_c;
       else
           Qs_c(i)=0;
       end
    end

    if down_bc==2 %Compute outgoing sediment flux at downstream boundary  
         Qs_b(N_b+1)=F*Qs_b(N_b)/((F-1)*(-dx/(L_b))+F);
         Qs_c(N_c+1)=F*Qs_c(N_c)/((F-1)*(-dx/(L_c))+F);
    end
    
%%%update bed elevations with exner (mass conservation)
    eta_b(1)=eta_b(1)+(1/(1-lamp))*(Qs_b(1)-Qs_b(2))*dt/(dx_up*B_b)+dt*(1/(1-lamp))*(Qsy)/(dx_up*B_b);%alpha*B_a reach
    eta_c(1)=eta_c(1)+(1/(1-lamp))*(Qs_c(1)-Qs_c(2))*dt/(dx_up*B_c)-dt*(1/(1-lamp))*(Qsy)/(dx_up*B_c);
    eta_b(2:N_b)=eta_b(2:N_b)+(1/(1-lamp))*dt*(Qs_b(2:N_b)-Qs_b(3:N_b+1))./(dx*B_b);%downstream reach
    eta_c(2:N_c)=eta_c(2:N_c)+(1/(1-lamp))*dt*(Qs_c(2:N_c)-Qs_c(3:N_c+1))./(dx*B_c); 
    
    if down_bc==1 %ensure that water levels are equal at downstream boundary
        eta_b(N_b)=eta_downstream_average+(eta_b(1)-eta_c(1))/2;
        eta_c(N_c)=2*eta_downstream_average-eta_b(N_b);
    end
    
    %%%apply differential subsidence
    eta_b(2:N_b)=eta_b(2:N_b)+tilt*aggradationrate*dt_years;
    eta_c(2:N_c)=eta_c(2:N_c)-tilt*aggradationrate*dt_years;

    
    
%%%save useful data
    Qb_save(j)=Q_b;
    Qc_save(j)=Q_c;
    Qsb_save(j)=Qs_b(2);
    Qsc_save(j)=Qs_c(2);
    Qsy_save(j)=Qsy;
    etab_save(j)=eta_b(1);
    etac_save(j)=eta_c(1);
    tausb_save(j)=taus_b(2);
    tausc_save(j)=taus_c(2);
    etab_branchav_save(j)=mean(eta_b(2:N_b));
    etac_branchav_save(j)=mean(eta_c(2:N_c));
    Hb_save(j)=H_b(2);
    Hc_save(j)=H_c(2);
    Sb_save(j)=S_b;
    Sc_save(j)=S_c;
    
    if variable_chezy==1
    Czb_save(j)=Cz_b(2);
    Czc_save(j)=Cz_c(2);
    end

%%%save profiles
     if rem(j,timesteps/Nprofiles)==0
        counter_profile=counter_profile+1;
        etab_profile(counter_profile,1:N_b)=(eta_b(1:N_b));
        etac_profile(counter_profile,1:N_c)=(eta_c(1:N_c));
        Qsb_profile(counter_profile,1:N_b+1)=Qs_b(1:N_b+1);
        Qsc_profile(counter_profile,1:N_c+1)=Qs_c(1:N_c+1);
        t_save(counter_profile)=t(j);
        j %provides update on how far the iteration has progressed
        (eta_b(1)-eta_c(1))/H_a %update on current state of the system
     end
%%%make movie
    if makemovie==1
    if rem(j,vid_interval)==0
        counter_movie=counter_movie+1;
        fig=figure;      
        subplot(2,1,1);
        hold on
        P1=plot((xb-alpha*B_a)/B_a,eta_b/B_a,'Color',[0    0.4470    0.7410]);
        P1.LineWidth=2;
        hold on
        P2=plot((xc-alpha*B_a)/B_a,eta_c/B_a,'Color',[0.8500    0.3250    0.0980],'LineStyle','--');
        P2.LineWidth=2;
        P3=line([-alpha (max(xb)-alpha*B_a)/B_a],[sealevel_b(j)/B_a sealevel_b(j)/B_a],'Color','blue','LineStyle',':');
        P3.LineWidth=1;
        axis([-alpha  (max(xb))/B_a -1/B_a 45/B_a])
        xlabel('x/B_a')
        ylabel('y/B_a')
        
        subplot(2,1,2);
        P1=plot(t_star(1:j),Qb_save(1:j)/Q_a,'Color',[0    0.4470    0.7410]);
        P1.LineWidth=2;
        hold on
        P2=plot(t_star(1:j),Qc_save(1:j)/Q_a,'Color',[0.8500    0.3250    0.0980],'LineStyle','--');
        P2.LineWidth=2;
        xlabel('t_*')
        ylabel('Q_i/Q_a')
        axis([0 t_star(end) 0 1])
     
        frame(counter_movie)=getframe(fig);
        close all
        hold off
    end
    end
    

end
%% finalize
%%%check mass balance
mass_fed=Qs_a*timesteps*dt;
mass_actual=(eta_b(1)*dx_up*B_b+eta_c(1)*dx_up*B_c+sum(eta_b(2:N_b))*B_b*dx+sum(eta_c(2:N_c)).*B_c*dx)*(1-lamp);%actual final mass
fractionretained=(mass_actual-mass_initial)/mass_fed
fractionretained_expected=2*(.5*Qs_a-F*Qsi0)/Qs_a %long-term average expected fraction retained for down_bc=2

if makemovie==1
writerObj = VideoWriter('downstreamsink_prograde.avi'); % Name it.
writerObj.FrameRate = 30; % How many frames per second.
writerObj.Quality =100;
open(writerObj); 
writeVideo(writerObj, frame);
close(writerObj)
end

%%%compute dimensionless variables
deta_star=(etab_save-etac_save)/B_a;
deta_starL=(etab_branchav_save-etac_branchav_save)/B_a;


