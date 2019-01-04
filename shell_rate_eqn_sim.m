%% This function evaluates the plasma rise time for given Ry densities
% I set n_min as lower decay boundry. States reaching levels <= n_min
% are considered to decay immediately

% den0   = arrays of Ry atom densities for all shells
% vol    = volume of each shell
% n0     = initial PQN
% t_final= final time in ns
% dt     = time steps

% t is time array from 0ns to t_final ns
% nden is matrix with densities of PQN n states.
 %      first 100 values are for first shell, next 100 values for next
 %      shell and so on
% eden - aray with electron densities
% deac - array with densities of deactivated states (n<10)
% Te   - array with electron temperatures

% i.e.: [t,nden,eden,aden,Te,~]=shell_rate_eqn_sim(0.3,50,50);

%save workspace save('den0=1 n0=50, t0=50.mat','t','nden','eden','aden','Te')

function [t,nden,eden,deac,Te,y0]=shell_rate_eqn_sim(den0, vol ,n0, dt, t_final) 

%%%%%%%%%%%%%%%%%%%%
% Penning fraction %
%%%%%%%%%%%%%%%%%%%%

    % This function calculates the inital seed electrons from Penning fraction

    % n are initial PQN-levels before Penning ionization
    % den is Ry-atom density in 10^12 pcc

    % eden is the # of electrons produced
    % rden are the remaining PQN-levels after Penning ionization

function [PenningFraction eden rden]=penningfraction(n,den)     % den in 10^12 pcc

    a0=5.2917721092e-5;     % bohr radius in um
    Rn0=n.^2*a0;            % radius of Rydb. atom by bohr model using semi-classical method
    Rmax=1.8*(Rn0*2);       % Robicheaux paper, within this distance, 90% penning ionize (~10ns)

    %PenningFraction=zeros(length(n),1);
    PenningFraction_analytic=zeros(length(n),1);

    % Calculates number of Penning partners based on Erlang distribution:
    for i=1:length(n)
    %PenningFraction(i)=quadv(@(r)4*pi*den*r.^2.*exp(-4*pi*den*r.^3/3),0,Rmax(i),1e-10); % proportion between 0 and Rmax
    PenningFraction_analytic(i)=1-exp(-4*pi*den*Rmax.^3/3); % integral is solved analytically
    end

    PenningFraction = PenningFraction_analytic;         % 90% ionization within certain time (assume rest non-interactive)
    eden=PenningFraction/2*den*.9;      % the dens. of electrons produce is half the proportion (1e- per partner)
    rden=(1-PenningFraction*0.9)*den;    % this is remaining density of rydbergs

end

%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions %
%%%%%%%%%%%%%%%%%%%%%%

N=length(den0); %number of shells

% Set constants and lower cut-off for PQN distribution:
kB = 1.3806504e-23;                 % #m2 kg s-2 K-1;J K-1; 8.62e-5 #eV K-1,
Ry = 2.179872e-18;                  % #J: Rydberg energy in J 

firstn=1;
n_min=10;
numlev=100;                         % this is the number of n levels initially considered
deac=0;                             % start with all Rydberg states
nl=(firstn:firstn+numlev-1)';       % array of accessible n levels

ns=length(nl);                      %number of n-states

NL=zeros(ns*N,1);
VOL=zeros(ns*N,1);
DEN0=zeros(ns*N,1);
NDEN=NL*0;
EDEN=zeros(N,1);
RDEN=zeros(N,1);
DEAC=zeros(N,1);
T_PENNING=zeros(1,N);

tspan=0:dt:t_final;%linspace(0,t_final,500);

% Sets initial conditions for electron and n-level distributions:
% no initial electrons -> calc. Penning seed electrons (look up Robicheaux)
for ii=1:N %loop though all shells
    [PF,eden,rden]=penningfraction(n0,den0(ii));
    % Redistributes the Penning partners over lower n's:
    f=@(x)5.95*x.^5;                % This is the penning fraction distribution
    np=firstn:fix(n0/sqrt(2));      % Array of n states allowed after Penn ion
    ind=1:length(np);               % This is the distribution of penning fraction
    nden=nl*0;
    nden(ind)=eden*f(np/n0)/sum(f(np/n0));  % dividing by the sum normalizes the function

    nden(nl==n0)=rden;              % set n0 to rden
    
    % Set initial temperature:    (Robicheaux 2005 JPhysB)
    T_PENNING(ii)=(-Ry*den0(ii)/n0^2 + Ry*rden/n0^2 + Ry*sum(nden(ind)./nl(ind).^2) )*1/(3/2*kB*eden); % by energy conservation

    deac=sum(nden(1:n_min));        % allow n<=n_min to decay to aden
    nden(1:n_min)=zeros(n_min,1); 
    
    NL(1+(ii-1)*ns:ii*ns)=nl;       % save values for this shell in arrays
    VOL(1+(ii-1)*ns:ii*ns)=vol(ii);
    DEN0(1+(ii-1)*ns:ii*ns)=den0(ii);
    NDEN(1+(ii-1)*ns:ii*ns)=nden;     
    EDEN(ii)=eden;
    DEAC(ii)=deac;
end

%the penning temperature is equal for all densities! (density scaling
%factor cancels)
T_penning = sum(T_PENNING.*vol.*den0)/sum(vol.*den0); %calculate equlibrated penning temperature by weigted average

%%%%%%%%%%%%%%%
% rate coeffs %
%%%%%%%%%%%%%%%

function [ni,nf,II,minn,maxn,diffsn]=buildns(nl)
    
    % Is needed to reevaluate the rate coeff.
    a=length(nl);
    II=ones(a,a); % will be matrix of ij=1 and ii=0
    ni=zeros(a,a);
    nf=zeros(a,a);
    minn=zeros(a,a);
    maxn=zeros(a,a);
    for i=1:a
        for j=1:a
            ni(i,j)=nl(i); % an array of initial states
            nf(i,j)=nl(j); % an array of final state
            minn(i,j)=min(ni(i,j),nf(i,j)); %find min of init and final state potential problem
            maxn(i,j)=max(ni(i,j),nf(i,j)); %find max of init and final state
            if i==j
            II(i,j)=0; % set to 0 ones with same init and final state
            end
        end
    end
    diffsn=abs(1./ni.^2-1./nf.^2); % difference in energy between the 2 states (no units)

end % help array

function ktbr=kTBR(n,T)     % three-body-recombination
    
    % TBR rates output units in um^6 ns-1
    kB = 1.3806504e-23;% #m2 kg s-2 K-1;J K-1; 8.62e-5 #eV K-1,
    emass = 9.1093822e-31;% #kg
    h = 6.6260690e-34;% #m2 kg / s 4.13e-15 #eV s
    Rydhc = 2.179872e-18;% #J: Ryd [cm-1] --> [J]

    epsi=Rydhc./(power(n,2.0)*kB*T);
    lmbd=1e6*h./sqrt(2.0*pi*emass*kB*T);

    ktbr=kION(n,T).*power(n,2.0).*power(lmbd,3.0).*exp(epsi); % unit um^6 ns-1   *0.1 to test effect of TBR suppression
    ktbr(isnan(ktbr))=0;       % take care of computation error for values matlab sees too small and turns NaN
    ktbr(~isfinite(ktbr))=0;    % take care of computation error for values matlab sees too large and turns inf
end

function kion=kION(n,T) % ionizing collisions
    
    % output unit is um^3
    kB = 1.3806504e-23;% #m2 kg s-2 K-1;J K-1; 8.62e-5 #eV K-1,
    Rydhc = 2.179872e-18;% #J: Ryd [cm-1] --> [J]
    epsi=Rydhc./(power(n,2)*kB*T); % find reduced initial energy
    kion=11*sqrt(Rydhc./(kB*T)).*kNOT(T).*exp(-epsi)./(power(epsi,2.33)+4.38*power(epsi,1.72)+1.32*epsi);

end

function out=knnp(ni,nf,II,minn,maxn,diffsn,T) % rate for transfer from n to n'
    
    % Rates calculated using PVS PRL 2008 
    % output units is um^3 ns^-1
    % use in conjunction with [ni,nf,II,minn,maxn,diffsn]=buildns(nl);

    kB = 1.3806504e-23;% #m2 kg s-2 K-1;J K-1; 8.62e-5 #eV K-1, % PVS PRL 2008
    Rydhc = 2.179872e-18;% #J: Ryd [cm-1] --> [J]

    eps_i=Rydhc./(power(ni,2.0)*kB*T);
    eps_f=Rydhc./(power(nf,2.0)*kB*T);
    max_eps_i=Rydhc./(power(minn,2)*kB*T);
    min_eps_i=Rydhc./(power(maxn,2)*kB*T);

    diffs=Rydhc.*diffsn./(kB*T); % scale dffsn properly

    out=II.*(kNOT(T).*power(eps_i,5/2).*power(eps_f,3/2)./power(max_eps_i,5/2))...
        .*exp(-(eps_i-min_eps_i)).*((22./(power(max_eps_i+0.9,7/3)))...
        +(4.5./(power(max_eps_i,2.5).*power(diffs+1-II,4/3)))); 

end    

function knot=kNOT(T)
    
    el = 1.6021765e-19;% #C
    kB = 1.3806504e-23;% #m2 kg s-2 K-1;J K-1; 8.62e-5 #eV K-1,
    epsilon = 8.854187817e-12;% #C2 J-1 m-1
    Rydhc = 2.179872e-18;% #J: Ryd [cm-1] --> [J]
    emass = 9.1093822e-31;% #kg
    
    knot=1e9*power(el,4.0)./(kB*T*sqrt(emass*Rydhc)*power(4*pi*epsilon,2)); % PRL 2008 PVS 
    % orginal units of knot = m^3/s
    % 1e9 converts it into um^3 / ns
    % dividing by (4\pi\epsilon)^2 converts a.u. --> s.i.
end

%%%%%%%%%%%%%%%%
% Calculations %
%%%%%%%%%%%%%%%%

y0=[NDEN;EDEN;DEAC;T_penning]; % set initial value for ODE

ncrit=@(T)round(sqrt(Ry/(kB*T)));  % to calculate n-max (got by fitting)

[ni,nf,II,minn,maxn,diffsn]=buildns(nl);   % is needed to calculate the rate coeffs

progress = waitbar(0,'Progress...');

function dy=eqrateode(t,y)
    
    % fprintf('%f\n',t);
    
    % Select valiables:
    nden=y(1:ns*N);       % pick out density over distribution of n
    eden=y(ns*N+1:(ns+1)*N);       % pick out electron density
    deac=y((ns+1)*N+1:(ns+2)*N);       % pick out density of radiatively decayed atoms
    T=y((ns+2)*N+1);          % pick out temperature
    % all of the above variables are arrays with data for all shells
    
    nc=ncrit(T);            % calculates n max with this temperature
        
    % Adjusts max allowed n: 
    if nc>=nl(1) && nc<nl(end) 
        index=find(nl==ncrit(T));
    elseif nc<nl(1) 
        index=1;
    else
        index=numlev;
    end
   
    %preallocate arrays for all shells
    D_NDEN=zeros(ns*N,1);
    D_NDEN_OLD=zeros(ns*N,1);
    D_EDEN=zeros(N,1);
    D_DEAC=zeros(N,1);

    
    % Evaluate the updated rate terms:
    for i=1:N
        
        h=1+(i-1)*ns;
        k=i*ns;
        
        d_tbr=zeros(numlev,1);
        d_tbr(1:index)=kTBR(nl(1:index),T)*eden(i)^3; % units [d_tbr] = um^-3 ns^-1
    
        d_ion=kION(nl,T).*nden(h:k)*eden(i); % units [d_ion] = um^-3 ns^-1
        
        % rate for transfer from n to n', unit [kn_np] = um^3 ns^-1
        k_n_np=knnp(ni,nf,II,minn,maxn,diffsn,T); 
        d_n_np=sum(k_n_np,2).*nden(h:k)*eden(i);            %[um^-3 ns^-1]
        
        % rate for transfer from n' to n, [knp_n] = ns^-1
        k_np_n=zeros(numlev,numlev);
        for j=1:index % only to levels <= nc
            k_np_n(j,1:numlev)=k_n_np(1:numlev,j).*(nden(h:h+numlev-1)); 
        end
        d_np_n=sum(k_np_n,2)*eden(i);                  %[um^-3 ns^-1]
    
        % transfer from n's above ncrit(T) to eden
        k_n_npion=zeros(numlev,1);
        if index<=numlev 
            k_n_npion(1:index)=sum(k_n_np(1:index,index+1:numlev),2).*nden(h:h+index-1); % [kn_npion] = ns^-1
        end
        d_n_npion=sum(k_n_npion)*eden(i);              %[um^-3 ns^-1]
        
        
        D_EDEN(i)=sum(d_ion-d_tbr)+d_n_npion; 

        d_nden=d_tbr-d_ion-d_n_np+d_np_n;
        D_NDEN_OLD(h:k)=d_nden; %array whithout radiative decay for temperature calculation
            % Implements radiative decay/PD for levels n <= n_min:
        D_DEAC(i)=sum(d_nden(1:n_min));     % aden is the number of Ry's decayed radiatively to the groundstate
        d_nden(1:n_min)=zeros(n_min,1); 
        
        D_NDEN(h:k)=d_nden;
    end

    %now calculate change in equilibrated temperature
    dT=(Ry*sum(D_NDEN_OLD./NL.^2.*VOL)-1.5*kB*T*sum(D_EDEN.*(vol')))/(1.5*kB*sum(eden.*(vol')));
  
    dy=[D_NDEN;D_EDEN;D_DEAC;dT];
    
    waitbar(t/t_final); 

end

%%%%%%%
% ODE %
%%%%%%%

options=odeset('reltol',1e-8);
[t,y]=ode23(@(t,y)eqrateode(t,y),tspan,y0);

nden=y(:,1:ns*N);       % pick out density over distribution of n
eden=y(:,ns*N+1:(ns+1)*N);       % pick out electron density
deac=y(:,(ns+1)*N+1:(ns+2)*N);       % pick out density of radiatively decayed atoms
Te=y(:,(ns+2)*N+1);   

close(progress);

end

