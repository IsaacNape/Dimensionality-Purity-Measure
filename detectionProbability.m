function Prob=detectionProbability(p, Kinp, theta, N, Spectype)
% p the purity of the state
% Kinp the dimensions the state
% Spectype 1 for gaussian and 2 for SPDC and 3 for maximally ent. (uniform)


% Analyser weightings for OAM modes (L) depending on N (number of analysers)
CN=@(L, oam, u, theta, N) 1i*exp( - 1i * (L.*((N-1)*pi + N.*theta))./N ) .*(exp(2*(1i*pi*u))).*sin(u*pi)./ (2*N.* (-oam+L-u)).* Al(L,N);
if Spectype == 1
    Spectrum = @(l, gamma) exp(-l.^2 ./ (gamma./2).^2);
    gamma=Kinp./1.2533; % mapping between K and gamma for gaussian dist.
    Dmax=ceil(Kinp); %% Set the max dimensions for to be the rounded up effective dimensions (Kinput)
    L=-Dmax:1:Dmax; % list of OAM indexes
elseif  Spectype == 2
    Spectrum = @(l, gamma) (2*gamma.^2./(1+2*gamma.^2)).^abs(2*l); 
    gamma = sqrt((Kinp-1)./4); % mapping between K and gamma SPDC dist.
    Dmax=ceil(Kinp); %% Set the max dimensions for to be the rounded up effective dimensions (Kinput)
    L=-Dmax:1:Dmax; % list of OAM indexes
elseif  Spectype == 3
    Spectrum = @(l, gamma) double(abs(l)>=0);
    gamma = (Kinp-1)./2; % mapping between K and gamma of uniform dist dist.
    Dmax = ceil((Kinp-1)./2); %% Set the max dimensions for to be the rounded up effective dimensions (Kinput)
    L=-Dmax:1:Dmax; % list of OAM indexes
    if mod(Kinp, 2)==0
        Dmax=Kinp./2;
        L=[-Dmax : 1 : Dmax-1];
    else
        Dmax=(Kinp-1)./2;
        L=[-Dmax:1:Dmax];
    end
else
    disp('Warning: Spectype must be 1 or 2 to select gaussian or SPDC spec, respectively!')
    return;
end 

%Angles to be considered in the projections. To be later divided by n (analyser no. of superpositions)
Nmaxtheta = length(theta);

% Set integer part of OAM
Mint=fix(N/2);

% Setting the OAM width depending on efefctive dimensions Kinp

Spec=Spectrum(L, gamma)./sum(Spectrum(L, gamma));

%Initialise probabilities (the pure and mixed parts)
probPure=zeros(Nmaxtheta,length(N));  %purepart
overlapMx=zeros(Nmaxtheta,length(N)); %mixed part

%Normalisation factors
temp = zeros(length(N));
for Nindex = 1:length(N)
        temp(Nindex) =  sum(abs(CN( L, Mint(Nindex), rem(N(Nindex)/2,1), 0, N(Nindex))).^2, 2);
end
NpureNorm =1./sqrt(temp);

%Calculate probabilities
for Nindex = 1:length(N)
    for theta_index = 1 : length(theta)
        M1 = CN( L, Mint(Nindex), rem(N(Nindex )/2,1), 0, N(Nindex ) ); %overlap photon 1
        M2 =  CN(-L, -Mint(Nindex), -rem(N(Nindex )/2,1), theta(theta_index), N(Nindex ));%overlap photon 2
        overlapA = M1.*M2.*sqrt(Spec); %overlap with spectrum of state (Gaussian spectrum)
        overlapMx(theta_index, Nindex) = abs(sum( (abs(M1).^2 .*NpureNorm(Nindex ).^2))).^2; %mixed part overlap probability
        probPure(theta_index, Nindex) =  abs(sum( (overlapA.*NpureNorm(Nindex ).*NpureNorm(Nindex )))).^2; %pure part overlap probability
    end
end
Prob = p.*probPure + overlapMx.*(1-p)./(Kinp).^2;
Prob(isnan(Prob))=0;