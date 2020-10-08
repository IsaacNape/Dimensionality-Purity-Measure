function Visi=PSPDCVisi(a, N)
%% This code calculates the visibility for an inputmixed state, rho(p=a(1), and K=a(2)).

p=a(1); %% the purity of the state
Kinp=a(2);%% the dimensions the state

%vAngles to be considered in the projections. To be later divided by n 
theta=[0, pi];
Nmaxtheta = length(theta);

%Spectype=1 for gaussian
Spectype=2;

%Initialise probabilities (the pure + mixed parts)
Prob=zeros(Nmaxtheta ,length(N));  %purepart

%Calculate probabilities
for Nindex = 1:length(N)
        Prob(:,Nindex) = detectionProbability(p, Kinp, theta./N(Nindex), N(Nindex), Spectype);
        Visi(Nindex)=abs(Prob(1,Nindex)-Prob(2,Nindex))./(Prob(1,Nindex)+Prob(2,Nindex));%
end
   %(Nindex)=abs(Prob(1,Nindex)-Prob(2,Nindex))./(Prob(1,Nindex)+Prob(2,Nindex));%(p*(DeltaP(Nindex))')./(p*(SumP(Nindex))' + 2*(1-p)./(Kinp.^2).*ProbMx(1, Nindex)); %visibility with consideration of mixed part

Visi(isnan(Visi))=0;