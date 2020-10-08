%%  This file reads the coincidence data and extracts the visibilities and stores both sets in a a .matfile
%% The data is stored as an array -> row == V_n (n=1,3,...2nmax-1)
%% The files are names as nSPP_n_#n_frac_0.5_int_#m. Here #n is n=1,3,5.... while #m is the interger part of #n/2

FolderDir='20200824';%parent directory
maxn=11; %% maximum number of visibilities measured (V_{n max})

for n=1:2:maxn
%read the txt file
filename = [ FolderDir '\nSPP_n_',num2str(n),'_frac_0.50_int_',num2str(fix(n/2)), '\']; %file directory
Coin=importdata([filename 'Coincidences.txt']); %read data

CoinMax = max(Coin); % get peak
CoinMin = min(Coin); % get trough

Visibilities((n+1)./2, :) = abs(CoinMax-CoinMin)./(CoinMax+CoinMin); % STORE VISIBILITIES
end 

save(['Data', FolderDir,'.mat'], 'Coincidence', 'Visibilities') %save visibilities in .mat file

