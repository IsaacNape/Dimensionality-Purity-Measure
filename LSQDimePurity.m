function [GuessedParameters,sigm, Rsq]=LSQDimePurity(SpectrumTypeIndex, VM, n)
% This code returns 2, 2D arrays containing the GuessedParameters=(p, K) and sigm=(sigmp, sigmK) and the R^2 goodness of fit 
%VM are the measured visibilities containing  with dimensions row==n measurements  , col == samples
%nIndices are the indexes of the analysers used
 
assumedSpectrum{1} = @PGaussVisi; %% Normal distribution model for pure part 
assumedSpectrum{2} = @PSPDCVisi; %%SPDC distribution model for pure part 
assumedSpectrum{3} = @PMaxENTVisi %%Max ent. ;
%Spectyptag = {'Gauss', 'SPDC', 'MaxEnt'};

if length(n)== length(VM)
    % calculate mean of measured samples 
    if sum(size(VM)==1)==0
        VMmean=zeros(length(n), 1);
        VMmean(:)=mean(VM, 2); %take mean of the visibilities row==n measurements  , col== sampless
        VMstad=std(VM');% extract the uncertainties from samples
    else
        disp('Not enough, sampled data points... I will continue with the calculation anyway :)')
    end
    
    %determining best fit
    [GuessedParameters,~,resid,~,~,~,J]=lsqcurvefit( assumedSpectrum{SpectrumTypeIndex}, [0.1, 5], n, VMmean', [0,2], [1, 100]);
    
    % error intervals of p and K
    sigm95 = nlparci(GuessedParameters,resid,'jacobian',J);% Calculation of 95% confidence intervals
    sigm=abs(GuessedParameters'-sigm95(:,1))./2./1.96; %% standard error
    sigm90 = sigm*1.645;%% standard error fo 90 % interval;
    
   % Goodness of fit Calculation
    SStot = sum((VMmean-mean(assumedSpectrum{SpectrumTypeIndex}(GuessedParameters, n))).^2);  % Total Sum-Of-Squares
    SSres = sum((VMmean(:)'-assumedSpectrum{SpectrumTypeIndex}(GuessedParameters, n)).^2);   % Residual Sum-Of-Squares
    Rsq = 1-SSres/SStot; %% goodness of fit -> R^2
else
    disp('Number of analysers (n) does not match number of visibility measurements (V_n)');
end
