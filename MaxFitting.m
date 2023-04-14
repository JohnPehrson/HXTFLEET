function [fitvariables] = MaxFitting(row_intensity,row_index,limits)
%This function fits a pseudo-voigt fit to FLEET data using a nonlinear
%least squares fit.

%% Options
options = optimset(@lsqcurvefit);
options.Display = 'off';
options.TolFun = 1e-5;
options.MaxFunEvals = 1e3;
options.MaxIter = 200;
options.FinDiffType = 'central';

%% Limits
lim_UB = limits(1,:);
lim_g = limits(2,:);
lim_LB = limits(3,:);

%% Fitting a single row of data with a pseudo-voigt fit
        %fitlims = [h,n,x0,sigma,R];
        x = row_index;
    fit_pseudo_Voigt = @(v) v(1).*(v(2).*exp(-((x-v(3)).^2)./ (2*v(4)^2)) + (1-v(2)).*(((v(5)/2).^2)./((x-v(3)).^2+(v(5)/2).^2)));
    err_fit_gauss_g = @(v) fit_pseudo_Voigt(v)-row_intensity;
    %fit
    [fitvariables] = lsqnonlin(err_fit_gauss_g,lim_g,lim_LB,lim_UB,options);
    %subtract the fit from the rest of data processing

    %Plot Fit
    figure;
    plot(row_index,row_intensity);
    hold on;
    plot(row_index,fit_pseudo_Voigt(fitvariables));
    title('Fitting with pseudo-voigt')


end

