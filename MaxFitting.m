function [fitvariables,outfit,fit_unc] = MaxFitting(row_intensity,row_index,limits)
%This function fits a pseudo-voigt fit to FLEET data using a nonlinear
%least squares fit.

%% Options
options = optimset(@lsqcurvefit);
options.Display = 'off';
options.TolFun = 1e-6;
options.MaxFunEvals = 1e3;
options.MaxIter = 400;
options.FinDiffType = 'central';

%% Limits
lim_UB = limits(1,:);
lim_g = limits(2,:);
lim_LB = limits(3,:);

%% Fitting a single row of data with a pseudo-voigt fit
        %fitlims = [h,n,x0,sigma,R,bkg];
        x = row_index;
    fit_pseudo_Voigt = @(v) v(1).*(v(2).*exp(-((x-v(3)).^2)./ (2*v(4)^2)) + (1-v(2)).*(((v(5)/2).^2)./((x-v(3)).^2+(v(5)/2).^2)))+v(6);
    err_fit_gauss_g = @(v) fit_pseudo_Voigt(v)-row_intensity;
    %fit
    [fitvariables,~,residual,~,~,~,jacobian] = lsqnonlin(err_fit_gauss_g,lim_g,lim_LB,lim_UB,options);

    var_err = nlparci(fitvariables,residual,'jacobian',jacobian);   %  95% confidence intervals for the fit coefficients
    fit_unc = (var_err(3,2)-var_err(3,1))/2;

%     %Plot Fit
%     figure;
%     plot(row_index,row_intensity);
%     hold on;
%     plot(row_index,fit_pseudo_Voigt(fitvariables));
%     title('Fitting with pseudo-voigt')

%     Fit to send out
    outfit = fit_pseudo_Voigt(fitvariables);

end

