function [sig] = elasticTrial(mpD,sig0,Del)
%f_elastic : function to calculate 1) a rotated stress to ensure objectivity
%of the stress tensor and, 2) an elastic predictor
%   Detailed explanation goes here
%% PREPROCESSING
omega      = -mpD.w(:,3)'                                                 ;% incremental rotation in (x,y) plane
cos2       = cos(omega).^2                                                ;% -
sin2       = sin(omega).^2                                                ;% -
sincos     = sin(omega).*cos(omega)                                       ;% -
%--------------------------------------------------------------------------%

%% OBJECTIVE STRESS ROTATION DUE TO NON-OBJECTIVE STESS RATE ds_{ij}/dt
sig0(1,:) = sig0(1,:).*cos2 + ...
            sig0(2,:).*sin2 - 2.*sig0(end,:).*sincos                      ;% primed sig_xx
sig0(2,:) = sig0(1,:).*sin2 + ...
            sig0(2,:).*cos2 + 2.*sig0(end,:).*sincos                      ;% primed sig_yy
sig0(4,:) =(sig0(1,:)-sig0(2,:)).*sincos + ...
            sig0(end,:).*(cos2-sin2)                                      ;% primed sig_xy/sig_yx       
%--------------------------------------------------------------------------%

%% STRESS-STRAIN CONSTITUTIVE RELATION
sig     = sig0+(Del*mpD.e)                                                ;% elastic constitutive relation
%--------------------------------------------------------------------------%

end

