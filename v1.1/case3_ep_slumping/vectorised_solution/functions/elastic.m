function [mpD] = elastic(mpD,Del)
%f_elastic : function to calculate 1) a rotated stress to ensure objectivity
%of the stress tensor and, 2) an elastic predictor
%   Detailed explanation goes here
%% PREPROCESSING
omega      = -mpD.w(:,3)'                                                 ;% incremental rotation in (x,y) plane
cos2       = cos(omega).^2                                                ;% -
sin2       = sin(omega).^2                                                ;% -
sincos     = sin(omega).*cos(omega)                                       ;% -
old_s      = mpD.s                                                        ;%
% OBJECTIVE STRESS ROTATION DUE TO NON-OBJECTIVE STESS RATE ds_{ij}/dt
mpD.s(1,:) = old_s(1,:).*cos2 + ...
             old_s(2,:).*sin2 - 2.*old_s(3,:).*sincos                     ;% primed sig_xx
mpD.s(2,:) = old_s(1,:).*sin2 + ...
             old_s(2,:).*cos2 + 2.*old_s(3,:).*sincos                     ;% primed sig_yy
mpD.s(3,:) =(old_s(1,:)-old_s(2,:)).*sincos + ...
             old_s(3,:).*(cos2-sin2)                                      ;% primed sig_xy/sig_yx       
%--------------------------------------------------------------------------%

%% STRESS-STRAIN CONSTITUTIVE RELATION
mpD.s      = mpD.s+(Del*mpD.e)                                            ;% elastic constitutive relation
clear omega cos2 sin2 sincos                                              ;% clear temporary variables  
%--------------------------------------------------------------------------%

end

