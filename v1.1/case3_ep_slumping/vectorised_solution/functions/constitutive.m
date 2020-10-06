function [mpD] = constitutive(mpD,Del,plasticity,dt,it,te)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% ELASTIC PREDICTOR: TRIAL ELASTIC STEP
[mpD] = elastic(mpD,Del)                                                  ;% - see function
%--------------------------------------------------------------------------%

%% PLASTIC CORRECTION: EXACT SOLUTION
if((plasticity)&&(dt*it>te))
    [mpD,dat] = plastic(mpD,Del)                                  ;% - see function
end
%--------------------------------------------------------------------------%

end

