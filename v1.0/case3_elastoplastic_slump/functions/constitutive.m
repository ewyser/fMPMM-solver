function [mpD] = constitutive(mpD,Del,Gc,Kc,nu,plasticity,dt,it,te)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% ELASTIC PREDICTOR: TRIAL ELASTIC STEP
[mpD] = elastic(mpD,Del)                                                  ;% - see function
%--------------------------------------------------------------------------%
%mpD.s(3,:) = nu.*(mpD.s(1,:)+mpD.s(2,:));
%% PLASTIC CORRECTION: EXACT SOLUTION
if((plasticity)&&(dt*it>te))
    %[mpD,dat] = plasticDP(mpD,Gc,Kc)                                  ;% - see function
    [mpD,D] = plastic(mpD,Del)                                            ;% - see function
end
%--------------------------------------------------------------------------%

end

