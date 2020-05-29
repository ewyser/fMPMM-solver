function [fext] = getFext(mpD,meD,g,l2g)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% CONTRIBUTION FROM p TO N
fe   = reshape([mpD.S.*0.0;mpD.S.*repmat(mpD.m,1,meD.nNe).*-g ],mpD.n*meD.nDoF(1),1)              ;% preprocesing of external force global vector
fext = accumarray(l2g(:),fe,[meD.nDoF(2) 1])                               ;% external force global vector
end

