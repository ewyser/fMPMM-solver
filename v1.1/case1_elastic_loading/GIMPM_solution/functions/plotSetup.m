function [pp] = plotSetup(xinf,xsup,yinf,ly,y,BCy)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
pp.LimX  =[xinf xsup]                                                 ;%
pp.LimY  =[yinf 2*ly]                                                 ;%
pp.ps    =3.0                                                         ;%
pp.BoundX=pp.LimX                                                     ;%
pp.BoundY=[y(BCy(1)) 10]                                          ;%
pp.cbclass = 20                                                       ;%
end

