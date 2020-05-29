function [] = f_clocktimer(time,char,bandwidth)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if(time>3600)
    fprintf([char ' %.2f hour(s),' '%.2f it/s \r'],time/3600,bandwidth);
elseif(time>60)
    fprintf([char ' %.2f minute(s),' '%.2f it/s \r'],time/60,bandwidth);
else
    fprintf([char ' %.2f second(s),' '%.2f it/s \r'],time,bandwidth);
end
end

