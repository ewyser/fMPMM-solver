function [vid1] = video(title)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    date    = datevec(now);
    name_im = [num2str(date(1)) '_' num2str(date(2)) '_' num2str(date(3)) '_' num2str(round(date(4))) 'h' num2str(date(5)) 'm' num2str(round(date(6))) 's' ];
    vid1         = VideoWriter(['./data/' title char '.avi']);
    vid1.Quality = 95;
    open(vid1);
end

