function p = nUpdateWaitbar(data, h, initialCount)
    persistent TOTAL COUNT H START
    if nargin == 2
        H = h;
        TOTAL = data;
        COUNT = 0;
        START = 0;
    elseif nargin == 3
        H = h;
        TOTAL = data;
        COUNT = initialCount;
        START = initialCount;
    else
        COUNT = 1 + COUNT;  
        p = COUNT/TOTAL;
        t = toc;
        tElap = elapsedTime(t);
        tRemain = estimateTimeRemaining(TOTAL - START, COUNT - START, t);
        waitStr = [num2str(COUNT) + "/" + num2str(TOTAL) + " | " + tElap, tRemain];
%         waitStr = "Completed points " + num2str(COUNT) + ...
%             ", elapsed time " + num2str(h) + ":" + num2str(m) + ...
%             ":" + num2str(s);
        waitbar(p, H, waitStr)
    end
end


function tElapsedStr = elapsedTime(t)
    [h,m,s] = hms(seconds(t));
    tElapsedStr = "Elapsed time "+ num2str(h) + ":" + num2str(m) + ":" + num2str(s);
end


function tRemainingStr = estimateTimeRemaining(total, count, t)
    tPerPoint = t/count;
    pointsRemaining = total - count;
    tRemaining = pointsRemaining*tPerPoint;
    [h,m,s] = hms(seconds(tRemaining));
    tRemainingStr = "Estimate time remaining " + num2str(h) + ":" + num2str(m) + ":" + num2str(s);
end