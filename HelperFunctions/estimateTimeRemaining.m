function tRemainingStr = estimateTimeRemaining(total, count, t)
    tPerPoint = t/count;
    pointsRemaining = total - count;
    tRemaining = pointsRemaining*tPerPoint;
    [h,m,s] = hms(seconds(tRemaining));
    tRemainingStr = "Estimate time remaining "+ num2str(h) + ":" + num2str(m) + ":" + num2str(s);
end
