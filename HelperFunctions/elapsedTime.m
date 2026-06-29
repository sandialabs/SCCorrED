function tElapsedStr = elapsedTime(t)
    [h,m,s] = hms(seconds(t));
    tElapsedStr = "Elapsed time "+ num2str(h) + ":" + num2str(m) + ":" + num2str(s);
end