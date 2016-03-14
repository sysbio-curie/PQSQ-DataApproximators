function inds = splitPointIntervalsFast1(x,intervals)

[A, inds] = histc(abs(x),intervals);