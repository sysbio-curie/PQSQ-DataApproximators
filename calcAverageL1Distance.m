function L1av = calcAverageL1Distance(x,xprojected)

ax = abs(x-xprojected);
ax = ax';
L1av = mean(sum(ax));

end