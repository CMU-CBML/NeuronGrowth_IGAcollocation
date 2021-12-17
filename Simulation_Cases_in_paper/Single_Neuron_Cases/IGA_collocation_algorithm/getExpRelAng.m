function [angle] = getExpRelAng(quantiles)

angle = johnsrnd(quantiles)*(2*(rand>.5)-1);
while abs(angle) >= quantiles(end)
    angle = johnsrnd(quantiles)*(2*(rand>.5)-1);
end

end

