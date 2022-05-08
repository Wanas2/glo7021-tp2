function out = gauss(x,s)
    out = exp(-(x.^2)/(2*s^2))./sqrt(2*pi*s^2);
end