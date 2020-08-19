function new = decimal(num, Ndecimals)

% This function round 'num' to 'Ndecimals' points

f=10.^Ndecimals;
new = round(f*num)/f;

end