function[PositivePly, NegativePly] = fibreindex(Phi,n)

PositivePly = [sum(Phi > 0); sum(Phi > 0)/n];
NegativePly = [sum(Phi < 0); sum(Phi < 0)/n];

end