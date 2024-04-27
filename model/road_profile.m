function zr = road_profile(t, p)
% quarter car model with trailing arm suspension on random road

% pseudo random road input
s = p.v*t; 
u = sum(p.Amp.*sin(p.Om*s+p.Psi)); 
zr = u;

end