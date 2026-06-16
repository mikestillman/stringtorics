mysterySurface = (S) -> (
     -- S should be a poly ring in 5 vars
     -- The ideal we will return is the embedding via quartics thru
     -- 10 points of P^2.
     kk := coefficientRing S;
     P2 = kk[x,y,z];
     pt1 := random(P2^10, P2^3);
     I = intersect apply(0..9, i -> trim minors(2,(vars P2) || pt1^{i}));
     trim ker map(P2,S,gens I)     
     )

end--
mysterySurface(S)


