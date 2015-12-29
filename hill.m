function hill_ll = hill(H,J,theta,P_hit_range)
hill_ll = 0;
for i = 1:length(P_hit_range)
    hill_ll = hill_ll + (H(i) - J(i))*log(1 - theta(1)*P_hit_range(i)^theta(3)/...
        (P_hit_range(i)^theta(3)+theta(2)^theta(3))) + J(i)*log(theta(1)*...
        P_hit_range(i)^theta(3)/(P_hit_range(i)^theta(3)+theta(2)^theta(3)));
end
hill_ll = -1*hill_ll;
end
