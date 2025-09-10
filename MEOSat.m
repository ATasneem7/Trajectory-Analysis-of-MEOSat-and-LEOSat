%% MEOSAT Capturing the Dynamics of Orbital Motion in a User-Defined Function "Satellite"

function s_dot = MEOSat(t, s, go, R_E)
s_dot = zeros(4, 1);
s_dot(1) = s(3);
s_dot(2) = s(4);
s_dot(3) = (-go*R_E^2)/(norm(s(1:2))^3)*s(1);
s_dot(4) = (-go*R_E^2)/(norm(s(1:2))^3)*s(2);
end