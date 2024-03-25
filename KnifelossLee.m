function Gd = KnifelossLee(v)

% Gd = 20 * log10(v);
Gd = 6.9 + 20 * log10(     ((v - 0.1)^2 +1) ^0.5 +v-0.1        );


end