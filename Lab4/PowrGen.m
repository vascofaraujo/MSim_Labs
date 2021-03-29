function [Powr, an] = PowrGen(spos,npos)
D = squareform(pdist([spos zeros(size(spos)) npos]'));
d = D(1,3:end);				% Source-anchor distances
an = D(2,3:end);			% Anchor norms

% Generate observations
Pow0 = 100;				% Source power
Powr = Pow0./(d.^2);				% Noiseless RSSI

stdev = 1e-1;				% Log-noise standard deviation
%stdev = 0;
Powr = Powr.*exp(stdev*randn(size(Powr)));	% Introduce noise
QP = 1e-2;
Powr = QP*round(Powr/QP);			% Quantize power measurements
end

