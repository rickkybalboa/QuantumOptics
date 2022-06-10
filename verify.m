clear all;

% Graphs thermal photon number distribution with mean photon number 'a' and
% N photons
figure(1);clf;

N=30;
a=8;
n=0:N;
%Distribution
pn = a.^n ./ (1+a).^(n+1); 
eta=1-cumsum(pn);

hold on;
bar(n,eta)
plot(n,a.^(n+1)./(1+a).^(n+1),'rx')


% By comparing eta = 1 - cumsum(pn)  (which we know is correct)
% to our sum expression for eta(derived from wolfram) we can verify 
% whether the latter expression is correct.