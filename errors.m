clear all;

% Graph N vs <n> for diff values of eta

% N=30;
% a=8;   %n_therm ('mean' photon number)
% n=0:N; %number of photons
% 
% %Distribution
% pn = a.^n ./ (1+a).^(n+1); 
% eta=1-cumsum(pn);

i=20;
a=0:i;

f1=figure(1);clf;set(f1,'units','centimeters','position',[0 3 12.5 6],...
    'color','w','paperpositionmode','auto');

ax(1)=axes('units','centimeters','position',[2 1 10 4.5],'box','on',...
    'linewidth',0.5,'fontsize',8,'nextplot','add');

xlabel(ax(1),'$\bar{n}$','fontsize',10,'interpreter','latex');
ylabel('N','fontsize',10,'interpreter','latex')

for eta=[0.001,0.01,0.1]
    N=log10(eta) ./ ( log10(a) - log10(1+a))   - 1; 
    hold on;
    box on;
    N(1)=1;
    plot(a,N)
 
 
end;

h = legend('$\mathcal{E}$=0.001','$\mathcal{E}$=0.01','$\mathcal{E}$=0.1',...
    'location','northwest');
set(h,'Interpreter','latex')
% legend('eta=0.1','eta=0.2','eta=0.3','eta=0.4','eta=0.5')