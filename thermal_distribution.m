clear all;

% Graphs thermal photon number distribution with mean photon number 'a' and
% N photons
figure(1);clf;

N=30;
a=8;
n=0:N;



pn = a.^n ./ (1+a).^(n+1); % probability to find photon in field

bar(n,pn);
title('Thermal photon number distribution')
xlabel('Photon number N')
ylabel('Probability')

hold on;

%Rescale
Pn=pn./sum(pn);

sum(Pn)

bar(n,Pn,'facecolor','r')
bar(n,pn);
% end;
% 
% clear all;
% % Number distribution over Renorm Factor A = 1 - error , also calculates
% % fidelity
% figure(2)
% for N=1:30
%     a=5;
%     x(N)=N;
%     p(N) = a.^N / (1+a).^(N+1);  % probability to find photon in field
%     eta(N) = (a.^(N+1) / (1+a).^(N+1)); %error in distribution
%     A(N) = 1 - eta(N);
%     y(N) = p(N) / A(N); %scaled probability
%     
%     Fid(N)=sqrt(A(N));  % Output fidelity
%     
%     bar(x,y);
%     title('Adjusted photon number distribution')
%     xlabel('Photon number N')
%     ylabel('Probability/A')
% end;




