k=1;
gamma=k;
eps = 0.01;
I = sqrt(-1);

%g2 as function of detuning, strong coupling
daA=linspace(-40,40)';
doA=daA';
dap = daA - I*k/2;
dop = doA - I*gamma/2;
g = 10;
C1gA = (eps*dop)./(g^2 - dap.*dop);
C2gA = (eps^2*((dap + dop).*dop + g^2))./(sqrt(2)*(dap.*dop - g^2).*((dap + dop).*dap - g^2));
g2A = (2*abs(C2gA).^2)./abs(C1gA).^4;
lg2A = log10(g2A);

%g2 as function of detuning, weak coupling
daB=linspace(-1,1)';
doB=daB';
dap = daB - I*k/2;
dop = doB - I*gamma/2;
g = 1/sqrt(2);
C1gB = (eps*dop)./(g^2 - dap.*dop);
C2gB = (eps^2*((dap + dop).*dop + g^2))./(sqrt(2)*(dap.*dop - g^2).*((dap + dop).*dap - g^2));
g2B = (2*abs(C2gB).^2)./abs(C1gB).^4;
lg2B = (g2B);

daC=linspace(-2,2)';
doC=daC';
dap = daC - I*k/2;
dop = doC - I*gamma/2;
g = 1;
C1gC = (eps*dop)./(g^2 - dap.*dop);
C2gC = (eps^2*((dap + dop).*dop + g^2))./(sqrt(2)*(dap.*dop - g^2).*((dap + dop).*dap - g^2));
g2C = (2*abs(C2gC).^2)./(abs(C1gC).^4);
lg2C = (g2C);

%plot g2
figure;
[c,h] = contourf(doA,daA,lg2A,100);
hold on;
plot(linspace(-40,-0.5), 10^2./linspace(-40,-0.5),'w--')
plot(linspace(0.5,40), 10^2./linspace(0.5,40),'w--')
plot(doA,doA,'k')
hold off;
set(h, 'edgecolor','none');
xlabel('\Delta_0 / \gamma')
ylabel('\Delta_a / \gamma')
title('log_{10} g^2(0), g=10\gamma')
xlim([-40,40])
ylim([-40,40])
colorbar;

figure;
lg2B(lg2B<0)=0; 
zmin=-2;
zmax=7;
zinc=0.001;                                      
zlev=zmin:zinc:zmax;
[cB,hB] = contourf(doB,daB,lg2B,1000);
set(hB, 'edgecolor','none');
title('g^2(0), g=2\surd{2}\gamma')
xlabel('\Delta_0 / \gamma')
ylabel('\Delta_a / \gamma')
yticks([-1 -0.5 0 0.5 1])
xticks([-1 -0.5 0 0.5 1])
colorbar;

figure;
lg2C(lg2C<0)=0;
zmin=0;
zmax=7;
zinc=0.001;                                      
zlev=zmin:zinc:zmax;
[cC,hC] = contourf(doC,daC,lg2C,200);
set(hC, 'edgecolor','none');
title('g^2(0), g=1\gamma')
xlabel('\Delta_0 / \gamma')
ylabel('\Delta_a / \gamma')
yticks([-2 -1 0 1 2])
xticks([-2 -1 0 1 2])
colorbar;

%transmission
% dc = 50:-0.5:-50;
% g = 25.3;
% k = 4.1;
% gamma = 2.6;
% 
% figure;
% detundiff = [0,11,23];
% for k = 1:length(detundiff)
%     da = detundiff(k) - dc;
%     T = (abs((k*(gamma+I*da)./((k+I*-dc).*(gamma+I*da)+g^2)))).^2;
%     plot(dc,T)
%     hold on;
% end
% legend('0 MHz', '11 MHz', '23 MHz')
