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



%g2 with fixed atom detuning
doD=linspace(-20,20);
dop = doD - I*gamma/2;
g = 10;
daD=[15 20 25 30];

figure;
for k=1:length(daD)
    daconst = daD(k);
    dap = daconst - I*k/2;
    C1gD = (eps*dop)./(g^2 - dap*dop);
    C2gD = (eps^2*((dap + dop).*dop + g^2))./(sqrt(2)*(dap*dop - g^2).*((dap + dop)*dap - g^2));
    g2D = (2*abs(C2gD).^2)./abs(C1gD).^4;
    lg2D = log10(g2D);
    plot(doD, lg2D)
    hold on;
end
yline(0,'--');
xlim([-20,10])
legend('\Delta_a = 15','\Delta_a = 20','\Delta_a = 25','\Delta_a = 30')

%g2 as function of detuning, weak coupling
daB=linspace(-1,1)';
doB=daB';
dap = daB - I*k/2;
dop = doB - I*gamma/2;
g = 1/sqrt(2);
C1gB = (eps*dop)./(g^2 - dap.*dop);
C2gB = (eps^2*((dap + dop).*dop + g^2))./(sqrt(2)*(dap.*dop - g^2).*((dap + dop).*dap - g^2));
g2B = (2*abs(C2gB).^2)./abs(C1gB).^4;
lg2B = log10(g2B);

daC=linspace(-2,2)';
doC=daC';
dap = daC - I*k/2;
dop = doC - I*gamma/2;
g = 1;
C1gC = (eps*dop)./(g^2 - dap.*dop);
C2gC = (eps^2*((dap + dop).*dop + g^2))./(sqrt(2)*(dap.*dop - g^2).*((dap + dop).*dap - g^2));
g2C = (2*abs(C2gC).^2)./abs(C1gC).^4;
lg2C = log10(g2C);

%plot g2
figure;
[c,h] = contourf(doA,daA,lg2A,100);
set(h, 'edgecolor','none');
xlabel('\Delta_0 / \gamma')
ylabel('\Delta_a / \gamma')
colorbar;

figure;
lg2B(lg2B<0)=0; %lg2B(lg2B<0)*0.1;
zmin=-2;
zmax=7;
zinc=0.001;                                      
zlev=zmin:zinc:zmax;
[cB,hB] = contourf(doB,daB,lg2B,zlev);
%hB.ZLocation = 0;
set(hB, 'edgecolor','none');
xlabel('\Delta_0 / \gamma')
ylabel('\Delta_a / \gamma')
colorbar;

figure;
lg2C(lg2C<0)=0;
zmin=0;
zmax=7;
zinc=0.001;                                      
zlev=zmin:zinc:zmax;
[cC,hC] = contourf(doC,daC,lg2C,zlev);
set(hC, 'edgecolor','none');
xlabel('\Delta_0 / \gamma')
ylabel('\Delta_a / \gamma')
colorbar;

%transmission
dc = 50:-0.5:-50;
g = 25.3;
k = 4.1;
gamma = 2.6;

figure;
detundiff = [0,11,23];
for k = 1:length(detundiff)
    da = detundiff(k) - dc;
    T = (abs((k*(gamma+I*da)./((k+I*-dc).*(gamma+I*da)+g^2)))).^2;
    plot(dc,T)
    hold on;
end
legend('0 MHz', '11 MHz', '23 MHz')
