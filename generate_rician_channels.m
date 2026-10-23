function [hBI_max,hIA_max,hEI_max,GBI,GIA,GEI] = generate_rician_channels( ...
    Mb,Ma,Me,P,NI_max,K,rho,LS,locB,locA,locI,locE)

if LS == 1
    [~,GBI,GIA] = func_path_gain(locB,locA,locI);
    [~,GEI,~] = func_path_gain(locE,locA,locI);
else
    GBI = 1;
    GIA = 1;
    GEI = 1;
end

aLoS = sqrt(K/(1+K));
aNLoS = sqrt(1/(1+K));

hBI_LoS = exp(1i*2*pi*rand(Mb,NI_max));
hBI_NLoS = (randn(Mb,NI_max)+1i*randn(Mb,NI_max))/sqrt(2);
hBI_max = sqrt(GBI)*(aLoS*hBI_LoS+aNLoS*hBI_NLoS);

hIA_LoS = exp(1i*2*pi*rand(NI_max,Ma));
hIA_NLoS = (randn(NI_max,Ma)+1i*randn(NI_max,Ma))/sqrt(2);
hIA_max = sqrt(GIA)*(aLoS*hIA_LoS+aNLoS*hIA_NLoS);

hEI_LoS = exp(1i*2*pi*rand(Me,NI_max));
innovation = (randn(Me,NI_max)+1i*randn(Me,NI_max))/sqrt(2);
hBI_shared = repmat(hBI_NLoS,Me,1);
hEI_NLoS = rho*hBI_shared+sqrt(1-rho^2)*innovation;
hEI_max = sqrt(GEI)*(aLoS*hEI_LoS+aNLoS*hEI_NLoS);
end
