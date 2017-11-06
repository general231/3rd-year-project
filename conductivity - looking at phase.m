global e
global hBar
global Ef
global tau
global thermalEnergy
e = 1.60217662e-19; %C charge of an electron
kB= 1.38064852e-23; % m^2 kg s^-2 K^-1 Boltzman constant
h = 6.62607004e-34; %J S Planck constant
hBar = h / (2*pi); %J s reducded planck constant
electronConcentration = 2.9e+12; % cm^-2
tau = 0.1e-12; %s carrier scattering timeat 293K
thermalEnergy = kB*293
Ef = (1.166e-7*sqrt(electronConcentration))*1.60218e-19 % Fermi Level
angularFrequency = logspace(13,15,10000)*2*pi;

intraConductivity = arrayfun(@intra, angularFrequency);
interConductivity = arrayfun(@inter, angularFrequency);
totalConductivity = intraConductivity + interConductivity;

subplot(2,2,1)
semilogx(angularFrequency,real(totalConductivity));
hold on
grid on
semilogx(angularFrequency,imag(totalConductivity));
legend('Real', 'Imaginary')
title('Subplot 1: Total Conductivity')
xlabel('Frequency(Hz)')
ylabel('Conductivity(S)')

subplot(2,2,2)
semilogx(angularFrequency, imag(intraConductivity))
hold on
grid on
semilogx(angularFrequency, imag(interConductivity))
semilogx(angularFrequency, abs(imag(intraConductivity)) - abs(imag(interConductivity)))
legend('Intra', 'Inter', 'differnece')
title('Subplot 2: Imaginary components of conductivity')
xlabel('Frequency(Hz)')
ylabel('Conductivity(S)')

subplot(2,2,3)
semilogx(angularFrequency, real(intraConductivity))
hold on
grid on
semilogx(angularFrequency, real(interConductivity))
semilogx(angularFrequency, abs(real(intraConductivity)) - abs(real(interConductivity)))
legend('Intra', 'Inter', 'differnece')
title('Subplot 3: Real components of conductivity')
xlabel('Frequency(Hz)')
ylabel('Phase(Radians)')

subplot(2,2,4)
semilogx(angularFrequency, radtodeg(angle(totalConductivity)))
grid on
title('Subplot 3: Phase of Conductivity')
xlabel('Frequency(Hz)')
ylabel('Phase(Degree)')

function conductivity = intra(w)
global e
global hBar
global Ef
global tau
global thermalEnergy

numerator = 2*i*e^2*thermalEnergy;

denominator = pi*hBar^2*(w+i/tau);
lnSection = log(2*cosh(Ef/(2*thermalEnergy)));

conductivity = (numerator/denominator)*lnSection;
end

function conductivity = inter(w)
global e
global hBar
global Ef
global thermalEnergy
hBarOmega = hBar*w;

if hBarOmega > 2*Ef;
    atanPart = atan((hBarOmega-2*Ef)/(2*thermalEnergy));
    numerator = (hBarOmega+2*Ef)^2;
    denominator = (hBarOmega-2*Ef)^2 + (2*thermalEnergy)^2;
    logPart = log(numerator/denominator);
    conductivity = ((e^2)/(4*hBar))*(0.5+(atanPart/pi)-((i/(2*pi))*logPart));
else
    conductivity = 0;
end
end

