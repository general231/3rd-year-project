global e
global hBar
global Ef
global tau
global thermalEnergy
e = 1.60217662e-19; %C charge of an electron
kB= 1.38064852e-23; % m^2 kg s^-2 K^-1 Boltzman constant
h = 6.62607004e-34; %J S Planck constant
hBar = h / (2*pi); %J s reducded planck constant
electronConcentration = 1.2e+13; % cm^-2
tau = 0.5e-12; %s carrier scattering timeat 293K
thermalEnergy = kB*293
Ef = (1.166e-7*sqrt(electronConcentration))*1.60218e-19 % Fermi Level
angularFrequency = logspace(8,15,10000)*2*pi;

intraConductivity = arrayfun(@intra, angularFrequency);
interConductivity = arrayfun(@inter, angularFrequency);
totalConductivity = intraConductivity + interConductivity;
cond = 4*e*e/hBar
subplot(2,2,1)
semilogx(angularFrequency,real(intraConductivity));
hold on
grid on
semilogx(angularFrequency,imag(intraConductivity));
semilogx(angularFrequency,abs(intraConductivity));
legend('Real', 'Imaginary', 'Magnitude')
title('Subplot 1: Intraband Conductivity')
xlabel('Frequency(Hz)')
ylabel('Conductivity(S)')

subplot(2,2,3)
semilogx(angularFrequency,real(interConductivity));
hold on
grid on
semilogx(angularFrequency,imag(interConductivity));
semilogx(angularFrequency,abs(interConductivity));
legend('Real', 'Imaginary', 'Magnitude')
title('Subplot 2: Interband Conductivity')
xlabel('Frequency(Hz)')
ylabel('Conductivity(S)')

subplot(2,2,2)
semilogx(angularFrequency,real(totalConductivity));
hold on
grid on
semilogx(angularFrequency,imag(totalConductivity));
semilogx(angularFrequency,abs(totalConductivity));
legend('Real', 'Imaginary', 'Magnitude')
title('Subplot 3: Total Conductivity')
xlabel('Frequency(Hz)')
ylabel('Conductivity(S)')

subplot(2,2,4)
semilogx(angularFrequency, radtodeg(angle(totalConductivity)))
grid on
legend('Total')
title('Subplot 4: Phase')
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

if true;
    atanPart = atan((hBarOmega-2*Ef)/(2*thermalEnergy));
    numerator = (hBarOmega+2*Ef)^2;
    denominator = (hBarOmega-2*Ef)^2 + (2*thermalEnergy)^2;
    logPart = log(numerator/denominator);
    conductivity = ((e^2)/(4*hBar))*(0.5+(atanPart/pi)-((i/(2*pi))*logPart));
else
    conductivity = 0;
end
end

