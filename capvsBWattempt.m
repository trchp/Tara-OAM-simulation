%% SPP: bandiwth vs capacity
 
xsRangeCenter = 0;
ysRangeCenter = xsRangeCenter;
xsRangeWidth = 0.06;
ysRangeWidth = xsRangeWidth;
matrixsSize = 80;
xsRange = linspace(xsRangeCenter - xsRangeWidth/2, xsRangeCenter + xsRangeWidth/2, matrixsSize);
ysRange = linspace(ysRangeCenter - ysRangeWidth/2, ysRangeCenter + ysRangeWidth/2, matrixsSize);
[Xs, Ys] = meshgrid(xsRange, ysRange);

l=1;

xeRangeCenter = 0;
yeRangeCenter = xeRangeCenter;
xeRangeWidth = 0.08;
yeRangeWidth = xeRangeWidth;
matrixeSize = 100;
xeRange = linspace(xeRangeCenter - xeRangeWidth/2, xeRangeCenter + xeRangeWidth/2, matrixeSize);
yeRange = linspace(yeRangeCenter - yeRangeWidth/2, yeRangeCenter + yeRangeWidth/2, matrixeSize);
[Xe, Ye] = meshgrid(xeRange, yeRange);

wavelength = .002; % meters
fc = 3*10^8/wavelength;
k = 2*pi/wavelength;
totalz = 1; % meters
n=1.6;

stepheight = wavelength/(n-1);
height = mod((atan2(Ys,Xs)*l), 2*pi) * stepheight / (2*pi);


bws = linspace(3000*10^9, 5000*10^9, 6);

capacitys = zeros(6, 1);

for i = 1:6
    bw = bws(i);
    frequencys = linspace(fc-bw/2, fc+bw/2, bw/10^9/500);
    captot = 0;
    for j = 1:bw/10^9/500
        frequency = frequencys(j);
        wl = 3*10^8/frequency;
        phase = 2*height*pi*(n-1)/wl;
        huygenSource = exp(1*j*phase);
        RxEfield1 = hf(matrixeSize, matrixeSize, Xs, Xe, Ys, Ye, totalz, huygenSource, wavelength);
        RxEfield1n = RxEfield1/max(RxEfield1, [], 'all');
        RxEfield = RxEfield1n.*exp(-1*j*l*atan2(Ye, Xe));
        s = 10; % dB
        RfieldHat = awgn(RxEfield, s);
            %0 db/pn = 50 db
            %pn = 0db/50db
            %10*log10(1)/10*log10(100,000);
        npraw = 1/(10^(s/10));
        sp = sum(abs(RfieldHat).^2, "all");
        np = abs(npraw)^2;
        snr = sp/np;
        capf = log2(1+snr);
        captot = captot+capf;
    end
    capacitys(i) = captot;
end

plot( bws, capacitys);


%% helper function hf

function Ee = hf(matrixeySize, matrixexSize, Xs, Xe, Ys, Ye, z, huygenSource, wavelength)
   Ee = zeros(matrixeySize, matrixexSize);
   k = 2 * pi/wavelength;
   for i = 1:matrixeySize
       for j = 1:matrixexSize
           allRs = ((Xs-Xe(i, j)).^2 + (Ys-Ye(i, j)).^2 + z.^2).^0.5;
           contributions = huygenSource.*exp(1j*k*allRs)./allRs;
           hfSum = sum(contributions, "all");
           Ee(i, j) = hfSum*(1/(1j*wavelength));
       end
   end
end
