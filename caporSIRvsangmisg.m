apw = 0.16;
aph = 0.08;
l1=1;
l2=1;
z=1;
txwidth = 0.06;
wavelength = 0.003;
% all descriptions from tx looking towards rx
% negative angles move tx to the left of the receiver
% tx one is on the left of tx2
phi11 = -pi/6;
phi12 = -pi/12;
phi13 = 0;
phi14 = pi/12;
phi15 = pi/6;
phi2space=linspace(-pi/4, pi/4, 7); % +
capacitys1 = getcapacitys(apw, aph, txwidth, phi11, phi2space, z, l1, l2, wavelength);
capacitys2 = getcapacitys(apw, aph, txwidth, phi12, phi2space, z, l1, l2, wavelength);
capacitys3 = getcapacitys(apw, aph, txwidth, phi13, phi2space, z, l1, l2, wavelength);
capacitys4 = getcapacitys(apw, aph, txwidth, phi14, phi2space, z, l1, l2, wavelength);
capacitys5 = getcapacitys(apw, aph, txwidth, phi15, phi2space, z, l1, l2, wavelength);

plot(phi2space, capacitys1, '-x', 'DisplayName', "phi1 = " + phi11);
hold on;
plot(phi2space, capacitys2, '-x','DisplayName', "phi1 = " + phi12);
plot(phi2space, capacitys3, '-x', 'DisplayName', "phi1 = " + phi13);
plot(phi2space, capacitys4, '-x', 'DisplayName', "phi1 = " + phi14);
plot(phi2space, capacitys5, '-x', 'DisplayName', "phi1 = " + phi15);
legend;

xlabel('Phi2 (angle of incidence from tx2)'); % Label for x-axis
ylabel('average Capacity at reception with steering and demodulation'); % Label for y-axis
title('Average Capacity vs angular misalignment'); % Title for the plot



%% helper function get capacitys

function capacitys = getcapacitys(apw, aph, txwidth, phi1, phi2space, z, l1, l2, wavelength, s)
capacitys = zeros(7,1);
for i=1:7
    
    phi2 = phi2space(i);
    
[Xs1, axs1, Ys1, ays1, afsz1, Xs2, axs2, Ys2, ays2, afsz2, Xe1, Ye1, Xe2, Ye2, ...
   matrixehalfxsize, matrixeysize, range1Center, range2Center, ...
   xeRange, yeRange] = generateAps(apw, aph, txwidth, phi1, phi2, z);

huygenSource1 = exp(1j * l1 * atan2(Ys1, Xs1-range1Center));
huygenSource2 = exp(1j * l2 * atan2(Ys2, Xs2-range2Center));

RxEfieldCat11 = hf(matrixeysize, matrixehalfxsize, axs1, Xe1, ays1, Ye1, afsz1, huygenSource1, wavelength);
RxEfieldCat12 = hf(matrixeysize, matrixehalfxsize, axs2, Xe1, ays2, Ye1, afsz2, huygenSource2, wavelength);
RxEfieldCat1 = RxEfieldCat11 + RxEfieldCat12;
RxEfieldCat22 = hf(matrixeysize, matrixehalfxsize, axs2, Xe2, ays2, Ye2, afsz2, huygenSource2, wavelength);
RxEfieldCat21 = hf(matrixeysize, matrixehalfxsize, axs1, Xe2, ays1, Ye2, afsz1, huygenSource1, wavelength);
RxEfieldCat2 = RxEfieldCat21 + RxEfieldCat22;

SteerCat11 = steer(phi1, RxEfieldCat11, Xe1, wavelength);
SteerCat12 = steer(phi1, RxEfieldCat12, Xe1, wavelength);
SteerCat21 = steer(phi2, RxEfieldCat21, Xe2, wavelength);
SteerCat22 = steer(phi2, RxEfieldCat22, Xe2, wavelength);
DemodSCat11 = SteerCat11.*exp(-1j*l1*atan2(Ye1, Xe1-range1Center));
DemodSCat12 = SteerCat12.*exp(-1j*l1*atan2(Ye1, Xe1-range1Center));
DemodSCat21 = SteerCat21.*exp(-1j*l2*atan2(Ye2, Xe2 - range2Center));
DemodSCat22 = SteerCat22.*exp(-1j*l2*atan2(Ye2, Xe2 - range2Center));

p11 = sum((abs(DemodSCat11)).^2, "all");
p12 = sum((abs(DemodSCat12)).^2, "all");
p21 = sum((abs(DemodSCat21)).^2, "all");
p22 = sum((abs(DemodSCat22)).^2, "all");

sir1 = p11/p12;
sir2 = p22/p21;

capacity1 = log2(1+sir1);
capacity2 = log2(1+sir2);

capacitys(i,1) = (capacity1+capacity2)/2;

end

end
 
%% helper function generate coords
function [Xs1, axs1, Ys1, ays1, afsz1, Xs2, axs2, Ys2, ays2, afsz2, Xe1, Ye1, Xe2, Ye2, matrixehalfxsize, matrixeysize, range1Center, range2Center, xeRange, yeRange]=generateAps(apwid, aph, txwidth, phi1, phi2, z)
    
    xeRangeCenter = 0;
    yeRangeCenter = 0;
    xeRangeWidth = apwid;
    yeRangeWidth = aph;
    matrixexSize = ceil(apwid*1500/10)*10;
    matrixeySize = ceil(aph*1500/10)*10;
    xeRange = linspace(xeRangeCenter-xeRangeWidth/2, xeRangeCenter+xeRangeWidth/2, matrixexSize);
    yeRange = linspace(yeRangeCenter-yeRangeWidth/2, yeRangeCenter+yeRangeWidth/2, matrixeySize);

    xe1RangeCenter = -1*xeRangeWidth/4;
    ye1RangeCenter = 0;
    xe1RangeWidth = apwid/2;
    ye1RangeWidth = aph;
    matrixe1xSize = ceil(xe1RangeWidth*1500/10)*10;
    matrixe1ySize = matrixeySize;
    xe1Range = linspace(xe1RangeCenter-xe1RangeWidth/2, xe1RangeCenter+xe1RangeWidth/2, matrixe1xSize);
    ye1Range = linspace(ye1RangeCenter-ye1RangeWidth/2, ye1RangeCenter+ye1RangeWidth/2, matrixe1ySize);
    [Xe1, Ye1] = meshgrid(xe1Range, ye1Range);

    xe2RangeCenter = xeRangeWidth/4;
    ye2RangeCenter = 0;
    xe2RangeWidth = apwid/2;
    ye2RangeWidth = aph;
    matrixe2xSize = matrixe1xSize;
    matrixe2ySize = matrixeySize;
    xe2Range = linspace(xe2RangeCenter-xe2RangeWidth/2, xe2RangeCenter+xe2RangeWidth/2, matrixe2xSize);
    ye2Range = linspace(ye2RangeCenter-ye2RangeWidth/2, ye2RangeCenter+ye2RangeWidth/2, matrixe2ySize);
    [Xe2, Ye2] = meshgrid(xe2Range, ye2Range);

    xs1RangeCenter = -1*xeRangeWidth/4;
    ys1RangeCenter = 0;
    xs1RangeWidth = txwidth;
    ys1RangeWidth = xs1RangeWidth;
    matrixs1Size = ceil(txwidth*1500/10)*10;
    xs1Range = linspace(xs1RangeCenter-xs1RangeWidth/2, xs1RangeCenter+xs1RangeWidth/2, matrixs1Size);
    ys1Range = linspace(ys1RangeCenter-ys1RangeWidth/2, ys1RangeCenter+ys1RangeWidth/2, matrixs1Size);
    [Xs1, Ys1] = meshgrid(xs1Range, ys1Range);

    [axs1, ays1, afsz1] = findActualCoords(phi1, 0, z, Xs1, Ys1, xs1RangeCenter);

    xs2RangeCenter = xeRangeWidth/4;
    ys2RangeCenter = 0;
    xs2RangeWidth = txwidth;
    ys2RangeWidth = xs2RangeWidth;
    matrixs2Size = ceil(txwidth*1500/10)*10;
    xs2Range = linspace(xs2RangeCenter-xs2RangeWidth/2, xs2RangeCenter+xs2RangeWidth/2, matrixs2Size);
    ys2Range = linspace(ys2RangeCenter-ys2RangeWidth/2, ys2RangeCenter+ys2RangeWidth/2, matrixs2Size);
    [Xs2, Ys2] = meshgrid(xs2Range, ys2Range);

    [axs2, ays2, afsz2] = findActualCoords(phi2, 0, z, Xs2, Ys2, xs2RangeCenter);

    range1Center = xs1RangeCenter;
    range2Center = xs2RangeCenter;
    matrixehalfxsize = matrixe1xSize;
    matrixeysize = matrixeySize;
end
%% helper function find actual coords
function [axs, ays, afsz] = findActualCoords(phi, phi2, azcenter, Xs, Ys, xcenter)
   fszc = azcenter*cos(phi)*cos(phi2); % distance from center of tx to plane of rx
   % Been Checked
   fsz = (Xs-xcenter) * sin(phi) + Ys * cos(phi) * sin(phi2);
   % Been Checked
   afsz = -1*fsz + fszc;
   xOffset = azcenter*sin(phi)*cos(phi2) + xcenter; % distance from tx x center to rx x center
   % Been Checked
   axs = xOffset + (Xs-xcenter)*cos(phi) + Ys*sin(phi)*sin(phi2);
   % Been Checked
   yOffset = azcenter*sin(phi2);
   % Been Checked
   ays = yOffset + Ys*cos(phi2);
end
%% hugyens fresnel
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
%% helper function plotPhase
function plotPhase(efield, xrange, yrange, fieldname)
   phase = angle(efield);
   imagesc(xrange, yrange, phase);
   colormap(jet);
   colorbar;
   xlabel('X Position (m)'); % Label for x-axis
   ylabel('Y Position (m)'); % Label for y-axis
   title('Phase of ' + fieldname); % Title for the plot
end
%% helper function plotIntensity
function plotIntensity(efield, xrange, yrange, fieldname)
   intensity = abs(efield).^2;
   imagesc(xrange, yrange, intensity);
   colormap(jet);
   colorbar;
   xlabel('X Position (m)'); % Label for x-axis
   ylabel('Y Position (m)'); % Label for y-axis
   title('Intensity of ' + fieldname); % Title for the plot
end
%% helper function beam steering
function steeringplatesurfacehs = steer(phi, Efield, Xe, wavelength)
 steering1aphase = -sin(phi)*2*pi*Xe/wavelength;
 steeringplatesurfacehs = Efield.*exp(-1j*steering1aphase);
end

%% helper function find purity

function purity = findPurity(Efield, l, Xe, Ye, xeRangeCenter)

    [Nr, Nphi] = size(Efield);
    r_max = sqrt((Xe(end)-xeRangeCenter)^2 + Ye(end)^2);
    r = linspace(0, r_max, Nr);
    phi = linspace(0, 2*pi, Nphi);
    [R, PHI] = meshgrid(r, phi);

    x_interp = R .* cos(PHI);
    y_interp = R .* sin(PHI);

    E_field_circ = interp2(Xe-xeRangeCenter, Ye, Efield, x_interp, y_interp, 'spline').';
    
    lrange = linspace(-20,20, 41);

    Al = 0;
    Alm = 0;
    for m = 1:41
        oam = lrange(m);
        cl = zeros(Nr, 1);
        % Issue is E field is in cartesian coordinates
        for r = 1:Nr
            azimuth = E_field_circ(r,:);
            cl(r) = 1/Nphi * sum(azimuth .* exp(-1j * oam * (0:Nphi-1) * (2*pi/Nphi)));
        end
        Alm = Alm + sum(abs(cl).^2);
        if oam == l 
            Al = sum(abs(cl).^2);
        end
    end
    purity = Al/Alm;
end   

