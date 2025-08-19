% angular misalignment animation 

loops = 51;


xsRangeCenter = 0;
ysRangeCenter = xsRangeCenter;
xsRangeWidth = 0.06;
ysRangeWidth = xsRangeWidth;
matrixsSize = 80;
xsRange = linspace(xsRangeCenter - xsRangeWidth/2, xsRangeCenter + xsRangeWidth/2, matrixsSize);
ysRange = linspace(ysRangeCenter - ysRangeWidth/2, ysRangeCenter + ysRangeWidth/2, matrixsSize);
[Xs, Ys] = meshgrid(xsRange, ysRange);

xeRangeCenter = 0;
yeRangeCenter = xeRangeCenter;
xeRangeWidth = 0.08;
yeRangeWidth = xeRangeWidth;
matrixeSize = 120;
xeRange = linspace(xeRangeCenter - xeRangeWidth/2, xeRangeCenter + xeRangeWidth/2, matrixeSize);
yeRange = linspace(yeRangeCenter - yeRangeWidth/2, yeRangeCenter + yeRangeWidth/2, matrixeSize);
[Xe, Ye] = meshgrid(xeRange, yeRange);

l = 1;
wavelength = 0.003; % meters
huygenSource = exp(1j * l * atan2(Ys, Xs));
azcenter = 1;


phi = 0;

[axs, ays, afsz] = findActualCoords(phi, 0, azcenter, Xs, Ys);

RxEfield = hf(matrixeSize, axs, Xe, ays, Ye, afsz, huygenSource, wavelength);

f=figure(1);

set(gca, 'nextplot', 'replacechildren');  

clim manual;
a = 12*10^11;
clim([0 a])
xlim([-0.04 0.04]);
ylim([-0.04 0.04]);

v = VideoWriter("angmisganimationtake8.avi");
open(v)

for j=2:52
    i = j-2;
    l = i/2;
    phi = (i)*(pi/360);
    [axs, ays, afsz] = findActualCoords(phi, 0, azcenter, Xs, Ys);
    RxEfield = hf(matrixeSize, axs, Xe, ays, Ye, afsz, huygenSource, wavelength);
    %phase = angle(RxEfield);
    %imagesc(xeRange, yeRange, phase);
    intensity = abs(RxEfield).^2;
    imagesc(xeRange, yeRange, intensity);
    colormap(jet);
    colorbar;
    xlabel('X Position (m)'); % Label for x-axis
    ylabel('Y Position (m)'); % Label for y-axis
    title("Amplitude at phi=" + l + " degrees"); % Title for the plot
    writeVideo(v, getframe(f));
    writeVideo(v, getframe(f));
    writeVideo(v, getframe(f));
    if j==52
        for k = 1:20
            writeVideo(v, getframe(f));
            writeVideo(v, getframe(f));
            writeVideo(v, getframe(f));
        end
    end
end

% this separation is only here for holding at a certain frame

for j = 54 : 62
    i = j-2;
    l = i/2;
    phi = (i)*(pi/360);
    [axs, ays, afsz] = findActualCoords(phi, 0, azcenter, Xs, Ys);
    RxEfield = hf(matrixeSize, axs, Xe, ays, Ye, afsz, huygenSource, wavelength);
    %phase = angle(RxEfield);
    %imagesc(xeRange, yeRange, phase);
    intensity = abs(RxEfield).^2;
    imagesc(xeRange, yeRange, intensity);
    colormap(jet);
    colorbar;
    xlabel('X Position (m)'); % Label for x-axis
    ylabel('Y Position (m)'); % Label for y-axis
    title("Phase at phi=" + l + " degrees"); % Title for the plot
    writeVideo(v, getframe(f));
    writeVideo(v, getframe(f));
    writeVideo(v, getframe(f));
end


close(v);







%% helper function hf
function Ee = hf(matrixeSize, Xs, Xe, Ys, Ye, z, huygenSource, wavelength)
    Ee = zeros(matrixeSize, matrixeSize);
    k = 2 * pi/wavelength;
    for i = 1:matrixeSize
        for j = 1:matrixeSize 
            allRs = ((Xs-Xe(i, j)).^2 + (Ys-Ye(i, j)).^2 + z.^2).^0.5;
            contributions = huygenSource.*exp(1j*k*allRs)./allRs;
            hfSum = sum(contributions, "all");
            Ee(i, j) = hfSum*(1/(1j*wavelength));
        end
    end
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
%% helper function plotPhase
function plotPhase(efield, xrange, yrange, fieldname)
    phase = angle(efield);
    imagesc(xrange, yrange, phase);
    colormap(jet);
    colorbar;
    xlabel('X Position (m)'); % Label for x-axis
    ylabel('Y Position (m)'); % Label for y-axis
    title('Phase'); % Title for the plot
end
%% helper function findPurity (hope is correct)
function purity = findPurity(Efield, l, Xe, Ye)

    [Nr, Nphi] = size(Efield);
    r_max = sqrt(Xe(end)^2 + Ye(end)^2);
    r = linspace(0, r_max, Nr);
    phi = linspace(0, 2*pi, Nphi);
    [R, PHI] = meshgrid(r, phi);

    x_interp = R .* cos(PHI);
    y_interp = R .* sin(PHI);

    E_field_circ = interp2(Xe, Ye, Efield, x_interp, y_interp, 'spline').';
    
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

%% helper function generate axs, ays, az

function [axs, ays, afsz] = findActualCoords(phi, phi2, azcenter, Xs, Ys)
    fszc = azcenter*cos(phi)*cos(phi2); % distance from center of tx to plane of rx
    % Been Checked
    fsz = Xs * sin(phi) + Ys * cos(phi) * sin(phi2);
    % Been Checked
    afsz = -1*fsz + fszc;
    xOffset = azcenter*sin(phi)*cos(phi2); % distance from tx x center to rx x center
    % Been Checked
    axs = xOffset + Xs*cos(phi) + Ys*sin(phi)*sin(phi2);
    % Been Checked
    yOffset = azcenter*sin(phi2);
    % Been Checked
    ays = yOffset + Ys*cos(phi2);
end

%% helper function plotPurityVsAngle
function plotPurityVsAngle(anglefixed, value, min, max, numas, azcenter, l, Xs, Ys, Xe, Ye, matrixeSize, wavelength)
    phi = 0;
    phi2 = 0;
    if anglefixed=="phi2"
        phi2 = value;
    else 
        phi = value;
    end
    varyingas = linspace(min, max, numas);
    startingHFSs = exp(1j * l * atan2(Ys, Xs));
    puritys = zeros(numas);
    for indx = 1 : numas
        if anglefixed=="phi2"
            phi = varyingas(indx);
        else 
            phi2 = varyingas(indx);
        end
        [axs, ays, afsz] = findActualCoords(phi, phi2, azcenter, Xs, Ys);
        EndingEfield = hf(matrixeSize, axs, Xe, ays, Ye, afsz, startingHFSs, wavelength);
        Efield = steer(phi, EndingEfield, Xe, wavelength);
        puritys(indx) = findPurity(Efield, l, Xe, Ye);
    end
    plot(varyingas, puritys);
    if anglefixed=="phi2"
        xlabel("phi (radian) with phi2 = " + value);
    else 
        xlabel("phi2 (radian) with phi = " + value);
    end
    ylabel('OAM Purity');
    title('Purity vs angular misalignment');
    grid on;
end

%% helper function steer

function steeringplatesurfacehs = steer(phi, Efield, Xe, wavelength)
    steering1aphase = -sin(phi)*2*pi*Xe/wavelength;
    steeringplatesurfacehs = Efield.*exp(-1j*steering1aphase);
end
