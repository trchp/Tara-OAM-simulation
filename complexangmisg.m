
%% Define transmit aperture dimensions

xsRangeCenter = 0;
ysRangeCenter = xsRangeCenter;
xsRangeWidth = 0.06;
ysRangeWidth = xsRangeWidth;
matrixsSize = 80;
xsRange = linspace(xsRangeCenter - xsRangeWidth/2, xsRangeCenter + xsRangeWidth/2, matrixsSize);
ysRange = linspace(ysRangeCenter - ysRangeWidth/2, ysRangeCenter + ysRangeWidth/2, matrixsSize);
[Xs, Ys] = meshgrid(xsRange, ysRange);


%% Define transmit aperture position relative to rx aperture


phi = pi/6; % affects x
phi2 = 0; % affects y (and also x)
azcenter = 0.4; % meters (actual propagation distance for ray transmitted 
% from center of tx

[axs, ays, afsz] = findActualCoords(phi, phi2, azcenter, Xs, Ys);

%fszc = azcenter*cos(phi)*cos(phi2); % distance from center of tx to plane of rx
%fsz = Xs * sin(phi) + Ys * cos(phi) * sin(phi2);
%afsz = -1*fsz + fszc;  % distance in the z direction from any tx point to the plane of rx
%xOffset = azcenter*sin(phi)*cos(phi2); % distance from tx x center to rx x center in the x direction
%axs = xOffset + Xs*cos(phi) + Ys*sin(phi)*sin(phi2); % distance from any tx point to rx x center in the x direction
%yOffset = azcenter*sin(phi2);
%ays = yOffset + Ys*cos(phi2); % actual y 

l = 1;
wavelength = 0.003; % meters

huygenSource = exp(1j * l * atan2(Ys, Xs));

%% Define Rx aperture

xeRangeCenter = 0;
yeRangeCenter = xeRangeCenter;
xeRangeWidth = 0.08;
%xeRangeWidth = 0.3;
yeRangeWidth = xeRangeWidth;
matrixeSize = 120;
%matrixeSize = 300;
xeRange = linspace(xeRangeCenter - xeRangeWidth/2, xeRangeCenter + xeRangeWidth/2, matrixeSize);
yeRange = linspace(yeRangeCenter - yeRangeWidth/2, yeRangeCenter + yeRangeWidth/2, matrixeSize);
[Xe, Ye] = meshgrid(xeRange, yeRange);

%%
RxEfield = hf(matrixeSize, axs, Xe, ays, Ye, afsz, huygenSource, wavelength);

%% attempts at steering

%s = atan2(cos(phi)*sin(phi2), 1);
%n = atan2(cos(phi2)/cos(s), 1);

%steering1aphase = -sin(phi)*2*pi*Xe/wavelength
%steeringphase = -2*pi/wavelength*sin(s)*(Xe*cos(n) + Ye*-1*sin(n));
%steerfieldHsource = RxEfield.*exp(-1j*steeringphase);
%RxEfield2 = hf(matrixeSize, Xe, Xe, Ye, Ye, 0.2, steerfieldHsource, wavelength);

plotPhase(RxEfield, xeRange, yeRange, "RxEfield");

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

%% running plots
%purity = findPurity(RxEfield, 1, xeRange, yeRange)

%plotPhase(RxEfield2, xeRange, yeRange, "RxEfield");
%plotPurityVsAngle("phi", pi/6, -2, 2, 33, 1, 1, Xs, Ys, Xe, Ye, matrixeSize, wavelength);

%plotPurityVsAngle("phi2", 0, -pi/6, pi/6, 9, 1, 1, Xs, Ys, Xe, Ye, matrixeSize, wavelength)

phi = 0;
    phi2 = 0;
           min = -pi/6;
           max = pi/6;
           numas = 9;
    
    varyingas = linspace(min, max, numas);
    startingHFSs = exp(1j * l * atan2(Ys, Xs));
    puritys = zeros(1, numas);
    for indx = 1 : numas
        
            phi = varyingas(indx);
        
        [axs, ays, afsz] = findActualCoords(phi, phi2, azcenter, Xs, Ys);
        EndingEfield = hf(matrixeSize, axs, Xe, ays, Ye, afsz, startingHFSs, wavelength);
        Efield = steer(phi, EndingEfield, Xe, wavelength);
        %Efield = EndingEfield;
        puritys(indx) = findPurity(Efield, l, Xe, Ye);
    end
    plot(varyingas, puritys, "LineWidth", 2);
    ylim([0 0.8]);
   
        xlabel("phi (radian) with phi2 = 0");
    
    ylabel('OAM Purity');
    title('Purity vs angular misalignment');
    grid on;
