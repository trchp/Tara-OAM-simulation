% SPP modeling with amplitude and wideband

%% defining Tx/SPP aperture

xsRangeCenter = 0;
ysRangeCenter = xsRangeCenter;
xsRangeWidth = 0.06;
ysRangeWidth = xsRangeWidth;
matrixsSize = 100;
xsRange = linspace(xsRangeCenter - xsRangeWidth/2, xsRangeCenter + xsRangeWidth/2, matrixsSize);
ysRange = linspace(ysRangeCenter - ysRangeWidth/2, ysRangeCenter + ysRangeWidth/2, matrixsSize);
[Xs, Ys] = meshgrid(xsRange, ysRange);

% OAM mode
l=2;

%% defining Rx aperture

xeRangeCenter = 0;
yeRangeCenter = xeRangeCenter;
xeRangeWidth = 0.1;
yeRangeWidth = xeRangeWidth;
matrixeSize = 130;
xeRange = linspace(xeRangeCenter - xeRangeWidth/2, xeRangeCenter + xeRangeWidth/2, matrixeSize);
yeRange = linspace(yeRangeCenter - yeRangeWidth/2, yeRangeCenter + yeRangeWidth/2, matrixeSize);
[Xe, Ye] = meshgrid(xeRange, yeRange);

% defining beam for which SPP is designed

wavelength = .002; % meters
k = 2*pi/wavelength;
totalz = 1; % propagation distance, meters
n=1.6; % refractive index of SPP material

%% defining SPP

stepheight = wavelength/(n-1);
height = mod((atan2(Ys,Xs)*l1), 2*pi) * stepheight / (2*pi);

dbloss = 0.41;
amplosspercm = 10^(dbloss/20);

%% empirically determined loss. 
% Based on desmos, these are very similar. Expect there are some matlab syntax errors with the first. 
% First is 1 / (empirical amp before / empirical amp after) ^ (height (AKA distance of propagation through material) / thickness of material used for measurements)
% Second is using amp after = amp before * exp (-1*attenuation coefficient * distance of propagation through material) to find attenuation coefficient empirically
% and implement that same equation


%amp = 1./(1.243^(height/0.002));
amp = 1*exp(-1*height1*108.848650032);

extfreq = 150*10^9; % 
extwl = 3*10^8/extfreq;

phase = 2*height*pi*(n-1)/extwl;
%phase = 2*height*pi*(n-1)/wavelength;

% determining the imposed phase, with an extreme frequency (literally anything) or intended wavelength

huygenSource = amp.*exp(1*j*phase);

%huygenSource = awgn(huygenSource, 1/1000); % adding noise if desired. not found to have much effect
    
RxEfield = hf(matrixeSize, matrixeSize, Xs, Xs, Ys, Ys, totalz, huygenSource, wavelength);

%demodEfield = RxEfield.*exp(-1*j*l*atan2(Ye, Xe)); % demodulation. Could also be modeled with a SPP

%plotIntensity(RxEfield, xeRange, yeRange, "RxEfield");
%plotPhase(huygenSource, xsRange, ysRange, "huygenSource");

%purity = findPurity(demodEfield, 0, Xe, Ye, 0)



function plotIntensity(efield, xrange, yrange, fieldname)
    intensity = abs(efield).^2;
    imagesc(xrange, yrange, intensity);
    colormap(jet);
    colorbar;
    xlabel('X Position (m)'); % Label for x-axis
    ylabel('Y Position (m)'); % Label for y-axis
    ax = gca;
    ax.YDir = 'normal'; % Set YDir to 'normal'
    title('Intensity of ' + fieldname); % Title for the plot
end

function plotPhase(efield, xrange, yrange, fieldname)
    phase = angle(efield);
    imagesc(xrange, yrange, phase);
    colormap(jet);
    colorbar;
    xlabel('X Position (m)'); % Label for x-axis
    ylabel('Y Position (m)'); % Label for y-axis
    ax = gca; % Get handle to current axes
    ax.YDir = 'normal'; % Set YDir to 'normal'
    title('Phase of ' + fieldname); % Title for the plot
end

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
