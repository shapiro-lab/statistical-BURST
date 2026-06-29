%% Import image dataset 
% This part is importing image dataset. Feel free to adjust this section to
% your own data

% Is your data intensity image (nonnegative real) or IQ (complex)?
is_IQ = true; % please enter true or false

basedir = 'enter your folder/burst';

if ~is_IQ
    file_start = 'image_block_*'; % nonnegative real-value image files
else 
    file_start = 'IQ_block_*'; % IQ files
end

fileList = dir(fullfile(basedir, file_start));
for i = 1:numel(fileList)
    
    fileName = fileList(i).name;
    filePath = fullfile(basedir, fileName);

    data = load(filePath); 
    if ~is_IQ
        image_seq(:,:,i) = data.RData;
    else
        image_seq(:,:,i) = complex(data.IQ(:,:,1), data.IQ(:,:,2)); % assuming that real and complex parts are stacked
    end

    % extract the voltage value from the folder name
    alpha = regexp(fileName, '([\d.]+).mat', 'tokens');
    a(i) = str2double(alpha{1}{1});

end

image_seq(:,:,a) = image_seq; % sorting based on the index
%% visualization of image sequence
% use this for checking collapse frame and when pure background starts
figure;
for i = 1:numel(fileList)

    imagesc(squeeze(abs(image_seq(:,:,i)))); colorbar; clim([0 max(abs(image_seq(:,:,3)), [], 'all')/5])
        colormap hot; colorbar
    title(num2str(i)); pause;

end

%% Parameter settings

start_high = 3; % start of high voltage 
num_collapse = 1; % number of frames in which GV collapses
start_backgr = 4; % frame number since which you are sure that there is no GV signal at all
end_backgr = 102; % end of the background frames
p_value = 1e-4; % desired p-value
IQ_processing = true; % when data is IQ and if you want IQ preprocessing (tissue background removal)

p_value = p_value / num_collapse; % Bonferroni correction

if is_IQ && ~IQ_processing
    
    image_seq = image_seq / 1e4; % IQ data sometimes have very large number, so it is more stable to compress
    Z = image_seq; % this is for mahalanobis later
    image_seq = abs(image_seq); % directly acquire real-value image

elseif is_IQ && IQ_processing

    image_seq = image_seq / 1e4; % IQ data sometimes have very large number, so it is more stable to compress
    Z = image_seq; % this is for mahalanobis later
    image_seq = abs(Z - mean(Z(:,:,start_backgr:end_backgr),3)); % remove background coponents in complex plane

end

%% subtraction-based BURST processing (conventional signal unmixing)
pAM2D_sub = image_seq - mean(image_seq(:,:,start_backgr:end_backgr),3);
pAM2D_sub = pAM2D_sub.*(pAM2D_sub>0);
pAM2D_fin = sum(pAM2D_sub(:,:,start_high + (0:num_collapse-1)), 3);

figure; imagesc(pAM2D_fin); colormap hot; axis image; title('subtraction-based'); colorbar; imcontrast

%% Correlation-based BURST processing (or tCNR-based)

for i = 1:num_collapse
    pAM2D_z = zscore(image_seq(:,:,[start_high+i-1 start_backgr:end_backgr]), [], 3);
    burst(:,:,i) = (pAM2D_z(:,:, 1) - mean(pAM2D_z(:,:,2:end), 3))/sqrt(size(pAM2D_z, 3));
    t_map(:,:,i) = burst(:,:,i).*sqrt((size(pAM2D_z, 3)-2)./(1-burst(:,:,i).^2));
end

p = 1 - p_value;
t = tinv(p, size(pAM2D_z, 3)-2);
threshold = t/sqrt(size(pAM2D_z, 3)-2+t^2);

burst_unmix_final = sum((image_seq(:,:,start_high + (0:num_collapse-1)) - mean(image_seq(:,:,start_backgr:end_backgr),3)).*(burst>threshold), 3);

figure; subplot(1,2,1); imagesc(t_map(:,:,1)); axis image; title('t_map 1st'); colorbar;
subplot(1,2,2); imagesc(burst_unmix_final); colormap hot; axis image; title('tCNR-based'); colorbar; imcontrast
%% GUI for changing p-value
% burstThresholdGUI(image_seq(:,:,start_high + (0:num_collapse-1)) - mean(image_seq(:,:,start_backgr:end_backgr),3),...
%     1 - tcdf(t_map, end_backgr-start_backgr));

function burstThresholdGUI(beta, p_map)

    [H, W, num_collapse] = size(beta);

    % Create figure
    hFig = figure('Name','Burst Threshold Explorer','NumberTitle','off', ...
        'MenuBar','none','ToolBar','none','Resize','on', ...
        'Units','normalized','Position',[0.2 0.1 0.3 0.5]);

    hAx = axes('Parent', hFig, 'Units','pixels', 'Position',[50 80 512 512]);
    colormap(hAx, 'hot');

    % Slider for p-value control
    hSlider = uicontrol('Style','slider', ...
        'Min',0,'Max',10,'Value',4, ...
        'Units','pixels','Position',[50 30 400 20], ...
        'Callback',@updateImage);

    % Label
    uicontrol('Style','text', ...
        'Position',[460 30 100 20], ...
        'String','p-value', ...
        'HorizontalAlignment','left');

    % Edit box
    hEdit = uicontrol('Style','edit', ...
        'String',num2str(hSlider.Value,'%.4f'), ...
        'Position',[460 55 80 20], ...
        'Callback',@editCallback);

    % Initial draw
    updateImage();

    function updateImage(~,~)
        p = 10^(-hSlider.Value);
        set(hEdit, 'String', sprintf('%.3f', log10(p)));

        % Use precomputed beta and stat_map
        mask = p_map < p/num_collapse;
        burst_unmix = sum(beta .* mask, 3);

        % Flatten and remove NaNs
        vals = beta.*(p_map < 1e-4/num_collapse);
        vals = vals(:);
        vals = vals(~isnan(vals));

        % Percentiles
        low  = prctile(vals, 0.1);
        high = prctile(vals, 99.9);

        % Display
        imagesc(burst_unmix, 'Parent', hAx);
        axis(hAx, 'image', 'on');
        title(hAx, sprintf('p = %.6f', p));
        colorbar('peer', hAx);
        hAx.CLim = [low high];
        drawnow;
    end

    function editCallback(src,~)
        val = str2double(src.String);
        if isnan(val) || val < hSlider.Min || val > hSlider.Max
            src.String = sprintf('%.4f', hSlider.Value);
            return;
        end
        hSlider.Value = val;
        updateImage();
    end
end

%% Nakagami-based BURST processing
m = nakagami_m_image(image_seq(:,:,start_backgr:end_backgr),1000);
m3 = repmat(m, 1, 1, max(1, num_collapse));

omega = mean(image_seq(:,:,start_backgr:end_backgr).^2, 3);
R = image_seq(:,:,start_high + (0:num_collapse-1)).^2./omega;
p_nakagami = 1-fcdf(R, 2*m3, 2*(end_backgr-start_backgr+1+num_collapse)*m3);
burst_nakagami = sum((image_seq(:,:,start_high + (0:num_collapse-1)) - mean(image_seq(:,:,start_backgr:end_backgr), 3)).*(p_nakagami < p_value), 3);

% for the intensity map, instead of F_collapse - average(F_post), you can
% also use F_collapse - mode value based on estimated Nakagami distribution
% background_estimate = sqrt(omega.*(m - 0.5)./m);
% burst_nakagami = sum((image_seq(:,:,start_high + (0:num_collapse-1)) - background_estimate).*(p_nakagami < p_value), 3);

figure; subplot(1,2,1); imagesc(R(:,:,1)); axis image; title('R_{BURST} 1st frame'); colorbar; imcontrast
subplot(1,2,2); imagesc(burst_nakagami); colormap hot; axis image; title('Nakagami-based'); colorbar; imcontrast

function m = nakagami_m_image(Y, nIter, eps0)
% per-pixel Nakagami m estimation
%
% Model:
%   X ~ Nakagami(m, Omega)
%   Y(:,:,t) = Xt            background
%
% Inputs
%   Y     : HxWxT array (amplitudes, nonnegative recommended)
%   nIter : Newton iterations for m (default 8)
%   eps0  : small floor to avoid log(0) (default 1e-12)
%
% Outputs (all HxW)
%   m     : HxW map of Nakagami shape parameter

    if nargin < 2 || isempty(nIter), nIter = 8; end
    if nargin < 3 || isempty(eps0),  eps0  = 1e-12; end
    
    Y = double(Y);  % psi/gammainc are happiest in double
    [H,W,T] = size(Y);
    assert(T >= 2, 'Need T>=3 to fit Nakagami baseline per pixel robustly.');

    % ---- Accumulate sufficient statistics over t=2..T (one loop over frames)
    s2   = zeros(H,W);   % sum x^2
    s4   = zeros(H,W);   % sum x^4  (for Var of x^2)
    slog = zeros(H,W);   % sum log(x)
    
    for t = 1:T
        x = max(Y(:,:,t), eps0);   % floor to avoid log(0)
        x2 = x.^2;
        s2   = s2   + x2;
        s4   = s4   + x2.^2;       % = x^4
        slog = slog + log(x);
    end

    Omega = s2 / T;                % E[X^2] = Omega
    Ez2   = s4 / T;                % E[(X^2)^2]
    Varz  = max(Ez2 - Omega.^2, 0);% Var(Z) where Z=X^2
    
    % ---- Initial m via method-of-moments on Z=X^2 ~ Gamma(m, scale=Omega/m)
    m = Omega.^2 ./ max(Varz, eps0);
    m = max(m, 1e-3);
    
    % ---- MLE refine m via Newton: 
    % log(m) - psi(m) = log(mean(Z)) - mean(log(Z))
    % mean(log(Z)) = mean(log(X^2)) = mean(2*log(X)) = 2*mean(log(X))
    mean_logZ = 2*(slog/T);
    rhs = log(max(Omega, eps0)) - mean_logZ;

    for k = 1:nIter
        f  = log(m) - psi(m) - rhs;
        fp = 1./m - psi(1,m);      % trigamma is psi(1,·)
        step = f ./ fp;
        m = max(m - step, 1e-6);
        % optional early stop:
        if max(abs(step(:))) < 1e-6, break; end
    end
end

%% Rayleigh-based BURST processing (recommend when IQ images are available)
bg_sum = sum(image_seq(:,:,start_backgr:end_backgr).^2, 3);
high_sq = image_seq(:,:,start_high + (0:num_collapse-1)).^2;

p_rayleigh = (bg_sum./(bg_sum + high_sq)).^(end_backgr-start_backgr+1);
burst_rayleigh = sum((image_seq(:,:,start_high + (0:num_collapse-1)) - mean(image_seq(:,:,start_backgr:end_backgr), 3)).*(p_rayleigh < p_value), 3);

% for the intensity map, instead of F_collapse - average(F_post), you can
% also use F_collapse - mode value based on estimated rayleigh distribution
% background_estimate = sqrt(bg_sum/2/(end_backgr-start_backgr+1)); % Note: this is a biased estiator of rayleigh parameter
% burst_rayleigh = sum((image_seq(:,:,start_high + (0:num_collapse-1)) - background_estimate).*(p_rayleigh < p_value), 3);

figure; subplot(1,2,1); imagesc(R(:,:,1)); axis image; title('R_{BURST} 1st frame'); colorbar; imcontrast;
subplot(1,2,2); imagesc(burst_rayleigh); colormap hot; axis image; title('Rayleigh-based'); colorbar; imcontrast
%% Mahalanobis distance (usable only when IQ iamges are available)
% Assume that image_seq is zero-centered based on post-collapse IQ frames

if is_IQ

    [nx, ny, nt] = size(Z);
    Zref = Z(:,:,start_backgr:end_backgr); % background frames

    % pixel-wise Mahalanobis distance calculation
    mur = mean(real(Zref), 3);
    mui = mean(imag(Zref), 3);
    
    xrc_ref = real(Zref) - mur;
    xic_ref = imag(Zref) - mui;
    
    xrc_all = real(Z) - mur;
    xic_all = imag(Z) - mui;
    
    den = end_backgr-start_backgr;
    a = sum(xrc_ref.^2, 3) / den;
    b = sum(xrc_ref .* xic_ref, 3) / den;
    c = sum(xic_ref.^2, 3) / den;
    
    eps_reg = 1e-8;
    a = a + eps_reg;
    c = c + eps_reg;
    
    detS = max(a .* c - b.^2, eps_reg);
    
    mahala = sqrt(max((c .* xrc_all.^2 - 2*b.*xrc_all.*xic_all + a.*xic_all.^2) ./ detS, 0)); % Mahalanobis distance

    m = den+1;
    threshold = sqrt(2*(m+1)*(m-1)/m/(m-2)*finv(1 - p_value, 2, m-2));
    burst_mahala = sum((image_seq(:,:,start_high + (0:num_collapse-1)) - mean(image_seq(:,:,start_backgr:end_backgr), 3))...
        .*(mahala(:,:,start_high + (0:num_collapse-1))>threshold),3);

    figure; subplot(1,2,1); imagesc(mahala(:,:,start_high)); axis image; title('Mahalanobis 1st frame'); colorbar; imcontrast;
    subplot(1,2,2); imagesc(burst_mahala); colormap hot; axis image; title('Mahalanobis-based'); colorbar; clim([0 max(burst_mahala(:))]); imcontrast
end
