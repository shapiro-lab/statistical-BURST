%% Import image dataset 
% This part is importing image dataset. Feel free to adjust this section to
% your own data

% Is your data intensity image (nonnegative real) or IQ (complex)?
is_IQ = false; % please enter true or false

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

image_seq(:,:,a) = image_seq;

burstModeGUI(image_seq, is_IQ, 3, 1, 4, 102)

function burstModeGUI(image_seq_raw, is_IQ, start_high, num_collapse, start_backgr, end_backgr, varargin)
% burstModeGUI
% Interactive BURST explorer based on the supplied processing code.
%
% Example:
%   burstModeGUI(image_seq, is_IQ, start_high, num_collapse, start_backgr, end_backgr, ...
%       'IQProcessing', IQ_processing, 'InitialPValue', 1e-4);
%
% Inputs
%   image_seq_raw  : HxWxT raw image sequence. Real-valued intensity images or
%                    complex IQ data, matching your original script.
%   is_IQ          : true for IQ data, false for intensity data.
%   start_high     : first collapse frame index.
%   num_collapse   : number of collapse frames.
%   start_backgr   : first pure-background frame index.
%   end_backgr     : last pure-background frame index.
%
% Name-value pairs
%   'IQProcessing' : true/false, default true.
%   'InitialPValue': default 1e-4.
%
% Notes on displayed left-panel maps (current assumption):
%   - subtraction (conventional): no statistical map, only pAM2D_fin shown.
%   - correlation or tCNR-based: t_map
%   - Nakagami: R = F_collapse.^2 / omega
%   - Rayleigh: -log10(p_rayleigh)
%   - Mahalanobis: Mahalanobis distance map
%
% Update those mappings if you prefer a different diagnostic image.

    p = inputParser;
    addParameter(p, 'IQProcessing', true, @(x)islogical(x) || isnumeric(x));
    addParameter(p, 'InitialPValue', 1e-4, @(x)isnumeric(x) && isscalar(x) && x > 0 && x < 1);
    parse(p, varargin{:});

    IQ_processing = logical(p.Results.IQProcessing);
    initial_p = p.Results.InitialPValue;

    % ---------------------------------------------------------------------
    % Preprocessing (matches the supplied script as closely as possible)
    % ---------------------------------------------------------------------
    image_seq_raw = double(image_seq_raw);
    Z = [];

    if ~is_IQ
        image_seq = image_seq_raw;
    else
        Z = image_seq_raw / 1e4;
        if IQ_processing
            image_seq = abs(Z - mean(Z(:,:,start_backgr:end_backgr), 3));
        else
            image_seq = abs(Z);
        end
    end

    beta = image_seq(:,:,start_high + (0:num_collapse-1)) - mean(image_seq(:,:,start_backgr:end_backgr), 3);

    results = computeAllModes(image_seq, Z, is_IQ, beta, start_high, num_collapse, start_backgr, end_backgr);
    preferredOrder = {'subtraction_conventional', 'correlation_tcNR', 'nakagami', 'rayleigh', 'mahalanobis'};
    modeNames = preferredOrder(isfield(results, preferredOrder));
    defaultMode = 'correlation_tcNR';
    if ~isfield(results, defaultMode)
        defaultMode = modeNames{1};
    end

    % ---------------------------------------------------------------------
    % GUI
    % ---------------------------------------------------------------------
    hFig = figure('Name', 'BURST Mode Explorer', ...
                  'NumberTitle', 'off', ...
                  'MenuBar', 'none', ...
                  'ToolBar', 'none', ...
                  'Units', 'normalized', ...
                  'Position', [0.08 0.08 0.82 0.80], ...
                  'Color', 'w');

    hPanelCtrl = uipanel('Parent', hFig, 'Units', 'normalized', ...
                         'Position', [0.02 0.02 0.96 0.16], 'Title', 'Controls');

    hPanelView = uipanel('Parent', hFig, 'Units', 'normalized', ...
                         'Position', [0.02 0.20 0.96 0.78], 'Title', 'Display');

    hAxLeft  = axes('Parent', hPanelView, 'Units', 'normalized', 'Position', [0.05 0.10 0.40 0.83]);
    hAxRight = axes('Parent', hPanelView, 'Units', 'normalized', 'Position', [0.55 0.10 0.40 0.83]);
    colormap(hAxLeft,  'hot');
    colormap(hAxRight, 'hot');
    % Mode selector
    uicontrol(hPanelCtrl, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.02 0.62 0.16 0.22], 'String', 'Mode', ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

    displayNames = cellfun(@(f)results.(f).displayName, modeNames, 'UniformOutput', false);
    defaultIdx = find(strcmp(modeNames, defaultMode), 1);
    hMode = uicontrol(hPanelCtrl, 'Style', 'popupmenu', 'Units', 'normalized', ...
        'Position', [0.02 0.34 0.28 0.30], 'String', displayNames, 'Value', defaultIdx, ...
        'Callback', @updateGUI);

    % Frame list
    uicontrol(hPanelCtrl, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.34 0.62 0.16 0.22], 'String', 'Collapsed frame', ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

    frameItems = arrayfun(@(k)frameLabel(k, start_high), 1:max(num_collapse, 1), 'UniformOutput', false);
    hFrame = uicontrol(hPanelCtrl, 'Style', 'listbox', 'Units', 'normalized', ...
        'Position', [0.34 0.22 0.24 0.38], 'String', frameItems, 'Value', 1, ...
        'Min', 0, 'Max', 1, 'Callback', @updateGUI);

    uicontrol(hPanelCtrl, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.34 0.10 0.14 0.10], 'String', 'Left contrast', ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
    hLeftContrastBtn = uicontrol(hPanelCtrl, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.34 0.02 0.13 0.10], 'String', 'Open imcontrast', ...
        'Callback', @openLeftContrastTool);
    hAutoContrastLeft = uicontrol(hPanelCtrl, 'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [0.48 0.02 0.10 0.10], 'String', 'Auto', 'Value', 1, ...
        'BackgroundColor', 'w', 'Callback', @updateGUI);

    % p-value slider (log10 scale: p = 10^(-value))
    uicontrol(hPanelCtrl, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.60 0.62 0.12 0.18], 'String', 'p-value', ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

    hP = uicontrol(hPanelCtrl, 'Style', 'slider', 'Units', 'normalized', ...
        'Position', [0.60 0.40 0.22 0.16], 'Min', 0, 'Max', 10, ...
        'Value', min(max(-log10(max(initial_p, realmin)), 0), 10), ...
        'SliderStep', [0.01 0.10], 'Callback', @updateGUI);

    hPEdit = uicontrol(hPanelCtrl, 'Style', 'edit', 'Units', 'normalized', ...
        'Position', [0.84 0.40 0.10 0.18], 'String', sprintf('%.4g', initial_p), ...
        'Callback', @pEditCallback);

    uicontrol(hPanelCtrl, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.62 0.20 0.14 0.12], 'String', 'Right contrast', ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

    hContrastBtn = uicontrol(hPanelCtrl, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.62 0.08 0.16 0.14], 'String', 'Open imcontrast', ...
        'Callback', @openContrastTool);

    hAutoContrast = uicontrol(hPanelCtrl, 'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [0.80 0.09 0.16 0.12], 'String', 'Auto contrast', 'Value', 1, ...
        'BackgroundColor', 'w', 'Callback', @updateGUI);

    hContrastInfo = uicontrol(hPanelCtrl, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.62 0.00 0.34 0.08], ...
        'String', 'Turn off auto to keep the same contrast range across modes', ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'left');

    setappdata(hFig, 'lastModeIdx', defaultIdx);
    updateGUI();

    % ---------------------------------------------------------------------
    % Nested callbacks
    % ---------------------------------------------------------------------
    function updateGUI(~, ~)
        prevLeftCLim = get(hAxLeft, 'CLim');
        prevRightCLim = get(hAxRight, 'CLim');
        modeIdx = get(hMode, 'Value');
        lastModeIdx = getappdata(hFig, 'lastModeIdx');
        if isempty(lastModeIdx)
            lastModeIdx = modeIdx;
        end
        if modeIdx ~= lastModeIdx
            set(hAutoContrastLeft, 'Value', 1);
            set(hAutoContrast, 'Value', 1);
            setappdata(hFig, 'lastModeIdx', modeIdx);
        end
        modeKey = modeNames{modeIdx};
        S = results.(modeKey);

                        if S.showOnlyFinal || isempty(S.leftMap)
            nFramesAvail = max(1, S.numFrames);
        else
            nFramesAvail = max(1, size(S.leftMap, 3));
        end
        frameVal = get(hFrame, 'Value');
        if isempty(frameVal)
            frameVal = 1;
        end
        frameIdx = clamp(round(frameVal(1)), 1, nFramesAvail);
                frameItems = arrayfun(@(k)frameLabel(k, start_high), 1:nFramesAvail, 'UniformOutput', false);
        set(hFrame, 'String', frameItems, 'Value', frameIdx);

        pUser = 10^(-get(hP, 'Value'));
        set(hPEdit, 'String', sprintf('%.4g', pUser));

        % LEFT PANEL
        cla(hAxLeft);
        if S.showOnlyFinal
            colorbar(hAxLeft, 'off');
            axis(hAxLeft, 'off');
            title(hAxLeft, '');
        else
            axis(hAxLeft, 'on');
            leftImg = S.leftMap(:,:,frameIdx);
            imagesc(leftImg, 'Parent', hAxLeft);
            axis(hAxLeft, 'image');
            colorbar(hAxLeft);
            title(hAxLeft, sprintf('%s | frame %d', S.leftLabel, start_high + frameIdx - 1));
            applyLeftContrast(prevLeftCLim, robustCLim(leftImg));
        end

        % RIGHT PANEL
        cla(hAxRight);
        if S.showOnlyFinal
            imagesc(S.finalImage, 'Parent', hAxRight);
            axis(hAxRight, 'image');
            colorbar(hAxRight);
            title(hAxRight, 'Final BURST image (subtraction)');
            applyRightContrast(prevRightCLim, zeroLowCLim(robustCLim(S.finalImage)));
        else
                        nFramesRight = min(size(S.pMap, 3), size(S.betaForDisplay, 3));
            pMapUse = S.pMap(:,:,1:nFramesRight);
            betaUse = S.betaForDisplay(:,:,1:nFramesRight);
            pThr = pUser / max(1, nFramesRight);
                        mask = pMapUse < pThr;
                        rightImg = sum(betaUse .* mask, 3);

                        vals = betaUse .* (pMapUse < 1e-4 / max(1, nFramesRight));
            vals = vals(:);
            vals = vals(~isnan(vals));
            if isempty(vals)
                cLim = robustCLim(rightImg);
            else
                low = prctile(vals, 0.1);
                high = prctile(vals, 99.9);
                cLim = [low, high];
                if ~all(isfinite(cLim)) || cLim(1) >= cLim(2)
                    cLim = robustCLim(rightImg);
                end
            end

            imagesc(rightImg, 'Parent', hAxRight);
            axis(hAxRight, 'image');
            colorbar(hAxRight);
            title(hAxRight, sprintf('%s | p = %.3g', S.rightLabel, pUser));
            applyRightContrast(prevRightCLim, zeroLowCLim(cLim));
        end

        drawnow;
    end

    function pEditCallback(src, ~)
        pVal = str2double(get(src, 'String'));
        if isnan(pVal) || pVal <= 0 || pVal > 1
            set(src, 'String', sprintf('%.4g', 10^(-get(hP, 'Value'))));
            return;
        end
        set(hP, 'Value', min(max(-log10(max(pVal, realmin)), get(hP, 'Min')), get(hP, 'Max')));
        updateGUI();
    end

    function openLeftContrastTool(~, ~)
        if ishghandle(hAxLeft)
            set(hAutoContrastLeft, 'Value', 0);
            set(hAxLeft, 'CLimMode', 'manual');
            imcontrast(hAxLeft);
        end
    end

    function applyLeftContrast(prevCLim, autoCLim)
        if get(hAutoContrastLeft, 'Value')
            set(hAxLeft, 'CLimMode', 'manual');
            caxis(hAxLeft, autoCLim);
        else
            if isempty(prevCLim) || any(~isfinite(prevCLim)) || prevCLim(1) >= prevCLim(2)
                set(hAxLeft, 'CLimMode', 'manual');
                caxis(hAxLeft, autoCLim);
            else
                set(hAxLeft, 'CLimMode', 'manual');
                caxis(hAxLeft, prevCLim);
            end
        end
    end

    function openContrastTool(~, ~)
        if ishghandle(hAxRight)
            set(hAutoContrast, 'Value', 0);
            set(hAxRight, 'CLimMode', 'manual');
            imcontrast(hAxRight);
        end
    end

    function applyRightContrast(prevCLim, autoCLim)
        if get(hAutoContrast, 'Value')
            set(hAxRight, 'CLimMode', 'manual');
            caxis(hAxRight, autoCLim);
        else
            if isempty(prevCLim) || any(~isfinite(prevCLim)) || prevCLim(1) >= prevCLim(2)
                set(hAxRight, 'CLimMode', 'manual');
                caxis(hAxRight, autoCLim);
            else
                set(hAxRight, 'CLimMode', 'manual');
                caxis(hAxRight, prevCLim);
            end
        end
    end
end

function cLim = zeroLowCLim(cLimIn)
    cLim = cLimIn;
    if ~isfinite(cLim(2)) || cLim(2) <= 0
        cLim = cLimIn;
        return;
    end
    cLim = [0, cLim(2)];
end

% ======================================================================
function results = computeAllModes(image_seq, Z, is_IQ, beta, start_high, num_collapse, start_backgr, end_backgr)

    nBg = end_backgr - start_backgr + 1;
    bgMean = mean(image_seq(:,:,start_backgr:end_backgr), 3);

    % ------------------------------------------------------------------
    % subtraction (conventional)
    % ------------------------------------------------------------------
    pAM2D_sub = image_seq - bgMean;
    pAM2D_sub = pAM2D_sub .* (pAM2D_sub > 0);
    pAM2D_fin = sum(pAM2D_sub(:,:,start_high + (0:num_collapse-1)), 3);

    results.subtraction_conventional = struct( ...
        'displayName', 'subtraction (conventional)', ...
        'showOnlyFinal', true, ...
        'leftMap', [], ...
        'leftLabel', '', ...
        'pMap', [], ...
        'betaForDisplay', [], ...
        'finalImage', pAM2D_fin, ...
        'rightLabel', 'Final BURST image', ...
        'numFrames', max(1, num_collapse));

    % ------------------------------------------------------------------
    % correlation or tCNR-based
    % ------------------------------------------------------------------
    burst = zeros(size(beta));
    t_map = zeros(size(beta));
    for i = 1:num_collapse
        pAM2D_z = zscore(image_seq(:,:,[start_high+i-1, start_backgr:end_backgr]), [], 3);
        burst(:,:,i) = (pAM2D_z(:,:,1) - mean(pAM2D_z(:,:,2:end), 3)) / sqrt(size(pAM2D_z, 3));
        t_map(:,:,i) = burst(:,:,i) .* sqrt((size(pAM2D_z, 3)-2) ./(1 - burst(:,:,i).^2));
    end
    df_t = size(pAM2D_z, 3) - 2;
    p_t = max(1 - tcdf(t_map, df_t), realmin);

    results.correlation_tcNR = struct( ...
        'displayName', 'correlation or tCNR-based (default)', ...
        'showOnlyFinal', false, ...
        'leftMap', t_map, ...
        'leftLabel', 't-map', ...
        'pMap', p_t, ...
        'betaForDisplay', beta, ...
        'finalImage', sum(beta .* (p_t < (1e-4 / num_collapse)), 3), ...
        'rightLabel', 'Final BURST image', ...
        'numFrames', max(1, num_collapse));

    % ------------------------------------------------------------------
    % Nakagami
    % ------------------------------------------------------------------
    m = nakagami_m_image(image_seq(:,:,start_backgr:end_backgr), 1000);
    omega = mean(image_seq(:,:,start_backgr:end_backgr).^2, 3);
    R = image_seq(:,:,start_high + (0:num_collapse-1)).^2 ./ max(omega, eps);
    m3 = repmat(m, 1, 1, max(1, num_collapse));
    p_nakagami = max(1 - fcdf(R, 2*m3, 2*(nBg + num_collapse)*m3), realmin);

    results.nakagami = struct( ...
        'displayName', 'Nakagami', ...
        'showOnlyFinal', false, ...
        'leftMap', R, ...
        'leftLabel', 'Nakagami statistic R', ...
        'pMap', p_nakagami, ...
        'betaForDisplay', beta, ...
        'finalImage', sum(beta .* (p_nakagami < (1e-4 / num_collapse)), 3), ...
        'rightLabel', 'Final BURST image', ...
        'numFrames', max(1, num_collapse));

    % ------------------------------------------------------------------
    % Rayleigh
    % ------------------------------------------------------------------
    bg_sum = sum(image_seq(:,:,start_backgr:end_backgr).^2, 3);
    high_sq = image_seq(:,:,start_high + (0:num_collapse-1)).^2;
    p_rayleigh = max((bg_sum ./ max(bg_sum + high_sq, eps)).^nBg, realmin);
    R_rayleigh = image_seq(:,:,start_high + (0:num_collapse-1)).^2 ./ max(mean(image_seq(:,:,start_backgr:end_backgr).^2, 3), eps);

    results.rayleigh = struct( ...
        'displayName', 'Rayleigh', ...
        'showOnlyFinal', false, ...
        'leftMap', R_rayleigh, ...
        'leftLabel', 'Rayleigh statistic R', ...
        'pMap', p_rayleigh, ...
        'betaForDisplay', beta, ...
        'finalImage', sum(beta .* (p_rayleigh < (1e-4 / num_collapse)), 3), ...
        'rightLabel', 'Final BURST image', ...
        'numFrames', max(1, num_collapse));

    % ------------------------------------------------------------------
    % Mahalanobis (IQ only)
    % ------------------------------------------------------------------
    if is_IQ && ~isempty(Z)
        Zref = Z(:,:,start_backgr:end_backgr);
        mur = mean(real(Zref), 3);
        mui = mean(imag(Zref), 3);

        xrc_ref = real(Zref) - mur;
        xic_ref = imag(Zref) - mui;

        xrc_all = real(Z) - mur;
        xic_all = imag(Z) - mui;

        den = end_backgr - start_backgr;
        a = sum(xrc_ref.^2, 3) / den;
        b = sum(xrc_ref .* xic_ref, 3) / den;
        c = sum(xic_ref.^2, 3) / den;

        eps_reg = 1e-8;
        a = a + eps_reg;
        c = c + eps_reg;
        detS = max(a .* c - b.^2, eps_reg);

        mahala = sqrt(max((c .* xrc_all.^2 - 2*b.*xrc_all.*xic_all + a.*xic_all.^2) ./ detS, 0));
        mahalaCollapse = mahala(:,:,start_high + (0:num_collapse-1));

        mden = den + 1;
        Fmahala = (mden * (mden - 2) ./ (2 * (mden + 1) * (mden - 1))) .* (mahalaCollapse.^2);
        p_mahala = max(1 - fcdf(Fmahala, 2, mden - 2), realmin);

        results.mahalanobis = struct( ...
            'displayName', 'Mahalanobis (only when IQ data are available)', ...
            'showOnlyFinal', false, ...
            'leftMap', mahalaCollapse, ...
            'leftLabel', 'Mahalanobis distance', ...
            'pMap', p_mahala, ...
            'betaForDisplay', beta, ...
            'finalImage', sum(beta .* (p_mahala < (1e-4 / num_collapse)), 3), ...
            'rightLabel', 'Final BURST image', ...
            'numFrames', max(1, num_collapse));
    end
end

% ======================================================================
function m = nakagami_m_image(Y, nIter, eps0)
% Per-pixel Nakagami m estimation (kept from the supplied code)
    if nargin < 2 || isempty(nIter), nIter = 8; end
    if nargin < 3 || isempty(eps0),  eps0  = 1e-12; end

    Y = double(Y);
    [H,W,T] = size(Y);
    assert(T >= 2, 'Need at least 2 frames to fit the Nakagami baseline per pixel.');

    s2   = zeros(H,W);
    s4   = zeros(H,W);
    slog = zeros(H,W);

    for t = 1:T
        x = max(Y(:,:,t), eps0);
        x2 = x.^2;
        s2   = s2   + x2;
        s4   = s4   + x2.^2;
        slog = slog + log(x);
    end

    Omega = s2 / T;
    Ez2   = s4 / T;
    Varz  = max(Ez2 - Omega.^2, 0);

    m = Omega.^2 ./ max(Varz, eps0);
    m = max(m, 1e-3);

    mean_logZ = 2 * (slog / T);
    rhs = log(max(Omega, eps0)) - mean_logZ;

    for k = 1:nIter
        f  = log(m) - psi(m) - rhs;
        fp = 1./m - psi(1, m);
        step = f ./ fp;
        m = max(m - step, 1e-6);
        if max(abs(step(:))) < 1e-6
            break;
        end
    end
end

% ======================================================================
function out = clamp(x, lo, hi)
    out = min(max(x, lo), hi);
end

function s = frameLabel(frameIdx, start_high)
    s = sprintf('index %d (global frame %d)', frameIdx, start_high + frameIdx - 1);
end

function cLim = robustCLim(img)
    vals = img(:);
    vals = vals(isfinite(vals));
    if isempty(vals)
        cLim = [0 1];
        return;
    end
    lo = prctile(vals, 0.1);
    hi = prctile(vals, 99.9);
    if ~isfinite(lo) || ~isfinite(hi) || lo >= hi
        lo = min(vals);
        hi = max(vals);
        if lo == hi
            hi = lo + eps;
        end
    end
    cLim = [lo hi];
end
