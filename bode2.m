%==========================================================================
%
% bode2  Bode frequency response of dynamic systems.
%
%   bode2(L)
%   bode2(L,opts)
%   
% See also bode.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-11-10
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   L       - (1×1 tf) open loop transfer function
%   opts    - (OPTIONAL) (struct) plot options structure
%       • axis_font_size  - (str) axis label font size
%       • color           - (char or 1×3 double) line color (can be 
%                           specified as a name, short name, or RGB 
%                           triplet)
%       • freq_label      - (str) frequency axis label
%       • line_style      - (char) line style
%       • line_width      - (1×1 double) line width
%       • mag_label       - (str) magnitude plot vertical axis label
%       • mag_units       - (char) 'abs' or 'dB'
%                               --> 'abs' = absolute magnitude (gain) 
%                                           plotted
%                               --> 'dB' = magnitude (gain) plotted in 
%                                          decibels
%       • position        - (1×4 double) plot position [x,y,l,w]
%       • phase_label     - (str) phase plot vertical axis label
%       • title_str       - (str) title
%       • title_font_size - (str) title font size
%   units   - (OPTIONAL) (char) 'abs' or 'dB'
%               --> 'abs' = absolute magnitude (gain) plotted
%               --> 'dB' = magnitude (gain) plotted in decibels
%
% -------
% OUTPUT:
% -------
%   GM   	- (1×1 double) gain margin [- or dB]
%   PM   	- (1×1 double) phase margin [deg]
%   wcg   	- (1×1 double) phase crossover frequency [rad/s]
%   wcp     - (1×1 double) gain crossover frequency [rad/s]
%   h       - 
%
% -----------
% REFERENCES:
% -----------
%   [1] https://www.mathworks.com/matlabcentral/fileexchange/94115-magnitude-and-phase-of-a-transfer-function-mag_phase
%   [2] https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w
%
%==========================================================================
%function [bode_plot,mag_plot,phase_plot] = bode2(L,opts)
function [GM,PM,wcg,wcp] = bode2(L,opts)
    
    % ------------------------------------
    % Sets (or defaults) plotting options.
    % ------------------------------------
    
    % sets magnitude units (defaults to dB)
    if (nargin < 2) || ~isfield(opts,'mag_units')
        mag_units = 'dB';
    else
        mag_units = opts.mag_units;
    end

    % sets axis font size (defaults to 18)
    if (nargin < 2) || ~isfield(opts,'axis_font_size')
        axis_font_size = 14;
    else
        axis_font_size = opts.axis_font_size;
    end
    
    % sets title font size (defaults to 18)
    if (nargin < 2) || ~isfield(opts,'title_font_size')
        title_font_size = 14;
    else
        title_font_size = opts.title_font_size;
    end

    % sets text interpreter for labels (defaults to LaTeX)
    if (nargin < 2) || ~isfield(opts,'interpreter')
        interpreter = 'latex';
    else
        interpreter = opts.interpreter;
    end

    % sets title (defaults to "Bode Plot")
    if (nargin < 2) || ~isfield(opts,'title_str')
        if strcmpi(interpreter,'latex')
            title_str = "\textbf{Bode Plot}";
        else
            title_str = "Bode Plot";
        end
    else
        title_str = opts.title_str;
    end

    % sets magnitude plot vertical axis label (defaults to "Magnitude" with
    % the corresponding units)
    if (nargin < 2) || ~isfield(opts,'mag_label')
        if strcmpi(interpreter,'latex')
            if strcmpi(mag_units,'dB')
                mag_label = "Magnitude $[\mathrm{dB}]$";
            else
                mag_label = "Magnitude $[\mathrm{abs}]$";
            end
        else
            if strcmpi(mag_units,'dB')
                mag_label = "Magnitude [dB]";
            else
                mag_label = "Magnitude [abs]";
            end
        end
    else
        mag_label = opts.mag_label;
    end

    % sets line color (defaults to default MATLAB color)
    if (nargin < 2) || ~isfield(opts,'color')
        color = [0,0.4470,0.7410];
    else
        color = opts.color;
    end

    % sets line width (defaults to 1.5)
    if (nargin < 2) || ~isfield(opts,'line_width')
        line_width = 1.5;
    else
        line_width = opts.line_width;
    end

    % sets line style (defaults to solid line)
    if (nargin < 2) || ~isfield(opts,'line_style')
        line_style = '-';
    else
        line_style = opts.line_style;
    end
    
    % ------------------------------
    % Calculate data for Bode plots.
    % ------------------------------
    
    % poles and zeros
    p = pole(L);
    z = zero(L);

    % break points
    b = [p;z];

    %
    f = zeros(length(b),2);
    for i = 1:length(b)
        if isreal(b(i))
            f(i,1) = abs(b(i)/5);
            f(i,2) = abs(b(i)*5);
        else
            x = real(b(i));
            y = imag(b(i));
            wn = abs(real(sqrt(x^2+1j*2*y-y^2)));
            zeta = abs(x/wn);
            f(i,1) = wn*exp(-pi*zeta/2);
            f(i,2) = wn*exp(pi*zeta/2);
        end
    end
    
    % break points
%     x = real(b);
%     y = imag(b);
%     f = real(sqrt(x.^2+1j*2*y-y.^2));

    % do not want to consider 0 as a critical frequency
    f(all(~f,2),:) = [];

    % maximum frequency to draw until
    w_min = min(f(:,1))
    w_max = max(f(:,2))

    % corresponding exponents
    %N_min = floor(log10(w_min))
    %N_max = ceil(log10(w_max))

    N_min = round(log10(w_min)-1)
    N_max = round(log10(w_max)+1)

    % logarithmically spaced domain
    w = logspace(N_min,N_max,1000);

    % preallocates arrays
    mag = zeros(size(w));
    phase = zeros(size(w));
    
    % calculates magnitude and phase
    for i = 1:length(w)
        [mag(i),phase(i)] = mag_phase(L,1j*w(i));
    end

    % converts magnitudes to dB if need
    if strcmpi(mag_units,'dB')
        mag = 20*log10(mag);
    end

    % -------------------------------------------------
    % Gain and phase margins and crossover frequencies.
    % -------------------------------------------------

    % gain (w_co_g), and crossover frequency for phase (w_co);
    [GM,PM,wcg,wcp] = margin(L);

    % converts gain margin to dB
    if strcmpi(mag_units,'dB')
        GM = 20*log10(GM);
    end

    % ---------------
    % Overall figure.
    % ---------------
    
    bode_plot = 0;
    axes_array = tight_subplot(2,1,0.03,0.1,0.1);

    % ---------------
    % Magnitude plot.
    % ---------------

    % initializes axes and stores corresponding handle
    axes(axes_array(1));
    mag_plot.axes = gca;

    % plots magnitude and stores corresponding handle
    if strcmpi(mag_units,'dB')
        mag_plot.line = semilogx(w,mag,'linewidth',line_width,'color',...
            color,'linestyle',line_style);
    else
        mag_plot.line = loglog(w,mag,'linewidth',line_width,'color',...
            color,'linestyle',line_style);
    end

    % plot formatting
    grid on;

    % vertical axis label
    ylabel(mag_label,'interpreter',interpreter,'fontsize',axis_font_size);
    
    % title of overall diagram (done on magnitude subplot)
    title(title_str,'interpreter',interpreter,'fontsize',title_font_size);
    
    % turns of horizontal axis tick labels (frequency will be labeled on
    % phase plot)
    set(mag_plot.axes,'XTickLabel',[]);


    % -----------
    % Phase plot.
    % -----------

    % initializes axes and stores corresponding handle
    axes(axes_array(2));
    phase_plot.axes = gca;

    % plots phase and stores corresponding handle
    phase_plot.line = semilogx(w,phase,'linewidth',line_width,'color',...
        color,'linestyle',line_style);
    
    % plot formatting
    grid on;

    % sets y axis for phase plot (limits plot to be between minimum and
    % maximum phases, rounded to the nearest multiple of 45°)
    min_phase = 45*round(min(phase)/45)
    max_phase = 45*round(max(phase)/45)
    ylim([min_phase,max_phase]);

    % tick marks spaced 45° apart
    yticks(min_phase:45:max_phase);

    % axis labels
    xlabel('Frequency [$\mathrm{rad/s}$]','interpreter',interpreter,...
        'fontsize',axis_font_size);
    ylabel('Phase [$^{\circ}$]','interpreter',interpreter,'fontsize',...
        axis_font_size);

    % wrapping to 360?

    %----------------------------------------------------------------------
    % mag_phase
    %
    % Magnitude and phase of a transfer function (i.e. linear system) at a
    % specific point in the frequency domain.
    %
    % Modified from:
    % https://www.mathworks.com/matlabcentral/fileexchange/94115-magnitude-and-phase-of-a-transfer-function-mag_phase
    %----------------------------------------------------------------------
    %
    % INPUT:
    %   sys     - (1×1 tf) continuous transfer function
    %   x       - (1×1 complex double) location in frequency domain
    %
    % -------
    % OUTPUT:
    % -------
    %   mag     - (1×1 double) magnitude of transfer function at s = x
    %   phase   - (1×1 double) phase of transfer function at s = x [deg]
    %
    %----------------------------------------------------------------------
    function [mag,phase] = mag_phase(sys,x)
        mag = norm(evalfr(sys,x));
        phase = atan2d(imag(evalfr(sys,x)),real(evalfr(sys,x)));
    end  

    %----------------------------------------------------------------------
    % tight_subplot
    %
    % Creates "subplot" axes with adjustable gaps and margins. 
    %
    % Modified from:
    % https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w
    %----------------------------------------------------------------------
    %
    % INPUT:
    %   Nh      - (1×1 double) number of axes in hight (vertical direction)
    %   Nw      - (1×1 double)number of axes in width (horizontaldirection)
    %   gap     - (1×1 double) gaps between the axes in normalized units
    %             (0...1) or [gap_h gap_w] for different gaps in height and
    %             width 
    %   marg_h  - (1×1 double) margins in height in normalized units 
    %             (0...1) or [lower upper] for different lower and upper 
    %             margins 
    %   marg_w  - (1×1 double) margins in width in normalized units (0...1)
    %             or [left right] for different left and right margins 
    %
    % OUTPUT:
    %   ha      - (Nh×Nw Axes) array of handles of the axes objects 
    %             starting from upper left corner, going row-wise as in
    %             subplot
    %   pos     - (Nh×Nw cell) positions of the axes objects
    %
    %----------------------------------------------------------------------
    function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
        if nargin<3; gap = .02; end
        if nargin<4 || isempty(marg_h); marg_h = .05; end
        if nargin<5; marg_w = .05; end
        if numel(gap)==1
            gap = [gap gap];
        end
        if numel(marg_w)==1
            marg_w = [marg_w marg_w];
        end
        if numel(marg_h)==1
            marg_h = [marg_h marg_h];
        end
        axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
        axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
        py = 1-marg_h(2)-axh; 
        ii = 0;
        for ih = 1:Nh
            px = marg_w(1);
            for ix = 1:Nw
                ii = ii+1;
                ha(ii) = axes('Units','normalized', ...
                    'Position',[px py axw axh], ...
                    'XTickLabel','', ...
                    'YTickLabel','');
                px = px+axw+gap(2);
            end
            py = py-axh-gap(1);
        end
        if nargout > 1
            pos = get(ha,'Position');
        end
        ha = ha(:);
    end

end