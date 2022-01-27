%
% Copyright (c) 2022 Bundesanstalt für Wasserbau
%
% This file is part of VIPER.
%
% VIPER is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% VIPER is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with VIPER.  If not, see <http://www.gnu.org/licenses/>.
%


function varargout = viper(varargin)
% VIPER M-file for viper.fig
% VIPER, by itself, creates a new VIPER or raises the existing
% singleton*.
% H = VIPER returns the handle to a new VIPER or the handle to
% the existing singleton*.
% VIPER('CALLBACK',hObject,eventData,handles,...) calls the local
%   function named CALLBACK in VIPER.M with the given input arguments.
%   VIPER('Property','Value',...) creates a new VIPER or raises 
%   the existing singleton*.  Starting from the left, property value pairs 
%   are applied to the GUI before adv_OpeningFunction gets called.  An
%   unrecognized property name or invalid value makes property 
%   application stop.All inputs are passed to viper_OpeningFcn via 
%   varargin.
%   *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%   instance to run (singleton)".
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help viper
% Last Modified by GUIDE v2.5 26-Jan-2022 15:00:44
%%%%%%%%%%% ----- Begin initialization code - DO NOT EDIT ---- %%%%%%%%%%%
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @viper_OpeningFcn, ...
                   'gui_OutputFcn',  @viper_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%%%%%%%%%%% ----- End initialization code - DO NOT EDIT  ----- %%%%%%%%%%%


% --- Executes just before viper is made visible.
function viper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% varargin   command line arguments to viper (see VARARGIN)
% Choose default command line output for viper
handles.output = hObject;

% -------------------------------------------------------------------------
%%% Initial dialog asking for basic user setting:
% -------------------------------------------------------------------------
input = questdlg(sprintf('%s\n\n%s',...
    'Staring VIPER...','Please specify GUI display units:'),...
    'VIPER defaults',...
    'm, m/s, ...','cm, cm/s, ...','cm, cm/s, ...');%
if strcmp(input,'m, m/s, ...')==1
    user_units = 1000;
elseif strcmp(input,'cm, cm/s, ...')==1
    user_units = 10;
else % strcmp(input,'Skip file')==1 || empty
    user_units = 10;
end

% -------------------------------------------------------------------------
%%% Initializing veriables...
% -------------------------------------------------------------------------
wbar = waitbar(0.01,'Initializing VIPER, please wait...',...
    'WindowStyle','modal');

% -------------------------------------------------------------------------
%%% Log and log convention
% -------------------------------------------------------------------------
handles.log = {};%size=[*growing*,1]
%--> new selection
%  > filename or name of task starting
%  ! error == task not done
%  # done task
%    * warning/information on task

% -------------------------------------------------------------------------
%%% Fields for filename of each time-series, size: [numfiles,1]
% -------------------------------------------------------------------------
%path name of time-series source file
handles.fpath    = {};
%source file name without .ext
%NOT USED to reread file but 
%USED AS root filename for exporting, such as .par file, image file,...
handles.fname    = {};% virtual files get (X,Y,Z)_sourcefilename
%fullfilenamewithpath:
%NOT USED to reread file but
%USED to check whether a files has already been loaded
%USED as file specifier in exported tables (mat)
handles.fullname = {};% virtual files get [fpath, (X,Y,Z)_sourcefilename]
%displayname
%USED only in listbox: [filename '   [' Pathname ']' ]
handles.dispname = {};% virtual files get [(X,Y,Z)_sourcefilename ...



% -------------------------------------------------------------------------
% Velocity time series DATASET
% * Content initially filled in: gui_addfile_import_ts_data
% -------------------------------------------------------------------------
handles.ts = struct([]);% size: [numfiles,1]
% ...with following fields:
% ts.intl          -- data used internally, see handles.ts_defaults.intl
% ts.hdr_nfo       -- .hdr data, see handles.ts_defaults.hdr_nfo
% *raw_tslength = raw ts-length after filling intermittent missing samples
% ts.raw_t         -- time axis for raw_veldata
% ts.raw_veldata   -- raw velocities, size [raw_tslength,5]: [u,v,w,w1,w2]
%           *missing samples: nan
% ts.raw_vals      -- 1==valid-in-raw (real sample), size [raw_tslength,1]
% ts.raw_snrcor    -- raw SNR & COR values for each beam and sample
% ts.raw_corsnr_xt -- size [raw_tslength,4]: [minCOR avgCOR minSNR avgSNR]
% ts.efilt_res     -- errfiltering results,see handles.ts_defaults.efilt_res
% ts.vals          -- 1==valid after error filtering, size [raw_tslength,1]
% ts.p             -- parameterset, see handles.ts_defaults.p
% ts.enable_coords -- marking whether internal? coordinates are enabled
% Data after error filtering - statistical values
% ts.valu_data     -- statistical values data, handles.ts_defaults.valu_data
% ts.valu_enable   -- marking whether values are enabled
% Data after error filtering - statistical functions
% Using only original samples
% ts.corsnr_dist   -- COR&SNR distribution of valid samples
% ts.hist          -- velocity histograms of valid data
% Using original and artificial samples:
% *seri_tslength = ts-length shortened (without initial and trailing erroneous samples)
% ts.seri_data   -- series data, see handles.ts_defaults.seri_data
% ts.seri_enable -- marking whether series are enabled
% ts.seri_x       -- independent (== x) variables for series data
%           ts.seri_x.t
%           ts.seri_x.lag
%           ts.seri_x.tlag
%           ts.seri_x.f



% -------------------------------------------------------------------------
% Defaults for such ts-field:
% - default data-PARAMETERSET
% - default (empty) fields for storing error filtering results
% - default hdr-info (Nortek ADV)
% - default fields for use within a session internally
% - default fields for results of statistical functions and values
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Time-series PARAMETERSET supported by present VIPER version!
% * each pargroup can have fields - BUT ONLY ONE LEVEL (except errfilt)!
% * variable classes allowed: char, double, cell-of-strings!!!
%parset_properties
handles.ts_defaults.p.parset_props.ver        = 1.3;
handles.ts_defaults.p.parset_props.id         = NaN;
%datfile_properties
handles.ts_defaults.p.dat_props.type        = 'vno';%'vno' | 'vec' | 'usr'
handles.ts_defaults.p.dat_props.column_t    = 0;
handles.ts_defaults.p.dat_props.column_u    = 3;
handles.ts_defaults.p.dat_props.frq          = 99;
handles.ts_defaults.p.dat_props.probe_cfg   = 'side-looking';
handles.ts_defaults.p.dat_props.probe_nvr   = {};
% handles.ts_defaults.p.dat_props.probe_vr_np = {};
%coord_parameters - note this also: handles.ts_defaults.intl.coordsys
handles.ts_defaults.p.coordsys.unit       = 'm';
handles.ts_defaults.p.coordsys.file_coords    = [NaN,NaN,NaN];
handles.ts_defaults.p.coordsys.measpt_offset  = [NaN,NaN,NaN];
handles.ts_defaults.p.coordsys.measpt_coords  = [NaN,NaN,NaN];
handles.ts_defaults.p.coordsys.dir_probe_u    = '+x';
% Error filtering parameters: (Default = pst, but lambda is ts-dependent!!)
handles.ts_defaults.p.errfilt.used_velcomp  = [1;1;1];
handles.ts_defaults.p.errfilt.used_types    = {'pst'};%empty == no filter
% * see supported types: errf_types_supported
handles.ts_defaults.p.errfilt.typ(1).name =...
    'Phase-space threshold despiking';
handles.ts_defaults.p.errfilt.typ(1).used_comp = [1;1;1];%used_comp
handles.ts_defaults.p.errfilt.typ(1).lambda = NaN;%for pst always nan!
handles.ts_defaults.p.errfilt.typ(1).hipass_frq = [];%hipass_frq []==no highpass (never NaN! or <=0)

% -------------------------------------------------------------------------
% Default fields storing errfilt results
handles.ts_defaults.efilt_res.valmarker     = [];
handles.ts_defaults.efilt_res.errmarker     = [];
handles.ts_defaults.efilt_res.eftyp         = cell(0);
handles.ts_defaults.efilt_res.badpos_cc_vct = [];
handles.ts_defaults.efilt_res.x_ell_vct     = [];
handles.ts_defaults.efilt_res.y_ell_vct     = [];
handles.ts_defaults.efilt_res.badpos_cc_pst = [];
handles.ts_defaults.efilt_res.x_ell_pst     = [];
handles.ts_defaults.efilt_res.y_ell_pst     = [];

% -------------------------------------------------------------------------
% Default field storing ADV HDR PARAMETERS used in present VIPER version!
handles.ts_defaults.hdr_nfo.existing      = false;
handles.ts_defaults.hdr_nfo.sampling_rate = '';
handles.ts_defaults.hdr_nfo.velocity_rnge = '';
handles.ts_defaults.hdr_nfo.probe_cfg     = '';
handles.ts_defaults.hdr_nfo.transmit_lgth = '';
handles.ts_defaults.hdr_nfo.sampling_volu = '';
handles.ts_defaults.hdr_nfo.power_level   = '';
handles.ts_defaults.hdr_nfo.coord_system  = '';
% handles.ts_defaults.hdr_nfo.sound_speed   = '';
% handles.ts_defaults.hdr_nfo.salinity      = '';
handles.ts_defaults.hdr_nfo.serial_no     = '';
handles.ts_defaults.hdr_nfo.ensemble_column  = '';
% handles.ts_defaults.hdr_nfo.num_measurements = [];

% -------------------------------------------------------------------------
% Default fields used for aktive viper session internally
handles.ts_defaults.intl.cmp_probe_uvw = {'up' 'vp' 'wp'};%cmp_probe_uvw col positions
handles.ts_defaults.intl.cmp_probe_np = {'n' 'n' 'p'};%par or norm to transmitter
handles.ts_defaults.intl.mean_correctn = [0 0 0];%correction after interp
%coordsys stored in gui_units: (.p.units is in meters)
%get .intl-values by multiplying .p-values by handles.gui.coordunitscale
handles.ts_defaults.intl.coordsys = handles.ts_defaults.p.coordsys;
%used segment of time series (for seri) after err.filt.
handles.ts_defaults.intl.used_ts_seg = [];

% -------------------------------------------------------------------------
% Default fields for statistical features (series)
% * see handles.field_props.seri for properties of seri-s ...
handles.ts_defaults.seri_data   = cell(8,1);
% Content:
% {1}: velocity fluctuation ts (contains nans at beginning and end)
% {3}: autocorr function
% {4}: cumulative mean function
% {8}: spectrum
handles.ts_defaults.seri_enable = false(8,1);%whether seri is enabled

% -------------------------------------------------------------------------
% Default fields for statistical features (values)
% * see handles.field_props.valu for properties of valu-s ...
% Caution: adjust size with .field_props below
handles.ts_defaults.valu_data   = nan(1,46);
% Content: see handles.field_props.valu_roots
handles.ts_defaults.valu_enable = false(1,46);



% -------------------------------------------------------------------------
% General properties of current Version of VIPER 
% - supported error filtering methods by viper
% - user-definable defaults
% - velocity ranges
% -------------------------------------------------------------------------
% Supported errfilt methods & which is enabled for which data_type
handles.viper_props.errf_types_supported = ...
    {'crop'   'Time-series cropping';...
     'pst'    'Phase-space threshold despiking';...
     'cormin' 'Correlation filter using minimum of beams';...
     'coravg' 'Correlation filter using average of beams'};
handles.viper_props.errf_types_enable_matrix = ...
    {'vno', [ true; true; true; true ];...
     'vec', [ true; true; true; true ];...
     'usr', [ true; true; false;false]};
% % How to get erffilt. methods that are enabled for dat_type:
% %  f_idx  = handles.sorted_index(get(handles.ui_listbox1,'Value'));
% %  ef_vals = strcmpi(...
% %      handles.viper_props.errf_types_enable_matrix(:,1),...
% %      handles.ts( f_idx ).p.dat_props.type );
% %  errf_types_enabled = ...
% %      handles.viper_props.errf_types_supported( ...
% %      handles.viper_props.errf_types_enable_matrix{ef_vals ,2}, 1);

% -------------------------------------------------------------------------
% USER-DEFINED DEFAULTS - changing during user interaction
% datfile properties default:
handles.viper_props.userdef.dat_props.column_t = 1;
handles.viper_props.userdef.dat_props.column_u = 4;
% Defaults for getting coordinates
handles.viper_props.userdef.cds_fmt = '(%x%;%y%;%z%)_20131022150823';
handles.viper_props.userdef.cds_unit = 1;%1000==m, 10==cm, 1==mm
% Defaults spectrum line add x-y pairs
handles.viper_props.userdef.plot_line1 = 0;
handles.viper_props.userdef.plot_line1_x = [1 1000];
handles.viper_props.userdef.plot_line1_y = [1000 1];
handles.viper_props.userdef.plot_line53 = 0;
handles.viper_props.userdef.plot_line53_x = [1 1000];
handles.viper_props.userdef.plot_line53_y = [1000 0.01];
handles.viper_props.userdef.plot_line3 = 0;
handles.viper_props.userdef.plot_line3_x = [1 10];
handles.viper_props.userdef.plot_line3_y = [1000 1];

% -------------------------------------------------------------------------
% Velocity ranges:
handles.viper_props.velrange_limits_vno = ...
    [0.03, 0.26, 0.08 ;...
     0.10, 0.44, 0.13 ;...
     0.30, 0.94, 0.27 ;...
     1.00, 1.88, 0.54 ;...
     2.50, 3.28, 0.94 ;...
     4.00, 5.25, 1.5 ];
handles.viper_props.velrange_limits_vec = ...
    [0.01, 0.05, 0.02 ;...
     0.10, 0.35, 0.10 ;...
     0.30, 0.81, 0.23 ;...
     1.00, 2.10, 0.60 ;...
     2.00, 3.50, 1.00 ;...
     4.00, 5.25, 1.50  ;...
     7.00, 8.75, 2.50 ];

% -------------------------------------------------------------------------
%%% Properties of statistical features field in this version
% * ALWAYS adjust the size of the fields at changes!
% -------------------------------------------------------------------------
% names of .valu_data content (title of columns at export)
% please place values describing u-v-w next to each other
handles.field_props.valu_name_roots = {...
    'u_mean',        'v_mean',        'w_mean',...%1-3
    'u_var',         'v_var',         'w_var',...%4-6
    'u_std',         'v_std',         'w_std',...%7-9
    'U',             'TKE',           'U_rms',...%10-12
    'r_uv',          'r_uw',          'r_vw',...%13-15
    'valid ratio',   'COR_avg',       'SNR_avg',...%16-18
    'COR_min',       'COR_q25',       'COR_median',...%19-21
    'COR_q75',       'COR_max',       'valid duration',...%22-24
    'T_I,u',         'T_I,v',         'T_I,w',...%25-27
    'my_4,u',        'my_4,v',        'my_4,w',...%28-30
    's_95%{u_mean}', 's_95%{v_mean}', 's_95%{w_mean}',...%31-33
    'e_95%{u_mean}', 'e_95%{v_mean}', 'e_95%{w_mean}',...%34-36
    's_95%{u_var}',  's_95%{v_var}',  's_95%{w_var}',...%37-39
    'e_95%{u_var}',  'e_95%{v_var}',  'e_95%{w_var}',...%40-42
    'I_u',           'I_v',           'I_w',...%43-45
    'I'};%46
handles.field_props.valu_units_cm = {...
    ', cm/s',   ', cm/s',    ', cm/s',...
    ', cm2/s2', ', cm2/s2',  ', cm2/s2',...
    ', cm/s',   ', cm/s',    ', cm/s',...
    ', cm/s',   ', cm2/s2',  ', cm/s',...
    ', cm2/s2', ', cm2/s2',  ', cm2/s2',...
    ', %',      ', %',       ', dB',...
    ', %',      ', %',       ', %',...
    ', %',      ', %',       ', s',...
    ', s',      ', s',       ', s',...
    ', cm4/s4', ', cm4/s4',  ', cm4/s4',...
    ', cm/s',   ', cm/s',    ', cm/s',...
    ', %',      ', %',       ', %',...
    ', cm2/s2', ', cm2/s2',  ', cm2/s2',...
    ', %',      ', %',       ', %',...
    ', -',      ', -',       ', -',...
    ', -'};
handles.field_props.valu_units_m = {...
    ', m/s',    ', m/s',     ', m/s',...
    ', m2/s2',  ', m2/s2',   ', m2/s2',...
    ', m/s',    ', m/s',     ', m/s',...
    ', m/s',    ', m2/s2',   ', m/s',...
    ', m2/s2',  ', m2/s2',   ', m2/s2',...
    ', %',      ', %',       ', dB',...
    ', %',      ', %',       ', %',...
    ', %',      ', %',       ', s',...
    ', s',      ', s',       ', s',...
    ', m4/s4',  ', m4/s4',   ', m4/s4',...
    ', m/s',    ', m/s',     ', m/s',...
    ', %',      ', %',       ', %',...
    ', m2/s2',  ', m2/s2',   ', m2/s2',...
    ', %',      ', %',       ', %',...
    ', -',      ', -',       ', -',...
    ', -'};
handles.field_props.valu_unitscalingpower = [...
    1, 1, 1,    2, 2, 2,    1, 1, 1,    1, 2, 1,    2, 2, 2,    0, 0, 0,...
    0, 0, 0,    0, 0, 0,    0, 0, 0,    4, 4, 4,    1, 1, 1,    0, 0, 0,...
    2, 2, 2,    0, 0, 0,    0, 0, 0,    0 ];%m, m/s ->1, m2/s2 -> 2, else: 0==not to scale
%Which vel component described: 1: u, 2: v, 3: w, 0 - not specific
handles.field_props.valu_uvw = [...
    1, 2, 3,    1, 2, 3,    1, 2, 3,    0, 0, 0,    0, 0, 0,    0, 0, 0,...
    0, 0, 0,    0, 0, 0,    1, 2, 3,    1, 2, 3,    1, 2, 3,    1, 2, 3,...
    1, 2, 3,    1, 2, 3,    1, 2, 3,    0 ];
%Sign of which vel component should be adjusted at vector rotation
handles.field_props.valu_adjrotsign = [...
    1, 1, 1,    0, 0, 0,    0, 0, 0,    0, 0, 0,    0, 0, 0,    0, 0, 0,...
    0, 0, 0,    0, 0, 0,    0, 0, 0,    0, 0, 0,    0, 0, 0,    0, 0, 0,...
    0, 0, 0,    0, 0, 0,    0, 0, 0,    0 ];
%names of groups describing u-v-w groups of .valu_uvw (= non-zero elements)
handles.field_props.valu_groups_roots = {...
    'mean';...
    'variance';...
    'st. deviation';...
    'integral time scale';...
    '4th central moment';...
    's_95%(mean)';...
    'e_95%(mean)';...
    's_95%(variance)';...
    'e_95%(variance)';...
    'turbulence intensity'};
handles.field_props.valu_groups_cm = {...
    ', cm/s';    ', cm2/s2';  ', cm/s';...
    ', s';       ', cm4/s4';  ', cm/s';...
    ', %';       ', cm2/s2';  ', %';...
    ', -'};
handles.field_props.valu_groups_m = {...
    ', m/s';    ', m2/s2';    ', m/s';...
    ', s';      ', m4/s4';    ', m/s';...
    ', %';      ', m2/s2';    ', %';...
    ', -'};
% -------------------------------------------------------------------------
% names of .seri_data content (title of columns at export)
% handles.field_props.seri_names = {...
%     'u'' , cm/s',     'v'' , cm/s',     'w'' , cm/s',...
%     '-',              '-',              '-',...
%     'R_uu , -',       'R_vv , -',       'R_ww , -',...
%     'u_mean , cm/s',  'v_mean , cm/s',  'w_mean , cm/s',...
%     '-',              '-',              '-',...
%     '-', '-',  '-',...
%     '-', '-',  '-',...
%     'S_uu , cm2/s',   'S_vv , cm2/s',   'S_ww , cm2/s'};%
%Sign of which component should be adjusted at vector rotation
handles.field_props.seri_adjrotsign = [...
    1, 0, 0, 0, 0, 0, 0, 0   ];



% -------------------------------------------------------------------------
% GUI Properties of plots - adjust number of rows!!!
% -------------------------------------------------------------------------
handles.plots.fig_h =  nan(11,1);%containers to store 11 object handles
%Names of plot: not used container: keep empty 
handles.plots.names =...
    {'Statistical values';...
     'COR & SNR distributions (line: avg. | bars: min. of beams)';...
     'Raw vs. error-filtered data';...
     'Histograms';...
     'Velocity spectra';...
     'Cumulative mean';...
     [];...
     [];...
     'Autocorrelation function';...
     [];...
     []};
%Suffix of plot at exporting: not used container: keep empty
handles.plots.expname_suffix =...
    {'_stat_table';...
     '_corsnrdist';...
     '_ts_errfilt';...
     '_histograms';...
     '_spectra';...
     '_cum_mean';...
     [];...
     [];...
     '_autocorr';...
     [];...
     []};
handles.plots.cmp_color=...
    {'b', 'r', 'm'};
handles.plots.cmp_names={...
     [],  [],  [];...
     [],  [],  [];...
    'u', 'v', 'w';...
    'u''','v''','w''';...
    'S_{uu}', 'S_{vv}', 'S_{ww}';...
    'u^m', 'v^m', 'w^m';...
     [],  [],  [];...
     [],  [],  [];...
    'r_{uu}', 'r_{vv}', 'r_{ww}';...
     [],  [],  [];...
     [],  [],  []};
%which seri is used for plot type: 0=raw, nan=no reference
handles.plots.reference_seri = ...
    [ nan; nan; 1; nan; 8; 4; nan; nan; 3; nan; nan];
%%% fixed figures:
handles.plots.fixed_figs = 1;
%Axes properties: %keep empty if no that many axes
handles.plots.axs_h = nan(11,4);
%Number of axes in plots
handles.plots.numaxs = ...
    [ 0; 3; 2; 3; 1; 1; 0; 0; 1; 0; 0];
%Userdefined Axislimit defaults in cm.
%WARNING: also used to determine number of axes in plots
handles.plots.axlims_usrdef ={...
    [], [], [], [];...
    [   0 100;   0 100], [ 60 90;0 100], [ 0 100; 0 25], [];...
    [   0 300;-200 200], [   0 300;-200 200], [], [];...
    [-100 100;   0  20], [-100 100;   0  20], [ -100 100;   0  20], [];...
    [  -4   4;  -4   4], [], [], [];...
    [   0 300;-100 100], [], [], [];...
    [], [], [], [];...
    [], [], [], [];...
    [   0  10;-0.5   1], [], [], [];...
    [], [], [], [];...
    [], [], [], []};%NON-nan defaults(!) and in cm, cm/s [x1 x2 ; y1 y2]!!
% Enable asx limits adjusting
handles.plots.adjustable_axlims = logical(...
    [ 0; 0; 1; 1; 1; 1; 0; 0; 1; 0; 0]);
% Axislimits to use: 0==auto, 1==usrdef, default==0
handles.plots.axlims_to_use = zeros(11,4);
% Initialize axs limits minmax values:
% for non-adjustable: use handles.plots.axlims_usrdef (default)
handles.plots.axlims_minmax = cell(11,4);%nan enabled
% for adjustables: use here nan, and calculate at data load/calc
for p_i=1:size( handles.plots.axs_h, 1)
    if handles.plots.adjustable_axlims(p_i)==1
        for a_i=1:size( handles.plots.axs_h, 2)
            if isempty(handles.plots.axlims_usrdef{p_i,a_i})==0
                handles.plots.axlims_minmax{p_i,a_i} =...
                    nan(size( handles.plots.axlims_usrdef{p_i,a_i} ));
            end
        end
    else
        handles.plots.axlims_minmax(p_i,:) =...
            handles.plots.axlims_usrdef(p_i,:);
    end
end
handles.plots.axlims_unitpower = {...
    [  0   0;  0   0], [], [], [];...
    [  0   0;  0   0], [], [  0   0;  0   0], [];...
    [  0   0;  1   1], [  0   0;  1   1], [], [];...
    [  1   1;  0   0], [  1   1;  0   0], [  1   1;  0   0], [];...
    [  0   0;  2   2], [], [], [];...
    [  0   0;  1   1], [], [], [];...
    [], [], [], [];...
    [], [], [], [];...
    [   0   0;   0   0], [], [], [];...
    [], [], [], [];...
    [], [], [], []};%[x1 x2 ; y1 y2]: 0==not scaled
%Axis Labels: %rows=plot, column=subplot
handles.plots.axs_names =...
    {'Statistics',        '-',              '-',           '-';...
     'COR hist.',  'COR quartiles',     'SNR hist.',  '-';...
     'Velocity raw',   'Velocity valid',           '-',           '-';...
     'Velocity u',    'Velocity v',    'Velocity w', '';...
     'Spectral power',  '-',            '-',           '-';...
     'Mean of velocity','-',            '-',           '-';...
     '-',                '-',            '-',           '-';...
     '-',                '-',            '-',           '-';...
     'Autocorrelation',     '-',            '-',           '-';...
     '-',                '-',            '-',           '-';...
     '-',              '-',              '-',           '-'};
handles.plots.axs_label_x_roots =...
    {'',              '',              '',           '';...
     'Relative frequency',    '',      'Relative frequency', '';...
     'Time',          'Time',          '',           '';...
     'Velocity u',    'Velocity v',    'Velocity w', '';...
     'Frequency',     '',              '',           '';...
     'Time',          '',              '',           '';...
     '',              '',              '',           '';...
     '',              '',              '',           '';...
     'Time lag',      '',              '',           '';...
     '',              '',              '',           '';...
     '',              '',              '',           ''};
handles.plots.axs_label_y_roots =...
    {'',                '',            '',           '';...
     'Correlation (COR)',  '',     'Signal-to-noise (SNR)', '';...
     'Velocity, raw',   'Velocity, valid',           '',           '';...
     'Relative frequency','Relative frequency','Relative frequency',  '';...
     'Spectral power',  '',            '',           '';...
     'Mean of velocity','',            '',           '';...
     '',                '',            '',           '';...
     '',                '',            '',           '';...
     'Correlation',     '',            '',           '';...
     '',                '',            '',           '';...
     '',                '',            '',           ''};
handles.plots.axs_label_x_units_cm =...
    {'',        '',        '',        '';...
     ', %',     '',        ', %',     '';...
     ', s',     ', s',     '',        '';...
     ', cm/s',  ', cm/s',  ', cm/s',  '';...
     ', Hz',    '',        '',        '';...
     ', s',     '',        '',        '';...
     '',        '',        '',        '';...
     '',        '',        '',        '';...
     ', s',     '',        '',        '';...
     '',        '',        '',        '';...
     '',        '',        '',        ''};
handles.plots.axs_label_x_units_m =...
    {'',        '',        '',        '';...
     ', %',     '',        ', dB',    '';...
     ', s',     ', s',     '',        '';...
     ', m/s',   ', m/s',   ', m/s',   '';...
     ', Hz',    '',        '',        '';...
     ', s',   '',        '',        '';...
     '',        '',        '',        '';...
     '',        '',        '',        '';...
     ', s',     '',        '',        '';...
     '',        '',        '',        '';...
     '',        '',        '',        ''};
handles.plots.axs_label_y_units_cm =...
    {'',        '',        '',        '';...
     ', %',     '',        ', dB',    '';...
     ', cm/s',  ', cm/s',  '',        '';...
     ', %',     ', %',     ', %',     '';...
     ', cm^2/s','',        '',        '';...
     ', cm/s',  '',        '',        '';...
     '',        '',        '',        '';...
     '',        '',        '',        '';...
     ', -',     '',        '',        '';...
     '',        '',        '',        '';...
     '',        '',        '',        ''};
handles.plots.axs_label_y_units_m =...
    {'',        '',        '',        '';...
     ', %',     '',        ', dB',    '';...
     ', m/s',   ', m/s',   '',        '';...
     ', %',     ', %',     ', %',     '';...
     ', m^2/s', '',        '',        '';...
     ', m/s',   '',        '',        '';...
     '',        '',        '',        '';...
     '',        '',        '',        '';...
     ', -',     '',        '',        '';...
     '',        '',        '',        '';...
     '',        '',        '',        ''};
% Enable asx labels adjusting
handles.plots.adjustable_axlabels = logical(...
    [ 0; 0; 1; 1; 1; 1; 0; 0; 1; 0; 0]);
% Enable asx labels adjusting
handles.plots.adjustable_axlegend = logical(...
    [ 0; 0; 1; 0; 1; 1; 0; 0; 1; 0; 0]);
% Enable asx legend on|off
handles.plots.axs_legend_on = logical(...
    [ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]);
%Axs width ratios (ratio of axes width to total width) & axes aspect ratios
handles.plots.siz.axs_wr =[...
    0.33, 0.67,   0,   0;...
    0.31, 0.07, 0.31,  0;...
       1,    1,   0,   0;...
    0.33, 0.33, 0.33,  0;...
       1,    0,   0,   0;...
       1,    0,   0,   0;...
       0,    0,   0,   0;...
       0,    0,   0,   0;...
       1,    0,   0,   0;...
       0,    0,   0,   0;...
       0,    0,   0,   0];
handles.plots.siz.axs_ar =[...
    1.16, 0.56,   0,   0;...
       1,  nan,   1,   0;...
    0.27, 0.27,   0,   0;...
       1,    1,   1,   0;...
       1,    0,   0,   0;...
    0.27,    0,   0,   0;...
       0,    0,   0,   0;...
       0,    0,   0,   0;...
    0.27,    0,   0,   0;...
       0,    0,   0,   0;...
       0,    0,   0,   0];
%Default sizes and font name:
handles.plots.ft_name   = 'Times New Roman';%'Arial'
handles.plots.siz.ft_pt = 11;%Font size
%Print: plot widths (the largest uses full width of siz.used_screen_wr)
handles.plots.siz.printwidth_mm =...
    [ 175; 175; 175; 175; 85; 175; 0; 0; 175; 0; 0];%Print width of 
%On-screen: Used width ratio of free screen width for plots
handles.plots.siz.used_screen_wr = 1;



waitbar( 1/10, wbar)

% -------------------------------------------------------------------------
% GUI-PROPERTIES:
% system-dependent properties (GUI & plot positioning, font size, etc)
% -------------------------------------------------------------------------
set(0,'Units','pixels');
handles.gui.pos_screen = get( 0, 'screensize');
% Warning: Do NOT USE here absolute positions of main_gui - they get
% changed later
% -------------------------------------------------------------------------
% MAIN GUI & SCREENSIZE
fig1out = get( handles.figure1, 'OuterPosition');
fig1pos = get( handles.figure1, 'Position');
% TASKBAR SIZE
% Automated - WITH Figure maximization  - untested
% create fig_plt here
% set( fig_plt, 'WindowState', 'maximized')
% pause(0.25)
% pos_plt_mx = get( fig_plt, 'Position');%!MAXIMIZED Plot-figure
% siz_taskbar = pos_plt_mx - 1;
% Manually - Without figure maximization - leave 60 pix at bottom
handles.gui.siz_taskbar = 60;%bottom
% GUI-borders: left - bottom - right - top: (no toolbar)
dp_gui = fig1out - fig1pos;
siz_borders_gui =...
    [-dp_gui(1) -dp_gui(2) dp_gui(3)-(-dp_gui(1)) dp_gui(4)-(-dp_gui(2))];
%space needed at top of window:
handles.gui.siz_ontop = siz_borders_gui(4) - siz_borders_gui(2);
%Maximal space for additional figures:  left bottom | width heihgt
handles.gui.lb_min = [ 1, 1 + handles.gui.siz_taskbar ];
handles.gui.wh_max = [...
    handles.gui.pos_screen(3) - fig1out(3),...
    handles.gui.pos_screen(4) - handles.gui.siz_taskbar - ...
    handles.gui.siz_ontop ];
% -------------------------------------------------------------------------
% GUI-PROPERTIES:
% user-defined properties:
% -------------------------------------------------------------------------
% Change basic GUI-units, variables (incl. labels) dependent on units
% handles.gui.units_id: 1==mm, 10==cm, 100==dm, 1000==m
handles.gui.units_id = user_units;
if handles.gui.units_id == 10%== cm, cm/s
    handles.gui.unit_length_str = 'cm';
    %get .intl-values by multiplying .p-values by:
    handles.gui.coordunitscale = 100;
    handles.ts_defaults.intl.coordsys.unit =...
        handles.gui.unit_length_str;
    %GUI text
    set( handles.ui_cds_unit_text, 'String',...
        sprintf('[ units: %s ]', handles.gui.unit_length_str) );
    %Variable names
    handles.field_props.valu_names = strcat(...
        handles.field_props.valu_name_roots,...
        handles.field_props.valu_units_cm);
    %Variable group names
    handles.field_props.valu_groups = strcat(...
        handles.field_props.valu_groups_roots,...
        handles.field_props.valu_groups_cm);
    %Plot labels:
    handles.plots.axs_label_x = strcat(...
        handles.plots.axs_label_x_roots,...
        handles.plots.axs_label_x_units_cm);
    handles.plots.axs_label_y = strcat(...
        handles.plots.axs_label_y_roots,...
        handles.plots.axs_label_y_units_cm);
    %Adjust axlims according to units
    % -> no change (they are in cm)
elseif handles.gui.units_id == 1000%== m, m/s
    handles.gui.unit_length_str = 'm';
    handles.gui.coordunitscale = 1;
    handles.ts_defaults.intl.coordsys.unit =...
        handles.gui.unit_length_str;
    %GUI text
    set( handles.ui_cds_unit_text, 'String',...
        sprintf('[ units: %s ]', handles.gui.unit_length_str) );
    %Variable names
    handles.field_props.valu_names = strcat(...
        handles.field_props.valu_name_roots,...
        handles.field_props.valu_units_m);
    %Variable group names
    handles.field_props.valu_groups = strcat(...
        handles.field_props.valu_groups_roots,...
        handles.field_props.valu_groups_m);
    %Plot labels:
    handles.plots.axs_label_x = strcat(...
        handles.plots.axs_label_x_roots,...
        handles.plots.axs_label_x_units_m);
    handles.plots.axs_label_y = strcat(...
        handles.plots.axs_label_y_roots,...
        handles.plots.axs_label_y_units_m);
    %Adjust axlims according to units
    % -> no change (they are in cm)
    for p_i=1:size( handles.plots.axs_h, 1)%Plot-types
    for s_i=1:size( handles.plots.axs_h, 2)%Subplots:
        %if subplot exists & those unit comtains length (==unitpower>0)!
        if isempty(handles.plots.axlims_usrdef{ p_i,s_i })==0 && ...
                isempty(find(handles.plots.axlims_unitpower{ p_i,s_i }))==0
            if p_i==5%spectra
                handles.plots.axlims_usrdef{ p_i,s_i }(2,:) =...
                    handles.plots.axlims_usrdef{ p_i,s_i }(2,:) +...
                    log10(10/handles.gui.units_id);
            else
                scaling = (10/handles.gui.units_id).^...
                    handles.plots.axlims_unitpower{ p_i,s_i };
                handles.plots.axlims_usrdef{ p_i,s_i } =...
                    handles.plots.axlims_usrdef{ p_i,s_i }.*scaling;
            end
        end
    end
    end
end
waitbar( 3/10, wbar)

% -------------------------------------------------------------------------
% GUI: Initialize plot properties
[ handles.plots.siz, handles.plots.figpos, handles.plots.axspos ] =...
    gui_plots_init_sizes( handles.plots.siz, 0, handles );

waitbar( 7/10, wbar)

% -------------------------------------------------------------------------
% GUI: default logbox variable:
handles.gui.fig_log = [];
handles.lbox = wbar;%init with an object handle that is deleted later

% -------------------------------------------------------------------------
% APPLY SWITCHES:
% Check input validitys:
% -------------------------------------------------------------------------
% Switch of filtering by defaults
if isempty( find( strcmp( 'nofilt', varargin )==1, 1))~=1
    handles.ts_defaults.p.errfilt.used_types    = {};%empty == no filter
    handles.ts_defaults.p.errfilt.typ = [];
end

waitbar( 10/10, wbar)
% -------------------------------------------------------------------------
handles.plotlimbox = [];
handles.log{end+1,1} = '';
handles.log{end+1,1} = sprintf('~~~ New VIPER session ~~~ ');
% Update handles structure
guidata(hObject, handles);
close(wbar);
% UIWAIT makes viper wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)



% --- Outputs from this function are returned to the command line.
function varargout = viper_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% Get default command line output from handles structure
varargout{1} = handles.output;
% -------------------------------------------------------------------------
%Final gui positioning
movegui( hObject, 'northeast');
fig1pos = get( handles.figure1, 'Position');
if fig1pos(4) > handles.gui.wh_max(2)
    sure = questdlg(...
        'GUI of viper is too large for this screen!',...
        'Screen too small!','Ignore','Close vetistool','Close vetistool');
    if strcmp(sure,'Close vetistool')==1
        figure1_CloseRequestFcn(hObject, eventdata, handles);
    end
end
% movegui( hObject, [ scnsiz(3)-fig1pos(3) scnsiz(4)-fig1pos(4)-30 ])



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isfield( handles, 'dvl')==1
    [ handles, ~ ] = viper_dvl('dvl_close_all', handles, {[]} );
end
if get(handles.ui_plot_fig1_checkbox,       'Value')==1
    guiexec_plots_close( [], [], handles.figure1, [], 1)
end
if get(handles.ui_plot_fig2_checkbox,        'Value')==1
    guiexec_plots_close( [], [], handles.figure1, [], 2)
end
if get(handles.ui_plot_fig3_checkbox,     'Value')==1
    guiexec_plots_close( [], [], handles.figure1, [], 3)
end
if get(handles.ui_plot_fig4_checkbox,        'Value')==1
    guiexec_plots_close( [], [], handles.figure1, [], 4)
end
if get(handles.ui_plot_fig5_checkbox,        'Value')==1
    guiexec_plots_close( [], [], handles.figure1, [], 5)
end
if get(handles.ui_plot_fig6_checkbox,        'Value')==1
    guiexec_plots_close( [], [], handles.figure1, [], 6)
end
if get(handles.ui_plot_fig9_checkbox,        'Value')==1
    guiexec_plots_close( [], [], handles.figure1, [], 9)
end
%If common limit box opern:
if isempty( handles.plotlimbox )~=1
    delete( handles.plotlimbox.fig )
end
%If log figures still open
if isempty( handles.gui.fig_log )~=1
    delete( handles.gui.fig_log )
end
% Hint: delete(hObject) closes the figure
delete(hObject);

function fcn_logbox_closerequest( src, eventdata, h_main_fig)
% Get handles of main GUI:
handles = guidata( h_main_fig );
% Changes:
set(handles.ui__main_logbox_togglebutton, 'Value', 0 )
delete( handles.gui.fig_log )
handles.gui.fig_log = [];
% Update handles structure
guidata(h_main_fig, handles);






% ------------------------------------------------------------------------
% Some recommendations:
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% FUNCTION Nomenclature:
% fcn_... Functions not reading and not writing to handles
% gui_... Functions     reading but not writing to handles
% guiexec_... Functions reading and writing to handles and saving handles
% ------------------------------------------------------------------------
% UPDATE ORDER:
% 1: Update  updates
% gui_update_ui_values(hObject, handles)
% gui_update_ui_enable(hObject, handles) %uses values!
% gui_update_figs(hObject, handles, 0 );
% 2: Update messagebox:
% if isvalid(handles.lbox)==1
% set(handles.lbox,'String',handles.log, 'Value',length(handles.log));
% end
























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 
function chapter_01_main_gui_callbacks
% 
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










function ui__main_addfile_pushbutton_Callback(hObject, eventdata, handles)
if isempty( handles.ts )==1%if first file in session
    curr_dir = pwd;
else
    %Get selected file index
    f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
    curr_dir = handles.fpath{ f_idx };
end
%Waitbar
wbar = waitbar(0.01,'Reading files, please wait...',...
    'WindowStyle','modal');
% ------------------------------------------------------------
% Select files:
% ------------------------------------------------------------
[Filenames,Pathname] = uigetfile('*.dat','Select time-series file(s)',...
    curr_dir,'MultiSelect','on');
    if isequal(Filenames,0)%for ex. in case Cancel was pushed
        %Closing waitbar:
        close(wbar);
        errordlg('No file selected!','!!! ERROR !!!');
        return
    elseif iscell(Filenames) == 0%In case 1 file was selected
        Filenames = cellstr(Filenames);%transform to cell array
    end
    test_ascii = false;
    try
        temp_data = dlmread( [Pathname,Filenames{1} ],'',[0 0 5 4]);
        %temp_data = load( [Pathname,Filenames{1} ], '-ascii' );
        test_ascii = true;
    end
%If first file -> initialize gui
if isempty( handles.ts )==1
    handles.log{end+1,1} =...
        sprintf('    Loading files from:');
    handles.log{end+1,1} =...
        sprintf('    %s',curr_dir );
    handles.log{end+1,1} =...
        sprintf('~~~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    if isvalid( handles.lbox )==1
        set(handles.lbox,'String',handles.log,'Value',length(handles.log));
    end
end
% ------------------------------------------------------------
% Initialize variables:
% ------------------------------------------------------------
num_selected = length(Filenames);%Number of files to add:
num_loaded = 0; %Files that has to be loaded (no error | not yet loaded)
mark_err_loading    = false(num_selected,1); %Lists files that were not possible to load(format)
mark_err_errfilting = false(num_selected,1); %Lists files that were not possible to errfilt
mark_parsloaded     = false(num_selected,1); %Lists files those par-files were loaded
mark_wrn_atload     = false(num_selected,1); %Lists files that gave error or warning at loading
mark_wrn_parsupport = false(num_selected,1); %Lists files with unsupported pars in par-file
mark_wrn_enucds     = false(num_selected,1); %Lists files that have ENU coordinates system
handles.log{end+1,1} =...
    sprintf('>>> Importing Velocity time series file(s)...');
if isvalid( handles.lbox )==1
    set(handles.lbox,'String',handles.log,'Value',length(handles.log));
end

% ------------------------------------------------------------
% CHECK AND LOAD DATA from files one by one:
% ------------------------------------------------------------
curr_log  = {};
if test_ascii==1 || length(Filenames)>1 || isfield( handles, 'dvl')~=1
raw_dsets   = struct([]);
num_steps = 2*num_selected;
for i=1:num_selected
    % Update log here - curr_log hasnt been stored in loop after "continue"
    if isempty( curr_log )==0
        handles.log = vertcat( handles.log, curr_log );
        curr_log = {};
    end
    if isvalid( handles.lbox )==1
        set(handles.lbox,...
            'String',handles.log,'Value',length(handles.log));
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If file was already loaded in session? -> Skip file!
    % else (File is new in the session) -> load available data
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(find(strcmp([Pathname,Filenames{i}],handles.fullname),1))==0
        curr_log{end+1,1} =...
            sprintf('!!! Skipped file %s: (already loaded)',...
            Filenames{i});
        curr_log{end+1,1} =...
            sprintf('    * parameters NOT reloaded');
    else
        % Gather data:
        [ temp_rawdset, temp_log, mark_err_loading(i), mark_parsloaded(i),...
            mark_wrn_atload(i),mark_wrn_parsupport(i),mark_wrn_enucds(i)]=...
            gui_addfile_check_load_files( Pathname, Filenames{i}, handles);
        %Update log:
        curr_log = vertcat( curr_log, temp_log );
        if mark_err_loading(i)==1
            continue
        else
            if isempty( raw_dsets )==1
                raw_dsets = temp_rawdset;
            else
                raw_dsets( end+1,1 ) = temp_rawdset;
            end
        end
        % Store as defaults:
        handles.viper_props.userdef.dat_props.column_t =...
            temp_rawdset.dat_props.column_t;
        handles.viper_props.userdef.dat_props.column_u =...
            temp_rawdset.dat_props.column_u;
    end
    waitbar(i/num_steps)
end
% Display message(s)
if isempty( curr_log )==0
    handles.log = vertcat( handles.log, curr_log );
    curr_log = {};
end
if isvalid( handles.lbox )==1
    set(handles.lbox,...
        'String',handles.log,'Value',length(handles.log));
end
else
    [ ~, raw_dsets ] = viper_dvl('dvl_loadfile',...
        handles, { [ Pathname,Filenames{1} ], handles.fullname });
end
% ------------------------------------------------------------
% IMPORT (== ADD) DATA as new ts-record (file contents are already read)
% At this point it is sure that loaded dataset are supported and get added
% ------------------------------------------------------------
num_steps = num_selected + length( raw_dsets );
for i=1:length( raw_dsets )
    % Update log here - curr_log hasnt been stored in loop after "continue"
    if isempty( curr_log )==0
        handles.log = vertcat( handles.log, curr_log );
        curr_log = {};
    end
    if isvalid( handles.lbox )==1
        set(handles.lbox,...
            'String',handles.log,'Value',length(handles.log));
    end
    %Import raw
    [ temp_ts, temp_fpath, temp_fname, temp_fullname, temp_dispname,...
        temp_log, mark_parsloaded(i), mark_err_errfilting(i) ] =...
        gui_addfile_import_ts_data(...
        raw_dsets(i), mark_parsloaded(i), handles );
    %Get next file index
    f_idx = length(handles.ts) + 1;%(Numfiles already loaded) + 1:
    %Store imported data
    if isempty( handles.ts )==1
        handles.ts = temp_ts;
    else
        handles.ts( f_idx ) = temp_ts;
    end
    %Store other data:
    handles.fpath    = vertcat( handles.fpath, temp_fpath);
    handles.fname    = vertcat( handles.fname, temp_fname);
    handles.fullname = vertcat( handles.fullname, temp_fullname);
    handles.dispname = vertcat( handles.dispname, temp_dispname);
    %Update log:
    curr_log   = vertcat( curr_log, temp_log );
    num_loaded = num_loaded+1;
    waitbar((num_selected+i)/num_steps)
end
% Display message(s)
if isempty( curr_log )==0
    handles.log = vertcat( handles.log, curr_log );
    curr_log = {};
end
if isvalid( handles.lbox )==1
    set(handles.lbox,...
        'String',handles.log,'Value',length(handles.log));
end
% ------------------------------------------------------------
% Save data:
% ------------------------------------------------------------
% Update handles structure
guidata(hObject, handles);
% ------------------------------------------------------------
% Closing:
% ------------------------------------------------------------
%%%% If some files could not be loaded due error (data not recognized)
if sum( mark_err_loading )~=0
    warndlg(['Following files could not be loaded (data not recognized):',...
        strjoin( Filenames( mark_err_loading ) ,', ')],...
        'Files not loaded!');
end
%%%% If some files err.filt_par were not valid
if sum( mark_err_errfilting )~=0
    warndlg(['Following files were not error filtered (invalid pars):',...
        strjoin( Filenames( mark_err_errfilting ) ,', '),...
        'For details please check log window'],...
        'Error filtering failed!');
end
%%%% If parameters of loaded files were loaded
if sum( mark_parsloaded )~=0
    warndlg(['Saved parameters were loaded for at least one time-series ',...
        'file. For concerned file(s) please check log window'],...
        'Parameters loaded!');
end
%%%% If some files loaded with warnings
if sum( mark_wrn_atload )~=0
    warndlg(['At least one time-series file was loaded with warnings. ',...
        'For concerned file(s) please check log window'],...
        'Warnings at loading!');
end
%%%% If some files loaded with warnings
if sum( mark_wrn_parsupport )~=0
    warndlg(['Par-file of at least one time-series contained unsupported',...
        ' parameters, which are ignored in this VIPER session',...
        'For concerned file(s) please check log window'],...
        'Unsupported parameters!');
end
%%%% If some files err.filt_par were not valid
if sum( mark_wrn_enucds )~=0
    warndlg(['Vectors of at least one time-series file was loaded in ENU',...
        ' coordinates! ENU coordinates will be handled as XYZ!',...
        'For concerned file(s) please check log window'],...
        'Coordinate system warning!');
end
%%%% If new files were added:
if num_loaded>0%length(handles.ts)~=0
    % Sorting in alphabetical order:
    [sorted_filenames,sorted_index] = sortrows(handles.fname);
    handles.sorted_index = sorted_index;
    %Last loaded file highlighted:
    lb_index = find( sorted_index== f_idx  );
    set(handles.ui_listbox1,...
        'String',handles.dispname(handles.sorted_index),...
        'Value',lb_index,'Max', f_idx );
    % Calc. new plot characteristics (axlims_minmax):
    [ handles.plots.data_info, handles.plots.axlims_minmax ] =...
        gui_plots_calc_axlims_minmax( handles);
    % Update file selection
    handles.log{end+1,1} =...
        sprintf('--> File selection: %s',...
        handles.fname{ f_idx } );
    % Update GUI:
    gui_update_ui_values(hObject, handles)
    gui_update_ui_enables(hObject, handles)
    gui_plots_update(hObject, handles, 0 );
end
% -----------------------------------------------------------------
% Update handles structure
guidata(hObject, handles);
%Update Log:
if isvalid( handles.lbox )==1
    set(handles.lbox,'String',handles.log, 'Value',length(handles.log));
end
%Closing waitbar:
close(wbar);





% ------------------------------------------------------------
% Log window on|off
% ------------------------------------------------------------
function ui__main_logbox_togglebutton_Callback(hObject, eventdata, handles)
% ------------------------------------------------------------
% If switching on:
% ------------------------------------------------------------
if get(handles.ui__main_logbox_togglebutton, 'Value')==1
    fig1pos = get( handles.figure1, 'Position');
    fig_H = 200;
    fig_H_brut = fig_H + handles.gui.siz_ontop;
    fig_W = fig1pos(3);
    %If enough space below:
    if handles.gui.wh_max(2)-fig1pos(4) > fig_H_brut
        pos_bo = handles.gui.lb_min(2);
        pos_le = fig1pos(1);
    else
        pos_bo = handles.gui.lb_min(2)+...
            handles.gui.wh_max(2) - fig_H;
        pos_le = handles.gui.lb_min(2);
    end
    handles.gui.fig_log = figure('Units','pixels',...
        'PaperPositionMode','auto',...
        'Color',[0.941 0.941 0.941],'MenuBar','none','Toolbar','none',...
        'HandleVisibility','callback','NumberTitle','Off',...
        'Name','Log window',...
        'Position', [ pos_le pos_bo fig_W fig_H]);
    set(handles.gui.fig_log,'CloseRequestFcn',...
        {@fcn_logbox_closerequest,handles.figure1});
    %logbox:
    handles.lbox = uicontrol('Style','listbox',...
        'Parent',handles.gui.fig_log,'BackgroundColor','w',...
        'Units','pixels',...
        'Position',[5 5 fig_W-10 fig_H-10]);
    if isvalid( handles.lbox )==1
        set(handles.lbox, 'FontName', 'FixedWidth')
    end
    %New log and display
    handles.log{end+1,1} =...
        sprintf('~~~ Switching on logbox ~~~ ');
    % Update handles structure
    guidata(hObject, handles);
    if isvalid( handles.lbox )==1
        set(handles.lbox,'String',handles.log,'Value',length(handles.log));
    end
else
    %If log figures still open
    if isempty( handles.gui.fig_log )~=1
        delete( handles.gui.fig_log )
    end
    handles.gui.fig_log = [];
    handles.lbox = [];
end





% ------------------------------------------------------------
% Help main
% ------------------------------------------------------------
function ui__main_help_pushbutton_Callback(hObject, eventdata, handles)
msgbox({...
'1. Load velocity time series file(s) with "Add file". Possible formats:',...
' > .dat file produced by Vectrino/Vector software - time column ',...
'    not required if Vectrino/Vector header file (.hdr) is available',...
' > .dat-file containing an arbitrary velocity time series with ',...
'   - a time column with equidistant discrete times (unit: seconds) and',...
'   - 3 adjoint velocity columns for u, v, w (units: m/s)',...
'',...
'2. During the loading process:',...
' > VIPER automatically gathers parameters:',...
'   - if available, from the VIPER parameter file (.par)',...
'   - else, if available, from the Vectrino/Vector header file (.hdr)',...
'   - otherwise from user input',...
' > VIPER automatically applies error filtering to ADV time series:',...
'   - if available, using paramters loaded from the VIPER .par file',...
'   - otherwise using defaults parameters **',...
'   ** You can modify default error-filtering parameters using "Config."',...
'',...
'3. After loading velocity time series you can:',...
' > *Inspect* data quality using diplayed data, plots and tables',...
' > Check and adjust *Error filtering* applied to time-series',...
' > *Post process* data by setting coordinate system parameters',...
'   - read coordinates from filenames and set measurement point offsets',...
'   - rotate vector direction by 90 degrees around z axis',...
'',...
'4. You can export :',...
'   - error-filtering & coordinate parameters of each file as ascii file',...
'   - statistical results of each time series as .xls or .mat',...
'   - plots as image files (Plots can be configured using "plot config")',...
'   - filtered time series as mat-file or as ascii-file',...
},...
    'Brief instructions...','help','modal')


% ------------------------------------------------------------
% File selection
% ------------------------------------------------------------
function ui_listbox1_Callback(hObject, eventdata, handles)
%Get selected file index
f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
if length( f_idx )==1
    handles.log{end+1,1} =...
    	sprintf('--> File selection: %s',...
        handles.fname{ f_idx } );
else
    handles.log{end+1,1} =...
    	sprintf('--> File selection: multiple files ...');
end
% Update handles structure
guidata(hObject, handles);
% Update GUI:
gui_update_ui_values(hObject, handles)
gui_update_ui_enables(hObject, handles)
gui_plots_update(hObject, handles, 0 );
%Update Log:
if isvalid( handles.lbox )==1
    set(handles.lbox,'String',handles.log, 'Value',length(handles.log));
end










% ------------------------------------------------------------
% Plot selection checkbox
% ------------------------------------------------------------
function ui_plot_fig1_checkbox_Callback(hObject, eventdata, handles)
wbar = waitbar(0.01,'Creating or closing figure, please wait...',...
    'WindowStyle','modal');
if get(handles.ui_plot_fig1_checkbox,       'Value')==1
    [ handles.plots ] = gui_plots_create_figs(handles.figure1,...
        handles, handles.plots, 1, []  );
    waitbar( 0.5, wbar)
    % Update handles structure
    guidata(hObject, handles);
    % Fill plot:
    gui_plots_update(hObject, handles, 1 );
else
    guiexec_plots_close( [], [], handles.figure1, [], 1)
end
close(wbar);

% ------------------------------------------------------------
% Plot selection checkbox
% ------------------------------------------------------------
function ui_plot_fig2_checkbox_Callback(hObject, eventdata, handles)
wbar = waitbar(0.01,'Creating or closing figure, please wait...',...
    'WindowStyle','modal');
if get(handles.ui_plot_fig2_checkbox,        'Value')==1
%     warndlg(...
%         'Noise level is measured only at the beginning of the sampling',...
%         'Warning','modal')
    [ handles.plots ] = gui_plots_create_figs(handles.figure1,...
        handles, handles.plots, 2, []  );
    waitbar( 0.5, wbar)
    % Update handles structure
    guidata(hObject, handles);
    % Fill plot:
    gui_plots_update(hObject, handles, 2 );
else
    guiexec_plots_close( [], [], handles.figure1, [], 2)
end
close(wbar);

% ------------------------------------------------------------
% Plot selection checkbox
% ------------------------------------------------------------
function ui_plot_fig3_checkbox_Callback(hObject, eventdata, handles)
wbar = waitbar(0.01,'Creating or closing figure, please wait...',...
    'WindowStyle','modal');
if get(handles.ui_plot_fig3_checkbox,     'Value')==1
    [ handles.plots ] = gui_plots_create_figs(handles.figure1,...
        handles, handles.plots, 3, []  );
    waitbar( 0.5, wbar)
    % Update handles structure
    guidata(hObject, handles);
    % Fill plot:
    gui_plots_update(hObject, handles, 3 );
else
    guiexec_plots_close( [], [], handles.figure1, [], 3)
end
close(wbar);

% ------------------------------------------------------------
% Plot selection checkbox
% ------------------------------------------------------------
function ui_plot_fig4_checkbox_Callback(hObject, eventdata, handles)
wbar = waitbar(0.01,'Creating or closing figure, please wait...',...
    'WindowStyle','modal');
if get(handles.ui_plot_fig4_checkbox,        'Value')==1
    [ handles.plots ] = gui_plots_create_figs(handles.figure1,...
        handles, handles.plots, 4, []  );
    waitbar( 0.5, wbar)
    % Update handles structure
    guidata(hObject, handles);
    % Fill plot:
    gui_plots_update(hObject, handles, 4 );
else
    guiexec_plots_close( [], [], handles.figure1, [], 4)
end
close(wbar);

% ------------------------------------------------------------
% Plot selection checkbox
% ------------------------------------------------------------
function ui_plot_fig5_checkbox_Callback(hObject, eventdata, handles)
wbar = waitbar(0.01,'Creating or closing figure, please wait...',...
    'WindowStyle','modal');
if get(handles.ui_plot_fig5_checkbox,        'Value')==1
    [ handles.plots ] = gui_plots_create_figs(handles.figure1,...
        handles, handles.plots, 5, []  );
    waitbar( 0.5, wbar)
    % Update handles structure
    guidata(hObject, handles);
    % Fill plot:
    gui_plots_update(hObject, handles, 5 );
else
    guiexec_plots_close( [], [], handles.figure1, [], 5)
end
close(wbar);

% ------------------------------------------------------------
% Plot selection checkbox
% ------------------------------------------------------------
function ui_plot_fig6_checkbox_Callback(hObject, eventdata, handles)
wbar = waitbar(0.01,'Creating or closing figure, please wait...',...
    'WindowStyle','modal');
if get(handles.ui_plot_fig6_checkbox,        'Value')==1
    [ handles.plots ] = gui_plots_create_figs(handles.figure1,...
        handles, handles.plots, 6, []  );
    waitbar( 0.5, wbar)
    % Update handles structure
    guidata(hObject, handles);
    % Fill plot:
    gui_plots_update(hObject, handles, 6 );
else
    guiexec_plots_close( [], [], handles.figure1, [], 6)
end
close(wbar);

% ------------------------------------------------------------
% Plot selection checkbox
% ------------------------------------------------------------
function ui_plot_fig9_checkbox_Callback(hObject, eventdata, handles)
wbar = waitbar(0.01,'Creating or closing figure, please wait...',...
    'WindowStyle','modal');
if get(handles.ui_plot_fig9_checkbox,        'Value')==1
    [ handles.plots ] = gui_plots_create_figs(handles.figure1,...
        handles, handles.plots, 9, [] );
    waitbar( 0.5, wbar)
    % Update handles structure
    guidata(hObject, handles);
    % Fill plot:
    gui_plots_update(hObject, handles, 9 );
else
    guiexec_plots_close( [], [], handles.figure1, [], 9)
end
close(wbar);

% ------------------------------------------------------------
% Plots: Add -5/3 line to loglog axes
% ------------------------------------------------------------
function ui_plot_add53tofig_togglebutton_Callback(hObject, eventdata,handles)
if get( handles.ui_plot_add53tofig_togglebutton ,'Value' )==1
    % ------------------------------------------------------------
    % Questions
    % ------------------------------------------------------------
    sure = questdlg(...
        'This draws a line of power -5/3 onto spectrum plots!',...
        'Sure?',...
        'Draw to default position',...
        'Select position to draw to','Cancel',...
        'Draw to default position');
    %Case Cancel
    if strcmp(sure,'Cancel')==1 || strcmp(sure,'')
        set( handles.ui_plot_add53tofig_togglebutton ,'Value', 0)
        return
    end
    % ------------------------------------------------------------
    % New user selected position - save to handles
    % ------------------------------------------------------------
    if strcmp(sure,'Select position to draw to')==1
        %User selection only possible on spectrum plot
        if isnan( handles.plots.fig_h(5) )==1
            errordlg('Spectrum plot has to exist to select position!',...
                'ERROR!','modal');
            set( handles.ui_plot_add53tofig_togglebutton ,'Value', 0)
            return
        end
        h = msgbox({...
            'The new position will be saved as default position!';...
            '1. Click OK, then';...
            '2. Select the approximate position of the line by mouse'},...
            'Information','help');
        waitfor(h)
        axes( handles.plots.axs_h(5,1) )
        hold on
        [xsel,ysel,c]=ginput(2);
        %Calculate -5/3 values for selected_x and their exponent:
        ycal=xsel.^(-5/3);%ycal
        ycal_exponent=log10(ycal);%
        %Calculate exponent of selected_y
        ysel_exponent=log10(ysel);
        %The nedded exponent to add:
        add_kit=ysel_exponent(1)-ycal_exponent(1);
        %The new exponent and new y values to show:
        yshw_exponent=ycal_exponent+add_kit;
        yshw=10.^yshw_exponent;
        % Save to handles:
        handles.viper_props.userdef.plot_line53_x = xsel;
        handles.viper_props.userdef.plot_line53_y = yshw;
    end
    handles.viper_props.userdef.plot_line53 = 1;
    % Update handles structure
    guidata(hObject, handles);
    % ------------------------------------------------------------
    % Draw to spectrum - if spectrum exist
    % ------------------------------------------------------------
    if isnan( handles.plots.fig_h(5) )~=1
        gui_plots_update(hObject, handles, 5 );
    end
else
    handles.viper_props.userdef.plot_line53 = 0;
    % Update handles structure
    guidata(hObject, handles);
    gui_plots_update(hObject, handles, 5 );
end

% ------------------------------------------------------------
% Plots: Add -3 line to loglog axes
% ------------------------------------------------------------
function ui_plot_add3tofig_togglebutton_Callback(hObject, eventdata, handles)
if get( handles.ui_plot_add3tofig_togglebutton ,'Value' )==1
    % ------------------------------------------------------------
    % Questions
    % ------------------------------------------------------------
    sure = questdlg(...
        'This draws a line of power -3 onto spectrum plots!',...
        'Sure?',...
        'Draw to default position',...
        'Select position to draw to','Cancel',...
        'Draw to default position');
    %Case Cancel
    if strcmp(sure,'Cancel')==1 || strcmp(sure,'')
        set( handles.ui_plot_add3tofig_togglebutton ,'Value', 0)
        return
    end
    % ------------------------------------------------------------
    % New user selected position - save to handles
    % ------------------------------------------------------------
    if strcmp(sure,'Select position to draw to')==1
        %User selection only possible on spectrum plot
        if isnan( handles.plots.fig_h(5) )==1
            errordlg('Spectrum plot has to exist to select position!',...
                'ERROR!','modal');
            set( handles.ui_plot_add3tofig_togglebutton ,'Value', 0)
            return
        end
        h = msgbox({...
            'The new position will be saved as default position!';...
            '1. Click OK, then';...
            '2. Select the approximate position of the line by mouse'},...
            'Information','help');
        waitfor(h)
        axes( handles.plots.axs_h(5,1) )
        hold on
        [xsel,ysel,c]=ginput(2);
        %Calculate -3 values for selected_x and their exponent:
        ycal=xsel.^(-3);%ycal
        ycal_exponent=log10(ycal);%
        %Calculate exponent of selected_y
        ysel_exponent=log10(ysel);
        %The nedded exponent to add:
        add_kit=ysel_exponent(1)-ycal_exponent(1);
        %The new exponent and new y values to show:
        yshw_exponent=ycal_exponent+add_kit;
        yshw=10.^yshw_exponent;
        % Save to handles:
        handles.viper_props.userdef.plot_line3_x = xsel;
        handles.viper_props.userdef.plot_line3_y = yshw;
    end
    handles.viper_props.userdef.plot_line3 = 1;
    % Update handles structure
    guidata(hObject, handles);
    % ------------------------------------------------------------
    % Draw to spectrum - if spectrum exist
    % ------------------------------------------------------------
    if isnan( handles.plots.fig_h(5) )~=1
        gui_plots_update(hObject, handles, 5 );
    end
else
    handles.viper_props.userdef.plot_line3 = 0;
    % Update handles structure
    guidata(hObject, handles);
    gui_plots_update(hObject, handles, 5 );
end


% ------------------------------------------------------------
% Plots: Add -1 line to loglog axes
% ------------------------------------------------------------
function ui_plot_add1tofig_togglebutton_Callback(hObject, eventdata, handles)
if get( handles.ui_plot_add1tofig_togglebutton ,'Value' )==1
    % ------------------------------------------------------------
    % Questions
    % ------------------------------------------------------------
    sure = questdlg(...
        'This draws a line of power -1 onto spectrum plots!',...
        'Sure?',...
        'Draw to default position',...
        'Select position to draw to','Cancel',...
        'Draw to default position');
    %Case Cancel
    if strcmp(sure,'Cancel')==1 || strcmp(sure,'')
        set( handles.ui_plot_add1tofig_togglebutton ,'Value', 0)
        return
    end
    % ------------------------------------------------------------
    % New user selected position - save to handles
    % ------------------------------------------------------------
    if strcmp(sure,'Select position to draw to')==1
        %User selection only possible on spectrum plot
        if isnan( handles.plots.fig_h(5) )==1
            errordlg('Spectrum plot has to exist to select position!',...
                'ERROR!','modal');
            set( handles.ui_plot_add1tofig_togglebutton ,'Value', 0)
            return
        end
        h = msgbox({...
            'The new position will be saved as default position!';...
            '1. Click OK, then';...
            '2. Select the approximate position of the line by mouse'},...
            'Information','help');
        waitfor(h)
        axes( handles.plots.axs_h(5,1) )
        hold on
        [xsel,ysel,c]=ginput(2);
        %Calculate -1 values for selected_x and their exponent:
        ycal=xsel.^(-1);%ycal
        ycal_exponent=log10(ycal);%
        %Calculate exponent of selected_y
        ysel_exponent=log10(ysel);
        %The nedded exponent to add:
        add_kit=ysel_exponent(1)-ycal_exponent(1);
        %The new exponent and new y values to show:
        yshw_exponent=ycal_exponent+add_kit;
        yshw=10.^yshw_exponent;
        % Save to handles:
        handles.viper_props.userdef.plot_line1_x = xsel;
        handles.viper_props.userdef.plot_line1_y = yshw;
    end
    handles.viper_props.userdef.plot_line1 = 1;
    % Update handles structure
    guidata(hObject, handles);
    % ------------------------------------------------------------
    % Draw to spectrum - if spectrum exist
    % ------------------------------------------------------------
    if isnan( handles.plots.fig_h(5) )~=1
        gui_plots_update(hObject, handles, 5 );
    end
else
    handles.viper_props.userdef.plot_line1 = 0;
    % Update handles structure
    guidata(hObject, handles);
    gui_plots_update(hObject, handles, 5 );
end

% ------------------------------------------------------------------------
% Plot component u
% ------------------------------------------------------------------------
function ui_plot_comp_u_checkbox_Callback(hObject, eventdata, handles)
if get(handles.ui_plot_comp_u_checkbox,'Value')==0 &&...
        get(handles.ui_plot_comp_v_checkbox,'Value')==0 &&...
        get(handles.ui_plot_comp_w_checkbox,'Value')==0
    errordlg('At least one component has to be checked! ',...
        '!!! ERROR !!!','modal');
    set(handles.ui_plot_comp_u_checkbox,'Value',1)
end
gui_plots_update(hObject, handles, 0 );

% ------------------------------------------------------------------------
% Plot component v
% ------------------------------------------------------------------------
function ui_plot_comp_v_checkbox_Callback(hObject, eventdata, handles)
if get(handles.ui_plot_comp_u_checkbox,'Value')==0 &&...
        get(handles.ui_plot_comp_v_checkbox,'Value')==0 &&...
        get(handles.ui_plot_comp_w_checkbox,'Value')==0
    errordlg('At least one component has to be checked! ',...
        '!!! ERROR !!!','modal');
    set(handles.ui_plot_comp_v_checkbox,'Value',1)
end
gui_plots_update(hObject, handles, 0 );

% ------------------------------------------------------------------------
% Plot component w
% ------------------------------------------------------------------------
function ui_plot_comp_w_checkbox_Callback(hObject, eventdata, handles)
if get(handles.ui_plot_comp_u_checkbox,'Value')==0 &&...
        get(handles.ui_plot_comp_v_checkbox,'Value')==0 &&...
        get(handles.ui_plot_comp_w_checkbox,'Value')==0
    errordlg('At least one component has to be checked! ',...
        '!!! ERROR !!!','modal');
    set(handles.ui_plot_comp_v_checkbox,'Value',1)
end
gui_plots_update(hObject, handles, 0 );










% ------------------------------------------------------------------------
% Plot limit setting tool
% ------------------------------------------------------------------------
function ui_plot_lims_set_togglebutton_Callback(hObject, eventdata,handles)
% ------------------------------------------------------------
% Create comlimbox:
% ------------------------------------------------------------
if get( handles.ui_plot_lims_set_togglebutton ,'Value' )==1
    
num_plots = sum( handles.plots.adjustable_axlims );
set(0,'Units','pixels') ;
f1_pos = get( handles.figure1, 'Position');
f1_out = get( handles.figure1, 'OuterPosition');
f1_brd = f1_out-f1_pos;
% Create GUI: Radiobuttons for plots + Box for plot limits
fig_H = num_plots*20+240+80;
fig_W = 240;
ui_W_full = 220;
handles.plotlimbox.fig = figure(...
    'Units','pixels','PaperPositionMode','auto',...
    'Color',[0.941 0.941 0.941],'MenuBar','none','Toolbar','none',...
    'HandleVisibility','callback','NumberTitle','Off',...
    'Name','Adjusting plot limits',...
    'Position',...
    [f1_out(1)-fig_W+f1_brd(3)/2 f1_pos(2)+f1_pos(4)-fig_H fig_W fig_H]);
set( handles.plotlimbox.fig, 'CloseRequestFcn',...
    {@guiexec_plotlimbox_closerequest,handles.figure1});
curr_bo = fig_H;
%--------------------------------------------------------------------------
%Create plot selection uicontrols:
curr_bo = curr_bo-30;
handles.plotlimbox.p_tx = uicontrol(...
    'Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 curr_bo 220 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String', 'Select plot type to edit:');%,...
handles.plotlimbox.p_rb = nan(size( handles.plots.adjustable_axlims ));
% posval = find( handles.plots.adjustable_axlims );
i_b = 0;
for i_p=1:length( handles.plots.adjustable_axlims )
if handles.plots.adjustable_axlims( i_p )==1
    i_b = i_b + 1;
    curr_bo = curr_bo-20;
    handles.plotlimbox.p_rb(i_p) = uicontrol(...
        'Parent',handles.plotlimbox.fig,...
        'HandleVisibility','callback','Enable','On','Units','pixels',...
        'Position',[ 30 curr_bo ui_W_full-20 21],'Style','radiobutton',...
        'String', handles.plots.names{ i_p, 1 } );%,...
    if i_b==1
        p_sel = i_p;
        set( handles.plotlimbox.p_rb(i_p), 'Value',1);
    end
    set( handles.plotlimbox.p_rb(i_p), 'Callback',...
        {@guiexec_plotlimbox_callbacks,...
        handles.figure1,handles.plotlimbox.p_rb(i_p)});
end
end
%Create axes selection uicontrols:
curr_bo = curr_bo-25;
handles.plotlimbox.a_tx = uicontrol(...
    'Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 curr_bo 220 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String', 'Select axes to edit:');%,...
handles.plotlimbox.a_rb = nan(size( handles.plots.axs_h, 2),1);
% posval = find( handles.plots.adjustable_axlims );
for i_a=1:size( handles.plots.axs_h, 2)
    curr_bo = curr_bo-20;
    handles.plotlimbox.a_rb(i_a) = uicontrol(...
        'Parent',handles.plotlimbox.fig,...
        'HandleVisibility','callback','Enable','On','Units','pixels',...
        'Position',[ 30 curr_bo 180 21],'Style','radiobutton',...
        'String', handles.plots.axs_names{ p_sel, i_a } );%,...
    if isempty( handles.plots.axlims_usrdef{ p_sel,i_a } )==0
        set( handles.plotlimbox.a_rb( i_a ), 'Enable','On');
    else
        set( handles.plotlimbox.a_rb( i_a ), 'Enable','Off');
    end
    if i_a==1
        a_sel = i_a;
        set( handles.plotlimbox.a_rb(i_a), 'Value',1);
    end
    set( handles.plotlimbox.a_rb(i_a), 'Callback',...
        {@guiexec_plotlimbox_callbacks,...
        handles.figure1,handles.plotlimbox.a_rb(i_a)});
end
%--------------------------------------------------------------------------
%Create axeslimits adjustment uicontrols:
%[X1 X2; Y1 Y2] => index positions: 1 3 2 4
curr_bo = curr_bo-30;
handles.plotlimbox.al_tx0 = uicontrol(...
    'Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 curr_bo 220 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'FontWeight','bold',...
    'String', 'Axis limits of selected axes:');%,...
curr_bo = curr_bo-20;
handles.plotlimbox.al_tx(1) = uicontrol('Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 35 curr_bo 40 21], 'Style', 'text',...
    'String', 'X_min:');%,...
handles.plotlimbox.al_ed(1) = uicontrol('Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 80 curr_bo 40 21], 'Style', 'edit',...
    'String', num2str(handles.plots.axlims_usrdef{ p_sel,a_sel }(1,1)) );
set( handles.plotlimbox.al_ed(1), 'Callback',...
    {@guiexec_plotlimbox_callbacks,...
    handles.figure1,handles.plotlimbox.al_ed(1)});
handles.plotlimbox.al_tx(2) = uicontrol('Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 130 curr_bo 40 21], 'Style', 'text',...
    'String', 'X_max:');%,...
handles.plotlimbox.al_ed(2) = uicontrol('Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 180 curr_bo 40 21], 'Style', 'edit',...
    'String', num2str(handles.plots.axlims_usrdef{ p_sel,a_sel }(1,2)) );
set( handles.plotlimbox.al_ed(2), 'Callback',...
    {@guiexec_plotlimbox_callbacks,...
    handles.figure1,handles.plotlimbox.al_ed(2)});
curr_bo = curr_bo-20;
handles.plotlimbox.al_tx(3) = uicontrol('Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 35 curr_bo 40 21], 'Style', 'text',...
    'String', 'Y_min:');%,...
handles.plotlimbox.al_ed(3) = uicontrol('Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 80 curr_bo 40 21], 'Style', 'edit',...
    'String', num2str(handles.plots.axlims_usrdef{ p_sel,a_sel }(2,1)) );
set( handles.plotlimbox.al_ed(3), 'Callback',...
    {@guiexec_plotlimbox_callbacks,...
    handles.figure1,handles.plotlimbox.al_ed(3)});
handles.plotlimbox.al_tx(4) = uicontrol('Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 130 curr_bo 40 21], 'Style', 'text',...
    'String', 'Y_max:');%,...
handles.plotlimbox.al_ed(4) = uicontrol('Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 180 curr_bo 40 21], 'Style', 'edit',...
    'String', num2str(handles.plots.axlims_usrdef{ p_sel,a_sel }(2,2)) );
set( handles.plotlimbox.al_ed(4), 'Callback',...
    {@guiexec_plotlimbox_callbacks,...
    handles.figure1,handles.plotlimbox.al_ed(4)});
curr_bo = curr_bo-30;
handles.plotlimbox.adj_tx = uicontrol('Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 20 curr_bo-5 70 21], 'Style', 'text',...
    'HorizontalAlignment','center',...
    'String', 'Edit or:');%,...
%curr_bo = curr_bo-20;
handles.plotlimbox.collect = uicontrol('Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Enable','Off','Units','pixels',...
    'Position',[ 80 curr_bo 150 21], 'Style', 'pushbutton' ,...
    'String', 'Load min. & max. over all files');%,...
set( handles.plotlimbox.collect, 'Callback',...
    {@guiexec_plotlimbox_callbacks,...
    handles.figure1,handles.plotlimbox.collect});
if isfield(handles, 'sorted_index')==1
    set( handles.plotlimbox.collect, 'Enable','On' );
end
curr_bo = curr_bo-30;
handles.plotlimbox.check = uicontrol('Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Enable','On','Units','pixels',...
    'Position',[ 10 curr_bo 220 21],'Style','checkbox',...
    'Value', handles.plots.axlims_to_use( p_sel,a_sel ),...
    'String', 'Use these axis limits instead of default limits' );%,...
set( handles.plotlimbox.check, 'Callback',...
    {@guiexec_plotlimbox_callbacks,...
    handles.figure1,handles.plotlimbox.check});
curr_bo = curr_bo-30;
handles.plotlimbox.apply = uicontrol('Parent',handles.plotlimbox.fig,...
    'HandleVisibility','callback','Enable','Off','Units','pixels',...
    'Position',[ 10 curr_bo 220 21], 'Style', 'pushbutton' ,...
    'String', 'Apply these settings for selected axes');%,...
set( handles.plotlimbox.apply, 'Callback',...
    {@guiexec_plotlimbox_callbacks,...
    handles.figure1,handles.plotlimbox.apply});
%--------------------------------------------------------------------------
%Store to handles:
guidata(handles.figure1, handles);% Update handles structure

% ------------------------------------------------------------
% Close comlimbox:
% ------------------------------------------------------------
else
    guiexec_plotlimbox_closerequest( [], [], handles.figure1 )
end










% ------------------------------------------------------------------------
% Plot configuration tool
% ------------------------------------------------------------------------
function ui_plot_config_pushbutton_Callback(hObject, eventdata, handles)
% ------------------------------------------------------------
% Create box:
% ------------------------------------------------------------
num_plots = sum( handles.plots.adjustable_axlabels );
set(0,'Units','pixels') ;
f1_pos = get( handles.figure1, 'Position');
f1_out = get( handles.figure1, 'OuterPosition');
f1_brd = f1_out-f1_pos;
% Create GUI:
fig_H = num_plots*20+510;
fig_W = 280;
ui_W_full = 260;
handles.plotsetupbox.fig = figure(...
    'Units','pixels','PaperPositionMode','auto',...
    'Color',[0.941 0.941 0.941],'MenuBar','none','Toolbar','none',...
    'HandleVisibility','callback','NumberTitle','Off',...
    'Name','Plot setup', 'WindowStyle','modal',...
    'Position',...
    [f1_out(1)-fig_W+f1_brd(3)/2 f1_pos(2)+f1_pos(4)-fig_H fig_W fig_H]);
set( handles.plotsetupbox.fig, 'CloseRequestFcn',...
    {@guiexec_plotsetupbox_closerequest,handles.figure1});
curr_bo = fig_H;
%--------------------------------------------------------------------------
%Create size uicontrols:
curr_bo = curr_bo-25;
handles.plotsetupbox.sz_tx = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 curr_bo ui_W_full 21], 'Style', 'text',...
    'HorizontalAlignment','left','FontWeight','bold',...
    'String', 'Settings common for all plots:');
%Figure onscreen width
curr_bo = curr_bo-20;
handles.plotsetupbox.dw_tx = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 curr_bo-7 ui_W_full 31], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String',...
    sprintf('On-screen width of standard plots, %% \n %s',...
        '(percent of screen width )'));
handles.plotsetupbox.dw_ed = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Enable','On','Units','pixels',...
    'Position',[ 220 curr_bo 50 21],'Style','edit',...
    'String', num2str(handles.plots.siz.used_screen_wr*100) );
set( handles.plotsetupbox.dw_ed, 'Callback',...
    {@guiexec_plotsetupbox_callbacks,...
    handles.figure1,handles.plotsetupbox.dw_ed});
%Create Font type uicontrols:
curr_bo = curr_bo-25;
gui_ft_types = {'Times New Roman' 'Arial'};
ft_sel = find(strcmp( handles.plots.ft_name, gui_ft_types));
handles.plotsetupbox.ft_tx = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 curr_bo-5 100 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String', 'Plot font type:');
handles.plotsetupbox.ft_pm = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Enable','On','Units','pixels',...
    'Position',[ 90 curr_bo 180 21],'Style','popupmenu',...
    'String', gui_ft_types,...
    'Value', ft_sel);%,...
set( handles.plotsetupbox.ft_pm, 'Callback',...
    {@guiexec_plotsetupbox_callbacks,...
    handles.figure1,handles.plotsetupbox.ft_pm});
%Figure export width
curr_bo = curr_bo-25;
handles.plotsetupbox.w_tx = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 curr_bo-7 ui_W_full 31], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String',...
    sprintf('In-print image width of standard plots, mm: \n %s',...
        '(spectra is of half width):'));
handles.plotsetupbox.w_ed = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Enable','On','Units','pixels',...
    'Position',[ 220 curr_bo 50 21],'Style','edit',...
    'String', num2str(handles.plots.siz.printwidth_mm(3)) );%3=standard fig
set( handles.plotsetupbox.w_ed, 'Callback',...
    {@guiexec_plotsetupbox_callbacks,...
    handles.figure1,handles.plotsetupbox.w_ed});
%Figure export font size
curr_bo = curr_bo-25;
handles.plotsetupbox.fs_tx = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 curr_bo-5 220 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String', 'In-print font size of plots in image file, pt:');%,...
handles.plotsetupbox.fs_ed = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Enable','On','Units','pixels',...
    'Position',[ 220 curr_bo 50 21],'Style','edit',...
    'String', num2str(handles.plots.siz.ft_pt) );%,...
set( handles.plotsetupbox.fs_ed, 'Callback',...
    {@guiexec_plotsetupbox_callbacks,...
    handles.figure1,handles.plotsetupbox.fs_ed});
%--------------------------------------------------------------------------
%Create apply commons uicontrol:
curr_bo = curr_bo-25;
handles.plotsetupbox.app_co = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Enable','Off','Units','pixels',...
    'Position',[ 10 curr_bo ui_W_full 21], 'Style', 'pushbutton' ,...
    'String', 'Apply for all plots');%,...
set( handles.plotsetupbox.app_co, 'Callback',...
    {@guiexec_plotsetupbox_callbacks,...
    handles.figure1,handles.plotsetupbox.app_co});
%--------------------------------------------------------------------------
%Create plots & axes uicontrols:
curr_bo = curr_bo-40;
handles.plotsetupbox.sp_tx = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 curr_bo ui_W_full 21], 'Style', 'text',...
    'HorizontalAlignment','left','FontWeight','bold',...
    'String', 'Plot specific settings:');
%Create plot selection uicontrols:
curr_bo = curr_bo-25;
handles.plotsetupbox.p_tx = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 curr_bo ui_W_full 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String', 'Select plot type to edit:');%,...
handles.plotsetupbox.p_rb = nan(size( handles.plots.adjustable_axlabels ));
i_b = 0;
for i_p=1:length( handles.plots.adjustable_axlabels )
if handles.plots.adjustable_axlabels( i_p )==1
    i_b = i_b + 1;
    curr_bo = curr_bo-20;
    handles.plotsetupbox.p_rb(i_p) = uicontrol(...
        'Parent',handles.plotsetupbox.fig,...
        'HandleVisibility','callback','Enable','On','Units','pixels',...
        'Position',[ 30 curr_bo 190 21],'Style','radiobutton',...
        'String', handles.plots.names{ i_p, 1 } );%,...
    if i_b==1
        p_sel = i_p;
        set( handles.plotsetupbox.p_rb(i_p), 'Value',1);
    end
    set( handles.plotsetupbox.p_rb(i_p), 'Callback',...
        {@guiexec_plotsetupbox_callbacks,...
        handles.figure1,handles.plotsetupbox.p_rb(i_p)});
end
end
%Create axes selection uicontrols:
curr_bo = curr_bo-25;
handles.plotsetupbox.a_tx = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 curr_bo ui_W_full 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String', 'Select axes to edit:');%,...
handles.plotsetupbox.a_rb = nan(size( handles.plots.axs_h, 2),1);
% posval = find( handles.plots.adjustable_axlabels );
for i_a=1:size( handles.plots.axs_h, 2)
    curr_bo = curr_bo-20;
    handles.plotsetupbox.a_rb(i_a) = uicontrol(...
        'Parent',handles.plotsetupbox.fig,...
        'HandleVisibility','callback','Enable','On','Units','pixels',...
        'Position',[ 30 curr_bo 190 21],'Style','radiobutton',...
        'String', handles.plots.axs_names{ p_sel, i_a } );%,...
    if isempty( handles.plots.axlims_usrdef{ p_sel,i_a } )==0
        set( handles.plotsetupbox.a_rb( i_a ), 'Enable','On');
    else
        set( handles.plotsetupbox.a_rb( i_a ), 'Enable','Off');
    end
    if i_a==1
        a_sel = i_a;
        set( handles.plotsetupbox.a_rb(i_a), 'Value',1);
    end
    set( handles.plotsetupbox.a_rb(i_a), 'Callback',...
        {@guiexec_plotsetupbox_callbacks,...
        handles.figure1,handles.plotsetupbox.a_rb(i_a)});
end
%--------------------------------------------------------------------------
%Create axeslimits adjustment uicontrols:
%[X1 X2; Y1 Y2] => index positions: 1 3 2 4
curr_bo = curr_bo-30;
handles.plotsetupbox.al_tx0 = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 curr_bo ui_W_full 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'FontWeight','bold',...
    'String', 'Axis labels of selected axes:');
curr_bo = curr_bo-25;
handles.plotsetupbox.al_tx(1) = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 30 curr_bo ui_W_full-20 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String', 'X-axis label (without unit):');%,...
curr_bo = curr_bo-15;
handles.plotsetupbox.al_ed(1) = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 30 curr_bo ui_W_full-20 21], 'Style', 'edit',...
    'HorizontalAlignment','left',...
    'String', handles.plots.axs_label_x_roots{ p_sel,a_sel } );%,...
set( handles.plotsetupbox.al_ed(1), 'Callback',...
    {@guiexec_plotsetupbox_callbacks,...
    handles.figure1,handles.plotsetupbox.al_ed(1)});
curr_bo = curr_bo-25;
handles.plotsetupbox.al_tx(2) = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 30 curr_bo ui_W_full-20 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String', 'Y-axis label (without unit):');%,...
curr_bo = curr_bo-15;
handles.plotsetupbox.al_ed(2) = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 30 curr_bo ui_W_full-20 21], 'Style', 'edit',...
    'HorizontalAlignment','left',...
    'String', handles.plots.axs_label_y_roots{ p_sel,a_sel } );%,...
set( handles.plotsetupbox.al_ed(2), 'Callback',...
    {@guiexec_plotsetupbox_callbacks,...
    handles.figure1,handles.plotsetupbox.al_ed(2)});
curr_bo = curr_bo-30;
handles.plotsetupbox.ag_tx0 = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 curr_bo ui_W_full 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'FontWeight','bold',...
    'String', 'Legend of selected plot:');
curr_bo = curr_bo-15;
handles.plotsetupbox.ag_cb = uicontrol(...
    'Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Enable','On','Units','pixels',...
    'Position',[ 30 curr_bo 190 21],'Style','checkbox',...
    'Value',handles.plots.axs_legend_on(p_sel),...
    'String', 'Show legend' );%,...
if handles.plots.adjustable_axlegend(p_sel)==1
    set( handles.plotsetupbox.ag_cb, 'Enable','On');
else
    set( handles.plotsetupbox.ag_cb, 'Enable','Off');
end
set( handles.plotsetupbox.ag_cb, 'Callback',...
    {@guiexec_plotsetupbox_callbacks,...
    handles.figure1,handles.plotsetupbox.ag_cb});
%--------------------------------------------------------------------------
%Create apply uicontrols:
curr_bo = curr_bo-35;
handles.plotsetupbox.app_lb = uicontrol('Parent',handles.plotsetupbox.fig,...
    'HandleVisibility','callback','Enable','Off','Units','pixels',...
    'Position',[ 10 curr_bo ui_W_full 21], 'Style', 'pushbutton' ,...
    'String', 'Apply legend & label setting for selected plot & axes');%,...
set( handles.plotsetupbox.app_lb, 'Callback',...
    {@guiexec_plotsetupbox_callbacks,...
    handles.figure1,handles.plotsetupbox.app_lb});
%--------------------------------------------------------------------------
%Store to handles:
guidata(handles.figure1, handles);% Update handles structure










% ------------------------------------------------------------
% Duration view: variable selection
% ------------------------------------------------------------
function ui_est_variable_popupmenu_Callback(hObject, eventdata, handles)
%Get selected file index
f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
gui_update_estimator_display( f_idx, handles )

% ------------------------------------------------------------
% Duration view: error selection
% ------------------------------------------------------------
function ui_est_edit_Callback(hObject, eventdata, handles)
input = str2double( regexprep( get(hObject,'String'), ',', '.') );
set(handles.ui_est_edit,'String', num2str(input) )
%Get selected file index
f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
gui_update_estimator_display( f_idx, handles )










% ------------------------------------------------------------------------
% Errorfiltering: select component u to use
% ------------------------------------------------------------------------
function ui_errf_cmp_u_checkbox_Callback(hObject,eventdata,handles)
%comp_u comp_v comp_w - at least one has to be checked
if get(handles.ui_errf_cmp_u_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_cmp_v_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_cmp_w_checkbox, 'Value')==0
    warndlg('At least one component has to be checked!',...
        'Component selection!');
    set(handles.ui_errf_cmp_u_checkbox, 'Value', 1)
    return
end
gui_update_ui_enables(hObject, handles)

% ------------------------------------------------------------------------
% Errorfiltering: select component v to use
% ------------------------------------------------------------------------
function ui_errf_cmp_v_checkbox_Callback(hObject,eventdata,handles)
%comp_u comp_v comp_w - at least one has to be checked
if get(handles.ui_errf_cmp_u_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_cmp_v_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_cmp_w_checkbox, 'Value')==0
    warndlg('At least one component has to be checked!',...
        'Component selection!');
    set(handles.ui_errf_cmp_v_checkbox, 'Value', 1)
    return
end
gui_update_ui_enables(hObject, handles)

% ------------------------------------------------------------------------
% Errorfiltering: select component w to use
% ------------------------------------------------------------------------
function ui_errf_cmp_w_checkbox_Callback(hObject,eventdata,handles)
%comp_u comp_v comp_w - at least one has to be checked
if get(handles.ui_errf_cmp_u_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_cmp_v_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_cmp_w_checkbox, 'Value')==0
    warndlg('At least one component has to be checked!',...
        'Component selection!');
    set(handles.ui_errf_cmp_w_checkbox, 'Value', 1)
    return
end
gui_update_ui_enables(hObject, handles)


% ------------------------------------------------------------------------
% Errorfiltering: select crop
% ------------------------------------------------------------------------
function ui_errf_crop_checkbox_Callback(hObject, eventdata, handles)
if isfield(handles, 'sorted_index')==1
    %Get selected file index
    f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
end
if get(handles.ui_errf_crop_checkbox,        'Value')==1
    %If editbox==NaN -> replace by first/last sample to simplify selection
    % (editbox had been set to NaN at change to this file, if checkbox==0)
    if isnan( str2double( regexprep(...
            get(handles.ui_errf_crop_1_edit,'String'), ',', '.') ) )==1
        set(handles.ui_errf_crop_1_edit,'String',...
            num2str( handles.ts(f_idx).raw_t(1) ) )
    end
    if isnan( str2double( regexprep(...
            get(handles.ui_errf_crop_2_edit,'String'), ',', '.') ) )==1
        set(handles.ui_errf_crop_2_edit,'String',...
            num2str( handles.ts(f_idx).raw_t(end) ) )
    end
%else no change, and next file change updates edit
end
gui_update_ui_enables(hObject, handles)

% ------------------------------------------------------------------------
% Errorfiltering: select crop parameter
% ------------------------------------------------------------------------
function ui_errf_crop_1_edit_Callback(hObject, eventdata, handles)
if isfield(handles, 'sorted_index')==1
    %Get selected file index
    f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
end
%Check validity
input = str2double(...
    regexprep( get(handles.ui_errf_crop_1_edit,'String'), ',', '.') );
set(handles.ui_errf_crop_1_edit,'String', num2str(input) )
if isnan(input)==1
    input = handles.ts(f_idx).raw_t(1);
end
if ( isfield(handles, 'sorted_index') && ...
        input < handles.ts(f_idx).raw_t(1)  )
    warndlg(sprintf(...
        'Value lower than discrete time of first time step!'));
%     set(handles.ui_errf_crop_1_edit,'String',...
%         num2str( handles.ts(f_idx).raw_t(1) ) )%first sample
end
if ( isfield(handles, 'sorted_index') && ...
        input >= handles.ts(f_idx).raw_t(end)  )
    errordlg(sprintf(...
        'Value higher than discrete time of last time step!'));
%     set(handles.ui_errf_crop_1_edit,'String',...
%         num2str( handles.ts(f_idx).raw_t(1) ) )%first sample
end
%Check whether 1<2:
% -> fcn_errfiltpars_check_validity

% ------------------------------------------------------------------------
% Errorfiltering: select crop parameter
% ------------------------------------------------------------------------
function ui_errf_crop_2_edit_Callback(hObject, eventdata, handles)
if isfield(handles, 'sorted_index')==1
    %Get selected file index
    f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
end
%Check validity
input = str2double(...
    regexprep( get(handles.ui_errf_crop_2_edit,'String'), ',', '.') );
set(handles.ui_errf_crop_2_edit,'String', num2str(input) )
if isnan(input)==1
    input = handles.ts(f_idx).raw_t(end);
end
if (    isfield(handles, 'sorted_index') && ...
        input <= handles.ts(f_idx).raw_t(1)  )
    errordlg(sprintf(...
        'Value lower than discrete time of first time step'));
%     set(handles.ui_errf_crop_2_edit,'String',...
%         num2str( handles.ts(f_idx).raw_t(end) ) )%last sample
end
if (    isfield(handles, 'sorted_index') && ...
        input > handles.ts(f_idx).raw_t(end)  )
    errordlg(sprintf(...
        'Value higher than discrete time of last time step!'));
%     set(handles.ui_errf_crop_2_edit,'String',...
%         num2str( handles.ts(f_idx).raw_t(end) ) )%last sample
end
%Check whether 1<2
% -> fcn_errfiltpars_check_validity


% ------------------------------------------------------------------------
% Errorfiltering: select PST-filter to use
% ------------------------------------------------------------------------
function ui_errf_pst_checkbox_Callback(hObject, eventdata, handles)
if isfield( handles, 'dvl')==1
    [ ~, ~ ] = viper_dvl('dvl_update_errfilt_selection', handles, {[]} );
end
gui_update_ui_enables(hObject, handles)
if get(handles.ui_errf_pst_cmp_u_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_pst_cmp_v_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_pst_cmp_w_checkbox, 'Value')==0
    warndlg(...
        'Check at least one component, othewise filter get switched off!',...
        'Component selection!');
end

% ------------------------------------------------------------------------
% Errorfiltering: select PST-filter component u to use
% ------------------------------------------------------------------------
function ui_errf_pst_cmp_u_checkbox_Callback(hObject,eventdata,handles)
%comp_u comp_v comp_w - at least one has to be checked
if get(handles.ui_errf_pst_cmp_u_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_pst_cmp_v_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_pst_cmp_w_checkbox, 'Value')==0
    warndlg(...
        'Check at least one component, othewise filter get switched off!',...
        'Component selection!');
end

% ------------------------------------------------------------------------
% Errorfiltering: select PST-filter component v to use
% ------------------------------------------------------------------------
function ui_errf_pst_cmp_v_checkbox_Callback(hObject,eventdata,handles)
%comp_u comp_v comp_w - at least one has to be checked
if get(handles.ui_errf_pst_cmp_u_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_pst_cmp_v_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_pst_cmp_w_checkbox, 'Value')==0
    warndlg(...
        'Check at least one component, othewise filter get switched off!',...
        'Component selection!');
end

% ------------------------------------------------------------------------
% Errorfiltering: select PST-filter component w to use
% ------------------------------------------------------------------------
function ui_errf_pst_cmp_w_checkbox_Callback(hObject,eventdata,handles)
%comp_u comp_v comp_w - at least one has to be checked
if get(handles.ui_errf_pst_cmp_u_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_pst_cmp_v_checkbox, 'Value')==0 &&...
        get(handles.ui_errf_pst_cmp_w_checkbox, 'Value')==0
    warndlg(...
        'Check at least one component, othewise filter get switched off!',...
        'Component selection!');
end

% ------------------------------------------------------------------------
% Errorfiltering: select PST-filter high pass data to use
% ------------------------------------------------------------------------
function ui_errf_pst_usehpf_checkbox_Callback(hObject, eventdata, handles)
% if get(handles.ui_errf_pst_usehpf_checkbox, 'Value')==0
%     set(handles.ui_errf_pst_hpfreq_edit,'String','NaN')
% end%not needed, as gui_collect changes it to NaN
gui_update_ui_enables(hObject, handles)

% ------------------------------------------------------------------------
% Errorfiltering: select PST-filter high pass frequency
% ------------------------------------------------------------------------
function ui_errf_pst_hpfreq_edit_Callback(hObject, eventdata, handles)
if isfield(handles, 'sorted_index')==1
    %Get selected file index
    f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
end
%Check
input = str2double(...
    regexprep( get(handles.ui_errf_pst_hpfreq_edit,'String'), ',', '.') );
set(handles.ui_errf_pst_hpfreq_edit,'String', num2str(input) )
if isnan(input)==1 || input<=0 || ...
        (  isfield(handles, 'sorted_index') && ...
        input>handles.ts( f_idx ).p.dat_props.frq  )
    errordlg(sprintf('Invalid input!'),'Error','modal');
    set(handles.ui_errf_pst_hpfreq_edit,'String','NaN')%<=0 not allowed
end

% ------------------------------------------------------------------------
% Errorfiltering: select COR filter COR-min to use
% ------------------------------------------------------------------------
function ui_errf_cor_min_checkbox_Callback(hObject, eventdata, handles)
if get(handles.ui_errf_cor_min_checkbox,        'Value')==1
    set(handles.ui_errf_cor_avg_checkbox,       'Value', 0);
end

% ------------------------------------------------------------------------
% Errorfiltering: select COR filter COR-avg to use
% ------------------------------------------------------------------------
function ui_errf_cor_avg_checkbox_Callback(hObject, eventdata, handles)
if get(handles.ui_errf_cor_avg_checkbox,        'Value')==1
    set(handles.ui_errf_cor_min_checkbox,       'Value', 0);
end

% ------------------------------------------------------------------------
% Errorfiltering: select COR filter lambda
% ------------------------------------------------------------------------
function ui_errf_cor_lambda_edit_Callback(hObject, eventdata, handles)
input = str2double(...
    regexprep( get(handles.ui_errf_cor_lambda_edit,'String'), ',', '.') );
set(handles.ui_errf_cor_lambda_edit,'String', num2str(input) )
if isnan(input)==1 || input<=0 || input>100
    errordlg(sprintf('Invalid input!'));
    set(handles.ui_errf_cor_lambda_edit,'String',num2str(NaN))
end





% ------------------------------------------------------------------------
% Errorfiltering: reaload from parfile
% ------------------------------------------------------------------------
function ui_errf_reload_par_pushbutton_Callback(hObject,eventdata,handles)
%Get selected file index
f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
Pathname = handles.fpath{ f_idx ,1};
name     = handles.fname{ f_idx ,1};
ask_loading_default_pars = 0;
% ------------------------------------------------------------
% Load .parfile - if exists
if exist([Pathname name '.par'],'file')==2
    pargroups = fieldnames( handles.ts_defaults.p );
    [ parfile_pars, parfile_pargroups, unknown_pars, log_msg_c ] =...
        fcn_parfile_load( [Pathname name '.par'], pargroups, ...
        handles.ts_defaults.p,...
        handles.viper_props.errf_types_supported(:,1));
    % ------------------------------------------------------------
    % If errfilt_pars existed within parfile:
    if sum(strcmpi(parfile_pargroups,'errfilt'))>0
        %unsupported pars:
        if isempty( unknown_pars )~=1 && isfield(unknown_pars,'errfilt')
            warndlg(['There are unsupported error-filtering parameters ',...
                'in par-file, which are ignored.'],...
                'WindowStyle','modal');
        end
        % Check whether errfilt pars are valid:
        [ err_msg ]= fcn_errfiltpars_check_validity(...
            parfile_pars.errfilt, [], isfield( handles, 'dvl') );
        % ------------------------------------------------------------
        % If valid errfilt_pars -> load to handles
        if isempty( err_msg )==1
            handles.ts( f_idx ).p.errfilt = parfile_pars.errfilt;
        % ELSE
        else
            q_msg{1} = 'Invalid error-filtering parameters in parfile!';
            q_msg(2:length( err_msg )+1) = err_msg;
            ask_loading_default_pars = 1;
        end
    % error in parfile or no errfilt_pars within parfile:
    else
        q_msg{1} ='Error loading error-filtering parameters from parfile!';
        q_msg(2:length( log_msg_c )+1) = log_msg_c;
        ask_loading_default_pars = 1;
    end
% No .parfile
else
    q_msg{1} = 'There exists no parfile to this file!';
    ask_loading_default_pars = 1;
end
% ------------------------------------------------------------
% Load defaults:
if ask_loading_default_pars==1
    q_msg{end+1} = '\nLoad default err. filtering parameters instead?';
    usr_input = questdlg(q_msg,'Load defaults?','Yes','No','Yes');
    if strcmp( usr_input,'Yes')~=1
        return
    end
    %Deafult errfilt pars
    handles.ts( f_idx ).p.errfilt = handles.ts_defaults.p.errfilt;
end
% Update GUI:
gui_update_ui_values(hObject, handles)
% Enable uicontrols:
gui_update_ui_enables(hObject, handles)





% ------------------------------------------------------------------------
% Errorfiltering: input edit or cancel
% ------------------------------------------------------------------------
function ui_errf_edit_n_cancel_toggle_b_Callback(hObject,eventdata,handles)
% ------------------------------------------------------------
% If switching EDIT on:
% ------------------------------------------------------------
if get(handles.ui_errf_edit_n_cancel_toggle_b, 'Value')==1
    % Change display name to "Cancel":
    set(handles.ui_errf_edit_n_cancel_toggle_b,'String','Cancel')
% ------------------------------------------------------------
% If cancelling == switching EDIT to off
% ------------------------------------------------------------
else
    % Change display name to "Edit":
    set(handles.ui_errf_edit_n_cancel_toggle_b,'String','Edit')
    % Cancel changes by loading stored parameters:
    gui_update_ui_values(hObject, handles)
end
% Enable uicontrols:
gui_update_ui_enables(hObject, handles)

% ------------------------------------------------------------------------
% Errorfiltering: input apply
% ------------------------------------------------------------------------
function ui_errf_apply_pushbutton_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Depending on which Mode is switched on:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------
% If Editing file Err. filtering switched on:
% ------------------------------------------------------------
if get(handles.ui_errf_edit_n_cancel_toggle_b, 'Value')==1
    %
    wbar = waitbar(0.01,...
        'Checking parameters & determining erroneous samples, please wait...',...
        'WindowStyle','modal');
    %Get selected file index
    f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
    % ------------------------------------------------------------
    % Final check of PST-Filter-GUI Valus: if all comp==0, pst->0
    % ------------------------------------------------------------
    if get(handles.ui_errf_pst_cmp_u_checkbox, 'Value')==0 &&...
            get(handles.ui_errf_pst_cmp_v_checkbox, 'Value')==0 &&...
            get(handles.ui_errf_pst_cmp_w_checkbox, 'Value')==0
        set(handles.ui_errf_pst_usehpf_checkbox, 'Value',0)
        set(handles.ui_errf_pst_checkbox, 'Value',0)
    end
    % ------------------------------------------------------------
    % Collect displayed err.filt parameters:
    % ------------------------------------------------------------
    [ disped_p_errfilt ] = gui_errfiltpars_collect( handles );
    % ------------------------------------------------------------
    % Check errfilt pars for validity == applicability to selected file:
    % ------------------------------------------------------------
    [ err_msg ] = fcn_errfiltpars_check_validity(...
        disped_p_errfilt,...
        handles.ts( f_idx ).p.dat_props, isfield(handles,'dvl') );
    if isempty( err_msg )~=1
        errordlg( err_msg, 'Error in selected paramters');
        close(wbar);
        return
    end
    % ------------------------------------------------------------
    % Collect stored parameters:
    % ------------------------------------------------------------
    stored_p_errfilt = handles.ts( f_idx ).p.errfilt;
    % ------------------------------------------------------------
    % Comparing stored and displayed - if mismatch => determine using new
    % ------------------------------------------------------------
    if isequalwithequalnans( stored_p_errfilt , disped_p_errfilt)~=1
        % ------------------------------------------------------------
        % Do error_filtering
        % ------------------------------------------------------------
        [ efilt_res, err_msg ] = fcn_errfilt__main(...
            0, disped_p_errfilt,...
            handles.ts( f_idx ).p.dat_props,...
            handles.ts( f_idx ).raw_veldata,...
            handles.ts( f_idx ).raw_vals,...
            handles.ts( f_idx ).raw_t,...
            handles.ts( f_idx ).raw_corsnr_xt,...
            isfield( handles, 'dvl'));
        if isempty( err_msg )~=1
            errordlg( err_msg, 'Error during error filtering');
            close(wbar);
            return
        end
        %Store new errfilt pars and efilt_res:
        handles.ts( f_idx ).p.errfilt = disped_p_errfilt;
        handles.ts( f_idx ).efilt_res = efilt_res;
        %Enabling elements: (everything off)% = set to 'Off'
        handles.ts( f_idx ).seri_enable = handles.ts_defaults.seri_enable;
        % -------------------------------------------------------------
        % Store NEW vals due to errorfiltering:
        % ------------------------------------------------------------
        handles.ts( f_idx ).vals =...
            handles.ts( f_idx ).efilt_res.valmarker(:,1);
        %Create & display message:
        handles.log{end+1,1} =...
            sprintf('  > Set NEW err.filt parameters & sample validity');
        if isvalid( handles.lbox )==1
            set(handles.lbox,'String', handles.log,'Value',length(handles.log));
        end
        % ------------------------------------------------------------
        % Calculate COR & SNR distribution for valid data
        % ------------------------------------------------------------
        if strcmp( handles.ts( f_idx ).p.dat_props.type, 'vno') ||...
                strcmp( handles.ts( f_idx ).p.dat_props.type, 'vec')
            %Calculate quality-plots:
            [ temp_corsnr, temp_boxp ] = fcn_calc_corsnrdist_vals_adv_sokoray(...
                handles.ts( f_idx ).vals,...
                handles.ts( f_idx ).raw_corsnr_xt );
            %Store results:
            %Correlation distribution+cumulative & SNR distribution+cumulative:
            handles.ts( f_idx ).corsnr_dist = temp_corsnr.corsnr_dist;
            %Boxplot data:
            handles.ts( f_idx ).valu_data( 19:23 ) = temp_boxp(1:5);
            %Enable
            handles.ts( f_idx ).valu_enable( 19:23 )=1;%COR-distrib boxplot
        end
        % ------------------------------------------------------------
        % Calculate velocity histograms
        % ------------------------------------------------------------
        [ handles.ts( f_idx ).hist.p_per_bin,...
            handles.ts( f_idx ).hist.vel_of_bin ] =...
            fcn_calc_histograms_of_vals(...
            handles.ts( f_idx ).raw_veldata,...
            handles.ts( f_idx ).vals, ...
            handles.ts( f_idx ).p.errfilt.used_velcomp );
        % -----------------------------------------------------------------
        % Calc. err.filtered fluctuation & basic statistics using .vals:
        % -----------------------------------------------------------------
        [ handles.ts( f_idx ).seri_data,...
            handles.ts( f_idx ).seri_enable,...
            handles.ts( f_idx ).intl.mean_correctn,...
            handles.ts( f_idx ).intl.used_ts_seg,...
            new_f_seri_x,...
            new_f_valu_data,...
            new_f_valu_enable,...
            temp_f_valuupdate ] =...
            fcn_calc_statistical_func_n_valu(...
            handles.ts( f_idx ).raw_veldata,...
            handles.ts( f_idx ).vals,...
            handles.ts( f_idx ).p,...
            handles.ts( f_idx ).raw_t,...
            handles.ts( f_idx ).raw_corsnr_xt,...
            handles.ts_defaults,...
            size( handles.ts( f_idx ).valu_data ));
        %Store additional seri_x
        handles.ts( f_idx ).seri_x.t    = new_f_seri_x.t;
        handles.ts( f_idx ).seri_x.lag  = new_f_seri_x.lag;
        handles.ts( f_idx ).seri_x.tlag = new_f_seri_x.tlag;
        handles.ts( f_idx ).seri_x.f    = new_f_seri_x.f;
        %Store newly calculated valu (non-nan values)
        handles.ts( f_idx ).valu_data( temp_f_valuupdate==1 ) =...
            new_f_valu_data( temp_f_valuupdate==1 );
        handles.ts( f_idx ).valu_enable( temp_f_valuupdate==1 ) =...
            new_f_valu_enable( temp_f_valuupdate==1 );
        %Update log
        if isempty( find( handles.ts( f_idx ).vals==0, 1 ))~=1
            % If bad data interpolation occured:
            handles.log{end+1,1} =...
                sprintf('    # %s [%2.2f, %2.2f, %2.2f] %s',...
                'velocity values corrected by',...
                handles.ts( f_idx ).intl.mean_correctn,...
                'to get zero mean after interpolation');
        end
        %Create message:
        handles.log{end+1,1} =...
            sprintf('    # statistical features recalculated');
        %Display message
        if isvalid( handles.lbox )==1
            set(handles.lbox,'String', handles.log,'Value',length(handles.log));
        end
        %DVL:
        if isfield( handles, 'dvl')==1
            [ handles, ~ ] = viper_dvl('dvl_calc_ser_n_val_reset',...
                handles, {[]});
        end
        % -----------------------------------------------------------------
        % Calc. new plot characteristics (axlims_minmax):
        % -----------------------------------------------------------------
        [ handles.plots.data_info, handles.plots.axlims_minmax ] =...
            gui_plots_calc_axlims_minmax( handles);
        % -----------------------------------------------------------------
        % Save data:
        % ------------------------------------------------------------
        % Update handles structure
        guidata(hObject, handles);
        gui_update_ui_values(hObject, handles)
        gui_plots_update(hObject, handles, 0 );
    end
    % ------------------------------------------------------------
    % If applying == switching EDIT to off
    % ------------------------------------------------------------
    % Change display name to "Edit":
    set(handles.ui_errf_edit_n_cancel_toggle_b,'String','Edit')
    set(handles.ui_errf_edit_n_cancel_toggle_b,'Value',0)
    gui_update_ui_enables(hObject, handles)
    
    close(wbar);
end
% ------------------------------------------------------------
% If Editing DEFAULT Err. filtering switched on:
% ------------------------------------------------------------
if get(handles.ui_errf_setdefaults_togglebutton, 'Value')==1
    wbar = waitbar(0.01,...
        'Checking parameters & saving as new defaults, please wait...',...
        'WindowStyle','modal');
    % ------------------------------------------------------------
    % Final check of PST-Filter: if all comp==0, pst->0
    % ------------------------------------------------------------
    if get(handles.ui_errf_pst_cmp_u_checkbox, 'Value')==0 &&...
            get(handles.ui_errf_pst_cmp_v_checkbox, 'Value')==0 &&...
            get(handles.ui_errf_pst_cmp_w_checkbox, 'Value')==0
        set(handles.ui_errf_pst_usehpf_checkbox, 'Value',0)
        set(handles.ui_errf_pst_checkbox, 'Value',0)
    end
    % Collect displayed err.filt parameters:
    [ disped_p_errfilt ] = gui_errfiltpars_collect( handles );
    % Store new errfilt pars as defaults:
    handles.ts_defaults.p.errfilt = disped_p_errfilt;
    %log:
    handles.log{end+1,1} =...
        sprintf('~~~ NEW DEFAULT err.filtering parameters saved ~~~');
    if isvalid( handles.lbox )==1
        set(handles.lbox,'String', handles.log,'Value',length(handles.log));
    end
    % Update handles structure
    guidata(hObject, handles);
    % Change display name to "Edit":
    set(handles.ui_errf_setdefaults_togglebutton,'String','Config.')
    set(handles.ui_errf_setdefaults_togglebutton,'Value',0)
    % Fill errf ui_controls with values:
    gui_update_ui_values(hObject, handles)
    % Enable uicontrols:
    gui_update_ui_enables(hObject, handles)
    
    close(wbar);
    msgbox(sprintf('%s\n%s',...
        'New default error-filtering parameters saved!',...
        'New defaults will be applied to files loaded from now on!'),...
        'Config err.filt.','help','modal')
    
end






function ui_errf_setdefaults_togglebutton_Callback(hObject, eventdata, handles)
% ------------------------------------------------------------
% If switching EDIT on:
% ------------------------------------------------------------
if get(handles.ui_errf_setdefaults_togglebutton, 'Value')==1
    msgbox(sprintf('%s%s\n\n%s%s\n%s',...
        'Now you can change the error-filtering paramater defaults ',...
        'by editing the parameters and clicking on ''Apply''.',...
        'New defaults will be applied to files imported after the ',...
        'change.',...
        'Note that if some parmeters are not compatible with a new file ',...
        '(e.g. to non-ADV data), then the file will be imported ',...
        'without error-filtering and can be error-filtered later.'),...
        'Config err.filt.','help','modal')
    % Change display name to "Cancel":
    set(handles.ui_errf_setdefaults_togglebutton,'String','Cancel')
% ------------------------------------------------------------
% If cancelling == switching EDIT to off
% ------------------------------------------------------------
else
    % Change display name to "Edit":
    set(handles.ui_errf_setdefaults_togglebutton,'String','Config.')
end
% Fill errf ui_controls with values:
gui_update_ui_values(hObject, handles)
% Enable uicontrols:
gui_update_ui_enables(hObject, handles)



% ------------------------------------------------------------------------
% Duration-error: switch to relative errors
% ------------------------------------------------------------------------
function ui_dure_rel_radiobutton_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    set(handles.ui_dure_abs_radiobutton,'Value',0)
else
    set(hObject,'Value',1)
end
% Update GUI:
gui_update_ui_values(hObject, handles)

% ------------------------------------------------------------------------
% Duration-error: switch to absolute errors
% ------------------------------------------------------------------------
function ui_dure_abs_radiobutton_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    set(handles.ui_dure_rel_radiobutton,'Value',0)
else
    set(hObject,'Value',1)
end
% Update GUI:
gui_update_ui_values(hObject, handles)










% ------------------------------------------------------------------------
% Coordinates: file coord x
% ------------------------------------------------------------------------
function ui_cds_fil_x_edit_Callback(hObject, eventdata, handles)
input = str2double( regexprep( get(hObject,'String'), ',', '.') );
%Fill final value, Value=1 and color mark new input
set(hObject,'String', num2str(input),...
    'Value',1, 'BackgroundColor',[1,0.6,0.6] )

% ------------------------------------------------------------------------
% Coordinates: file coord y
% ------------------------------------------------------------------------
function ui_cds_fil_y_edit_Callback(hObject, eventdata, handles)
input = str2double( regexprep( get(hObject,'String'), ',', '.') );
set(hObject,'String', num2str(input),...
    'Value',1, 'BackgroundColor',[1,0.6,0.6] )

% ------------------------------------------------------------------------
% Coordinates: file coord z
% ------------------------------------------------------------------------
function ui_cds_fil_z_edit_Callback(hObject, eventdata, handles)
input = str2double( regexprep( get(hObject,'String'), ',', '.') );
set(hObject,'String', num2str(input),...
    'Value',1, 'BackgroundColor',[1,0.6,0.6] )

% ------------------------------------------------------------------------
% Coordinates: mpt offset x
% ------------------------------------------------------------------------
function ui_cds_ofs_x_edit_Callback(hObject, eventdata, handles)
input = str2double( regexprep( get(hObject,'String'), ',', '.') );
set(hObject,'String', num2str(input),...
    'Value',1, 'BackgroundColor',[1,0.6,0.6] )

% ------------------------------------------------------------------------
% Coordinates: mpt offset y
% ------------------------------------------------------------------------
function ui_cds_ofs_y_edit_Callback(hObject, eventdata, handles)
input = str2double( regexprep( get(hObject,'String'), ',', '.') );
set(hObject,'String', num2str(input),...
    'Value',1, 'BackgroundColor',[1,0.6,0.6] )

% ------------------------------------------------------------------------
% Coordinates: mpt offset z
% ------------------------------------------------------------------------
function ui_cds_ofs_z_edit_Callback(hObject, eventdata, handles)
input = str2double( regexprep( get(hObject,'String'), ',', '.') );
set(hObject,'String', num2str(input),...
    'Value',1, 'BackgroundColor',[1,0.6,0.6] )

% ------------------------------------------------------------------------
% Coordinates: mpt coord x
% ------------------------------------------------------------------------
function ui_cds_mpt_x_edit_Callback(hObject, eventdata, handles)
input = str2double( regexprep( get(hObject,'String'), ',', '.') );
set(hObject,'String', num2str(input),...
    'Value',1, 'BackgroundColor',[1,0.6,0.6] )

% ------------------------------------------------------------------------
% Coordinates: mpt coord y
% ------------------------------------------------------------------------
function ui_cds_mpt_y_edit_Callback(hObject, eventdata, handles)
input = str2double( regexprep( get(hObject,'String'), ',', '.') );
set(hObject,'String', num2str(input),...
    'Value',1, 'BackgroundColor',[1,0.6,0.6] )

% ------------------------------------------------------------------------
% Coordinates: mpt coord z
% ------------------------------------------------------------------------
function ui_cds_mpt_z_edit_Callback(hObject, eventdata, handles)
input = str2double( regexprep( get(hObject,'String'), ',', '.') );
set(hObject,'String', num2str(input),...
    'Value',1, 'BackgroundColor',[1,0.6,0.6] )

% ------------------------------------------------------------------------
% Coordinates: probe-u-direction +x
% ------------------------------------------------------------------------
function ui_cds_dir_posx_radiobutton_Callback(hObject, eventdata, handles)
%posx negx posy negy - at least one has to be checked
if get(hObject, 'Value')==1
%     set(handles.ui_cds_dir_posx_radiobutton, 'Value',0)
    set(handles.ui_cds_dir_negx_radiobutton, 'Value',0)
    set(handles.ui_cds_dir_posy_radiobutton, 'Value',0)
    set(handles.ui_cds_dir_negy_radiobutton, 'Value',0)
else
    set(hObject, 'Value',1)
end

% ------------------------------------------------------------------------
% Coordinates: probe-u-direction -x
% ------------------------------------------------------------------------
function ui_cds_dir_negx_radiobutton_Callback(hObject, eventdata, handles)
%posx negx posy negy - at least one has to be checked
if get(hObject, 'Value')==1
    set(handles.ui_cds_dir_posx_radiobutton, 'Value',0)
%     set(handles.ui_cds_dir_negx_radiobutton, 'Value',0)
    set(handles.ui_cds_dir_posy_radiobutton, 'Value',0)
    set(handles.ui_cds_dir_negy_radiobutton, 'Value',0)
else
    set(hObject, 'Value',1)
end

% ------------------------------------------------------------------------
% Coordinates: probe-u-direction +y
% ------------------------------------------------------------------------
function ui_cds_dir_posy_radiobutton_Callback(hObject, eventdata, handles)
%posx negx posy negy - at least one has to be checked
if get(hObject, 'Value')==1
    set(handles.ui_cds_dir_posx_radiobutton, 'Value',0)
    set(handles.ui_cds_dir_negx_radiobutton, 'Value',0)
%     set(handles.ui_cds_dir_posy_radiobutton, 'Value',0)
    set(handles.ui_cds_dir_negy_radiobutton, 'Value',0)
else
    set(hObject, 'Value',1)
end

% ------------------------------------------------------------------------
% Coordinates: probe-u-direction -x
% ------------------------------------------------------------------------
function ui_cds_dir_negy_radiobutton_Callback(hObject, eventdata, handles)
%posx negx posy negy - at least one has to be checked
if get(hObject, 'Value')==1
    set(handles.ui_cds_dir_posx_radiobutton, 'Value',0)
    set(handles.ui_cds_dir_negx_radiobutton, 'Value',0)
    set(handles.ui_cds_dir_posy_radiobutton, 'Value',0)
%     set(handles.ui_cds_dir_negy_radiobutton, 'Value',0)
else
    set(hObject, 'Value',1)
end





% ------------------------------------------------------------------------
% Coordinates: input Edit OR Cancel
% ------------------------------------------------------------------------
function ui_cds_edit_n_cancel_togglebutton_Callback(hObject, eventdata, handles)
%Get selected file index
f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
% ------------------------------------------------------------
% If switching on EDIT:
% ------------------------------------------------------------
if get(handles.ui_cds_edit_n_cancel_togglebutton, 'Value')==1
    % Change display name to "Cancel":
    set(handles.ui_cds_edit_n_cancel_togglebutton,'String','Cancel')
    % Warnings
    h3 = warndlg(sprintf(...
        '%s%s%s',...
        'Best practice is to first use ''Read from filename'', because ',...
        'values of ''Measurement Point coords'' are automatically ',...
        'recalculated and overwritten after using it!'),...
        'Warning!','modal');
    waitfor(h3)
    if length( f_idx )>1
        h1 = warndlg(sprintf(...
            '%s\n\n%s%s\n\n%s%s%s\n\n',...
            'Changes will be applied to EACH selected file as follows:',...
            '1. New values will OVERWRITE old values of same type AND ',...
            'same component for each selected file!',...
            '2. When changing values of ''File coords'' or ''offset'', ',...
            'then values of ''Measurement Point coords'' are ',...
            'AUTOMATICALLY recalculated and overwritten!'),...
            'Warning!','modal');
        waitfor(h1)
    else
        h2 = warndlg(sprintf(...
            '%s%s%s',...
            'When changing values of ''File coords'' or ''offset'', ',...
            'then values of ''Measurement Point coords'' are ',...
            'AUTOMATICALLY recalculated and overwritten!'),...
            'Warning!','modal');
        waitfor(h2)
    end
% ------------------------------------------------------------
% If cancelling == switching off EDIT
% ------------------------------------------------------------
else
    %Reset editboxes
    set( handles.ui_cds_fil_x_edit,...
        'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
    set( handles.ui_cds_fil_y_edit,...
        'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
    set( handles.ui_cds_fil_z_edit,...
        'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
    set( handles.ui_cds_ofs_x_edit,...
        'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
    set( handles.ui_cds_ofs_y_edit,...
        'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
    set( handles.ui_cds_ofs_z_edit,...
        'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
    set( handles.ui_cds_mpt_x_edit,...
        'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
    set( handles.ui_cds_mpt_y_edit,...
        'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
    set( handles.ui_cds_mpt_z_edit,...
        'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
    % Change display name to "Edit":
    set(handles.ui_cds_edit_n_cancel_togglebutton,'String','Edit')
    % Cancel changes by loading stored parameters:
    gui_update_ui_values(hObject, handles)
end
% Enable uicontrols:
gui_update_ui_enables(hObject, handles)





% ------------------------------------------------------------------------
% Coordinates: input apply
% ------------------------------------------------------------------------
function ui_cds_apply_pushbutton_Callback(hObject, eventdata, handles)
%Get selected file index
f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
%Create temporary containers:
input_cds_fil = [NaN,NaN,NaN];
input_cds_ofs = [NaN,NaN,NaN];
input_cds_mpt = [NaN,NaN,NaN];
%Get entered values (NaN input is not registered as newly entered)
if get( handles.ui_cds_fil_x_edit, 'Value')==1
    input_cds_fil(1) = str2double( ...
        regexprep( get(handles.ui_cds_fil_x_edit ,'String'), ',', '.') );
end
if get( handles.ui_cds_fil_y_edit, 'Value')==1
    input_cds_fil(2) = str2double( ...
        regexprep( get(handles.ui_cds_fil_y_edit ,'String'), ',', '.') );
end
if get( handles.ui_cds_fil_z_edit, 'Value')==1
    input_cds_fil(3) = str2double( ...
        regexprep( get(handles.ui_cds_fil_z_edit ,'String'), ',', '.') );
end
if get( handles.ui_cds_ofs_x_edit, 'Value')==1
    input_cds_ofs(1) = str2double( ...
        regexprep( get(handles.ui_cds_ofs_x_edit ,'String'), ',', '.') );
end
if get( handles.ui_cds_ofs_y_edit, 'Value')==1
    input_cds_ofs(2) = str2double( ...
        regexprep( get(handles.ui_cds_ofs_y_edit ,'String'), ',', '.') );
end
if get( handles.ui_cds_ofs_z_edit, 'Value')==1
    input_cds_ofs(3) = str2double( ...
        regexprep( get(handles.ui_cds_ofs_z_edit ,'String'), ',', '.') );
end
if get( handles.ui_cds_mpt_x_edit, 'Value')==1
    input_cds_mpt(1) = str2double( ...
        regexprep( get(handles.ui_cds_mpt_x_edit ,'String'), ',', '.') );
end
if get( handles.ui_cds_mpt_y_edit, 'Value')==1
    input_cds_mpt(2) = str2double( ...
        regexprep( get(handles.ui_cds_mpt_y_edit ,'String'), ',', '.') );
end
if get( handles.ui_cds_mpt_z_edit, 'Value')==1
    input_cds_mpt(3) = str2double( ...
        regexprep( get(handles.ui_cds_mpt_z_edit ,'String'), ',', '.') );
end
pos_new_cds_fil = find(isnan(input_cds_fil)==0);
pos_new_cds_ofs = find(isnan(input_cds_ofs)==0);
pos_new_cds_mpt = find(isnan(input_cds_mpt)==0);
%Create string:
cds_fil_str = num2str(input_cds_fil);
cds_fil_str = regexprep( cds_fil_str , 'NaN', 'n.a.');
cds_ofs_str = num2str(input_cds_ofs);
cds_ofs_str = regexprep( cds_ofs_str , 'NaN', 'n.a.');
cds_mpt_str = num2str(input_cds_mpt);
cds_mpt_str = regexprep( cds_mpt_str , 'NaN', 'n.a.');
if get(handles.ui_cds_dir_posx_radiobutton,'Value')==1
    new_dir = '+x';
elseif get(handles.ui_cds_dir_negx_radiobutton,'Value')==1
    new_dir = '-x';
elseif get(handles.ui_cds_dir_posy_radiobutton,'Value')==1
    new_dir = '+y';
elseif get(handles.ui_cds_dir_negy_radiobutton,'Value')==1
    new_dir = '-y';
end
% -------------------------------------------------------------------------
% Test and confirm on first filename:
% -------------------------------------------------------------------------
%Ask to overwrite stored values:
sure = questdlg(sprintf(...
    '%s\n%s\n\n%s %s \n%s %s \n%s %s',...
    'This input overwrites the following values for EACH selected file:',...
    '(Please consider the currently used units of the GUI!)',...
    'File coords: ', cds_fil_str,...
    'MP offset:   ', cds_ofs_str,...
    'MP coords:   ', cds_mpt_str ),...
    'Apply new coordinates?',...
    'Continue','Cancel','Cancel');
if strcmp(sure,'Cancel')==1 || strcmp(sure,'')
    return
end
%Log
handles.log{end+1,1} =...
    sprintf('>>> Setting new coord. system properties manually:');
handles.log{end+1,1} =...
    sprintf('    * File coords input: [ %s ]', cds_fil_str);
handles.log{end+1,1} =...
    sprintf('    * MP. offset input: [ %s ]', cds_ofs_str);
handles.log{end+1,1} =...
    sprintf('    * MP. coords input: [ %s ]', cds_mpt_str);
handles.log{end+1,1} =...
    sprintf('    * Direction of ''Probe-X'': [ %s ]', new_dir);
if isvalid( handles.lbox )==1
    set(handles.lbox,'String',handles.log,'Value',length(handles.log));
end
%%%% Apply coordinates - at his point already checked for validity
% Applying to all selected files:
wbar = waitbar(0.01,'Applying coordinates to all files, please wait...',...
    'WindowStyle','modal');
list_err_files  = {};
for fik=1:length(f_idx)
    %Log
    handles.log{end+1,1} =...
        sprintf('  > %s', handles.fname{ f_idx(fik) } );
    if isvalid( handles.lbox )==1
        set(handles.lbox,'String',handles.log,'Value',length(handles.log));
    end
    
    %Store new coordinates (=if non-nan) - else do nothing
    % -----------------------------------------------------------------
    %If file_coords got manual input
    if isempty( pos_new_cds_fil )~=1
        handles.ts( f_idx(fik) ).intl.coordsys.file_coords(pos_new_cds_fil) =...
            input_cds_fil(pos_new_cds_fil);
    end
    % -----------------------------------------------------------------
    %If mpt_offset got manual input
    if isempty( pos_new_cds_ofs )~=1
        handles.ts( f_idx(fik) ).intl.coordsys.measpt_offset(pos_new_cds_ofs)=...
            input_cds_ofs(pos_new_cds_ofs);
    end
    % -----------------------------------------------------------------
    % Recalculate cds_mpt ANYWAY AND
    %    if cds_fil OR cds_ofs changed and are non-NAN -> 
    %    == mark cds_mpt as recalculated
    calc_cds_mpt =...
        handles.ts( f_idx(fik) ).intl.coordsys.file_coords +...
        handles.ts( f_idx(fik) ).intl.coordsys.measpt_offset;
    %old value:
    %old_cds_mpt = handles.ts( f_idx(fik) ).intl.coordsys.measpt_coords;
    if isempty( pos_new_cds_fil )~=1 || isempty( pos_new_cds_ofs )~=1
        pos2cal = find(...
            isnan(input_cds_fil)==0 | isnan(input_cds_ofs)==0);
    else
        pos2cal = [];
    end
    %mark as to_calculate
    marker_newcal_cds_mpt = nan(1,3);
    marker_newcal_cds_mpt(pos2cal) = 1;
    % Conflict check and overwrite cds_mpt using calculated cds_mpt
    % -----------------------------------------------------------------
    % if mpt_coords did not get manual input -> 
    % -> use old or non-nan calculated
    if isempty( pos_new_cds_mpt )==1
        pos_val = isnan( calc_cds_mpt )==0;
        handles.ts( f_idx(fik) ).intl.coordsys.measpt_coords(pos_val) =...
            calc_cds_mpt(pos_val);
        %Log
        handles.log{end+1,1} =...
            sprintf('    # set coords & offsets' );
    % -----------------------------------------------------------------
    % if mpt_coords got manual input:
    else
        %Conflicts to check:
        %(mpt positions with manual input fits to calculated mpt)
        pos2check = find(...
            isnan( input_cds_mpt )==0 & isnan( calc_cds_mpt )==0 );
        %If there is a conflict:
        if isempty(pos2check)~=1 &&...
                isequal(...
                calc_cds_mpt(pos2check), input_cds_mpt(pos2check) )~=1
            list_err_files{1,end+1} = handles.fname{ f_idx(fik) };
            % -> do not apply mpt_coords, but use recalculated!
            handles.ts( f_idx(fik) ).intl.coordsys.measpt_coords =...
                calc_cds_mpt;
            %Log:
            handles.log{end+1,1} =...
                sprintf('    ! NOT applying entered MP coords' );
        else%There is no conflict -> apply manually-set coords
            handles.ts( f_idx(fik) ).intl.coordsys.measpt_coords(pos_new_cds_mpt)=...
                input_cds_mpt(pos_new_cds_mpt);
            %Log
            handles.log{end+1,1} =...
                sprintf('    # set coords & offsets' );
        end
    end
    %update .p.coordsys by adjusting units:
    handles.ts( f_idx(fik) ).p.coordsys.file_coords =...
        handles.ts( f_idx(fik) ).intl.coordsys.file_coords /...
        handles.gui.coordunitscale;
    handles.ts( f_idx(fik) ).p.coordsys.measpt_offset =...
        handles.ts( f_idx(fik) ).intl.coordsys.measpt_offset /...
        handles.gui.coordunitscale;
    handles.ts( f_idx(fik) ).p.coordsys.measpt_coords =...
        handles.ts( f_idx(fik) ).intl.coordsys.measpt_coords /...
        handles.gui.coordunitscale;
    %Log
    if isvalid( handles.lbox )==1
        set(handles.lbox,'String',handles.log,'Value',length(handles.log));
    end
    % -----------------------------------------------------------------
    % velocity direction input - only if there is a change
    if strcmp( handles.ts( f_idx(fik) ).intl.coordsys.dir_probe_u,...
            new_dir)~=1
        dir_old = handles.ts( f_idx(fik) ).intl.coordsys.dir_probe_u;
        only_raw = 0;
        %rotate u-v vector in affected variables
            %  - rotate raw_veldata_uvw & intl
            %  - rotate seri & valu data
            %  - rotate parameters and arrfilt res
            %  - also copies dir_probe_u to .p
        [ handles.ts( f_idx(fik) ).raw_veldata,...
            handles.ts( f_idx(fik) ).seri_data,...
            handles.ts( f_idx(fik) ).valu_data,...
            handles.ts( f_idx(fik) ).intl,...
            handles.ts( f_idx(fik) ).p,...
            handles.ts( f_idx(fik) ).efilt_res, alpha ] =...
            fcn_calc_rotate_uv(...
              only_raw, dir_old, new_dir,...
              handles.ts( f_idx(fik) ).raw_veldata,...
              handles.ts( f_idx(fik) ).seri_data,...
              handles.ts( f_idx(fik) ).valu_data,...
              handles.ts( f_idx(fik) ).intl,...
              handles.ts( f_idx(fik) ).p,...
              handles.ts( f_idx(fik) ).efilt_res,...
              handles.ts( f_idx(fik) ).seri_enable, handles.field_props );
        %recalculate hist:
        [ handles.ts( f_idx(fik) ).hist.p_per_bin,...
          handles.ts( f_idx(fik) ).hist.vel_of_bin ] =...
          fcn_calc_histograms_of_vals(...
            handles.ts( f_idx(fik) ).raw_veldata,...
            handles.ts( f_idx(fik) ).vals, ...
            handles.ts( f_idx(fik) ).p.errfilt.used_velcomp );
        %Recalc. plot characteristics (axlims_minmax):
        [ handles.plots.data_info, handles.plots.axlims_minmax ] =...
            gui_plots_calc_axlims_minmax( handles);
        %Log
        handles.log{end+1,1} =...
            sprintf('    # set probe u-direction' );
        %rotate in DVL 
        if isfield( handles, 'dvl')==1
            [ handles, ~ ] = viper_dvl('dvl_calc_rotate_uv',...
                handles, { f_idx(fik), alpha } );
        end
    end
    if isempty(find(...
            isnan(handles.ts( f_idx(fik) ).p.coordsys.measpt_coords),1))==1
        handles.ts( f_idx(fik) ).enable_coords = {'On'};
    else
        handles.ts( f_idx(fik) ).enable_coords = {'Off'};
    end
    if isvalid( handles.lbox )==1
        set(handles.lbox,'String',handles.log,'Value',length(handles.log));
    end
    waitbar( fik/length(f_idx), wbar)
end
%List errors - if any
if isempty(list_err_files)~=1
    warndlg(...
        ['Not applied entered ''MP coords'' for following files ',...
        '(Conflicting MP coords input & calculations)',...
        list_err_files],'Skipped file(s)!');
end
%Closing waitbar:
close(wbar);
% Update handles structure
guidata(hObject, handles);
%Reset editboxes
set( handles.ui_cds_fil_x_edit,...
    'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
set( handles.ui_cds_fil_y_edit,...
    'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
set( handles.ui_cds_fil_z_edit,...
    'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
set( handles.ui_cds_ofs_x_edit,...
    'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
set( handles.ui_cds_ofs_y_edit,...
    'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
set( handles.ui_cds_ofs_z_edit,...
    'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
set( handles.ui_cds_mpt_x_edit,...
    'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
set( handles.ui_cds_mpt_y_edit,...
    'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
set( handles.ui_cds_mpt_z_edit,...
    'Value',0,'BackgroundColor',[0.941,0.941,0.941] );
%Switch off edit mode
set(handles.ui_cds_edit_n_cancel_togglebutton,'String','Edit')
set(handles.ui_cds_edit_n_cancel_togglebutton,'Value',0)
gui_update_ui_enables(hObject, handles)
%Update GUI:
gui_update_ui_values(hObject, handles)
gui_plots_update(hObject, handles, 0 );









% ------------------------------------------------------------------------
% Coordinates: Read coordinates from filename
% ------------------------------------------------------------------------
function ui_cds_getfromfile_pushbutton_Callback(hObject,eventdata,handles)
h3 = warndlg(sprintf(...
    '%s%s%s',...
    'When changing new values for ''File coords'' or ''offset'', ',...
    'then values of ''Measurement Point coords'' are ',...
    'automatically recalculated and overwritten!'),...
    'Warning!','modal');
waitfor(h3)
%Get selected file index
f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
% -------------------------------------------------------------------------
% User input: coords format string and unit:
% -------------------------------------------------------------------------
usr_input_str{1} = '';
usr_input_str{2} = '';
repeat_input = 0;
while isempty(usr_input_str)==1 || isempty(usr_input_str{1})==1 || ...
        isempty(usr_input_str{2})==1 || repeat_input == 0
    usr_input_str = inputdlg(...
        { sprintf(...
        'Please replace %s , %s (Filename example: %s)',...
        'x- y- z-coordinate strings by ''%x%'', ''%y%'' and ''%z%'' ',...
        'and measurement point offset by ''%xm%'', ''%ym%'' and ''%zm%''.',...
        handles.fname{ f_idx(1) }),...
        'Units of these coordinates (use: 1000==m, 10==cm, 1==mm)'},...
        'Get coordinates from filename',1,...
        { handles.viper_props.userdef.cds_fmt,...
        num2str(handles.viper_props.userdef.cds_unit)} );
    repeat_input = 1;
    % -------------------------------------------------------------------------
    % Check user inputs:
    if isempty( usr_input_str )==1%==Cancelled
        return
    elseif isempty( usr_input_str{1} )==1
        dlg_h = errordlg('Missing input: coordinate-format string!',...
            'Error');
        waitfor(dlg_h);
        repeat_input = 0;
        continue
    elseif isempty( usr_input_str{2} )==1
        dlg_h = errordlg('Missing input for unit!','Error');
        waitfor(dlg_h);
        repeat_input = 0;
        continue
    end
    cds_fmt_str = usr_input_str{1};
    cds_unit = str2double( regexprep( usr_input_str{2}, ',', '.') );
    if isnan( cds_unit )==1
        if strcmp(usr_input_str{2},'mm')
            cds_unit=1;
        elseif strcmp(usr_input_str{2},'cm')
            cds_unit=10;
        elseif strcmp(usr_input_str{2},'m')
            cds_unit=1000;
        else
            dlg_h = errordlg('Incorrect input for unit!', 'Error');
            waitfor(dlg_h);
            repeat_input = 0;
            continue
        end
    elseif (cds_unit==1 || cds_unit==10 || cds_unit==1000)~=1
        dlg_h = errordlg('Incorrect input for unit!', 'Error');
        waitfor(dlg_h);
        repeat_input = 0;
        continue
    end
end
% -------------------------------------------------------------------------
% File-coordinate format string:
if  ~contains(cds_fmt_str, '%x%' )==0 &&...
        ~contains(cds_fmt_str, '%y%' )==0 &&...
        ~contains(cds_fmt_str, '%z%' )==0
    pos_cds(1) = strfind(cds_fmt_str, '%x%' );
    pos_cds(2) = strfind(cds_fmt_str, '%y%' );
    pos_cds(3) = strfind(cds_fmt_str, '%z%' );
    h1 = [];
elseif   ~contains(cds_fmt_str, '%X%' )==0 &&...
        ~contains(cds_fmt_str, '%Y%' )==0 &&...
        ~contains(cds_fmt_str, '%Z%' )==0
    pos_cds(1) = strfind(cds_fmt_str, '%X%' );
    pos_cds(2) = strfind(cds_fmt_str, '%Y%' );
    pos_cds(3) = strfind(cds_fmt_str, '%Z%' );
    h1 = [];
else
    h1 = warndlg('No format-string input for file-coordinates!',...
        'Error');
    pos_cds = [];
end
if isempty( pos_cds )~=1 && isempty(find(diff(pos_cds)<=3,1))~=1
    errordlg('The order of coordinates within the string is not x-y-z!',...
        'Error');
    pos_cds = [];
end
% -------------------------------------------------------------------------
% Meas.-point offset format string:
if  ~contains(cds_fmt_str, '%xm%' )==0 &&...
        ~contains(cds_fmt_str, '%ym%' )==0 &&...
        ~contains(cds_fmt_str, '%zm%' )==0
    pos_ofs(1) = strfind(cds_fmt_str, '%xm%' );
    pos_ofs(2) = strfind(cds_fmt_str, '%ym%' );
    pos_ofs(3) = strfind(cds_fmt_str, '%zm%' );
    h2 = [];
elseif ~contains(cds_fmt_str, '%Xm%' )==0 &&...
        ~contains(cds_fmt_str, '%Ym%' )==0 &&...
        ~contains(cds_fmt_str, '%Zm%' )==0
    pos_ofs(1) = strfind(cds_fmt_str, '%Xm%' );
    pos_ofs(2) = strfind(cds_fmt_str, '%Ym%' );
    pos_ofs(3) = strfind(cds_fmt_str, '%Zm%' );
    h2 = [];
else
    h2 = warndlg('No format-string input for measurement-point offset!',...
        'Error');
    pos_ofs = [];
end
if isempty( pos_ofs )~=1 && isempty(find(diff(pos_ofs)<=3,1))~=1
    errordlg('The order of offset within the string is not x-y-z!',...
        'Error');
    pos_ofs = [];
end
if isempty( pos_cds ) ==1 && isempty( pos_ofs ) ==1
    return
end
%Coords unitscale to get in hui_units:
cds_unitscale = cds_unit/handles.gui.units_id;
% -------------------------------------------------------------------------
% Test and confirm on first filename:
% -------------------------------------------------------------------------
%Get stored values for first file in selection - as default!
file_coords = handles.ts( f_idx(1) ).intl.coordsys.file_coords;
measpt_offset = handles.ts( f_idx(1) ).intl.coordsys.measpt_offset;
status_fil_cds = 'no change';
status_mpt_ofs = 'no change';
status_mpt_cds = 'no change';
%Read values from filename
if isempty( pos_cds )~=1
    % Get file-coordinates from filename of current file:
    read_coords = fcn_cds_read_from_str_sokoray(...
        pos_cds, cds_fmt_str, handles.fname{ f_idx(1) } );
%     read_coords(2)
    new_coords = read_coords*cds_unitscale;%
    %If worked -> set as editbox value AND ask doing it for the other files
    if isempty(find(isnan( new_coords ),1)) ~=0
        file_coords = new_coords;
        status_fil_cds = 'NEW values';
        status_mpt_cds = 'NEW values';
    end
end
%Read values from filename
if isempty( pos_ofs )~=1
    % Get meas. point offset from filename of current file:
    read_ofset = fcn_cds_read_from_str_sokoray(...
        pos_ofs, cds_fmt_str, handles.fname{ f_idx(1) } );
    new_offset = read_ofset*cds_unitscale;%
    if isempty(find(isnan( new_offset ),1)) ~=0
        measpt_offset = new_offset;
        status_mpt_ofs = 'NEW values';
        status_mpt_cds = 'NEW values';
    end
end
% return
%Calculate meas. point coordinates:
calc_measpt_coords = file_coords + measpt_offset;
waitfor(h1)
waitfor(h2)
%Ask to overwrite stored values:
sure = questdlg(sprintf(...
    '%s\n%s\n\n%s %.1f %.1f %.1f \n%s %.1f %.1f %.1f \n%s %.1f %.1f %.1f',...
    'This input overwrites values for EACH selected file as follows:',...
    '(The following examples apply to the first selected file)',...
    ['File coords: (' status_fil_cds ')'], file_coords,...
    ['MP offset:   (' status_mpt_ofs ')'], measpt_offset,...
    ['MP coords:   (' status_mpt_cds ')'], calc_measpt_coords ),...
    'Apply new coordinates?',...
    'Continue','Cancel','Cancel');
if strcmp(sure,'Cancel')==1 || strcmp(sure,'')
    return
end
%Store format as default:
handles.viper_props.userdef.cds_fmt = cds_fmt_str;
handles.viper_props.userdef.cds_unit = cds_unit;%1000==m, 10==cm, 1==mm
%Log
handles.log{end+1,1} =...
    sprintf('>>> Reading coordinates from filename(s) using:');
handles.log{end+1,1} =...
    sprintf('    * format: %s', cds_fmt_str);
handles.log{end+1,1} =...
    sprintf('    * coordinate units: %d', cds_unit);
if isvalid( handles.lbox )==1
    set(handles.lbox,'String',handles.log,'Value',length(handles.log));
end
% ------------------------------------------------------------
% Applying to all selected files:
% ------------------------------------------------------------
wbar = waitbar(0.01,'Reading coordinates from filenames, please wait...',...
    'WindowStyle','modal');
list_err_files  = {};
for fik=1:length(f_idx)
    %Log
    handles.log{end+1,1} =...
        sprintf('  > %s', handles.fname{ f_idx(fik) } );
    if isvalid( handles.lbox )==1
        set(handles.lbox,'String',handles.log,'Value',length(handles.log));
    end
    %Read stored values:
    file_coords = handles.ts( f_idx(fik) ).intl.coordsys.file_coords;
    measpt_offset = handles.ts( f_idx(fik) ).intl.coordsys.measpt_offset;
    %Read values from filename
    if isempty( pos_cds )~=1
        % Get file-coordinates from filename of current file:
        read_coords = fcn_cds_read_from_str_sokoray(...
            pos_cds, cds_fmt_str, handles.fname{ f_idx(fik) } );
        new_coords = read_coords*cds_unitscale;%
        %If worked -> set as editbox value AND ask doing it for other files
        if isempty(find(isnan( new_coords ),1)) ~=0
            file_coords = new_coords;
        else
            list_err_files{1,end+1} = handles.fname{ f_idx(fik) };
        end
    end
    %Read values from filename
    if isempty( pos_ofs )~=1
        % Get meas. point offset from filename of current file:
        read_ofset = fcn_cds_read_from_str_sokoray(...
            pos_ofs, cds_fmt_str, handles.fname{ f_idx(fik) } );
        new_offset = read_ofset*cds_unitscale;%
        if isempty(find(isnan( new_offset ),1)) ~=0
            measpt_offset = new_offset;
        else
            %if not yet on list - else no need
            if strcmp(list_err_files{1,end},handles.fname{ f_idx(fik) })~=1
                list_err_files{1,end+1} = handles.fname{ f_idx(fik) };
            end
        end
    end
    %if no error:
    if isempty( list_err_files )==1 ||...
        strcmp(list_err_files{1,end},handles.fname{ f_idx(fik) })~=1
        %Calculate meas. point coordinates:
        calc_measpt_coords = file_coords + measpt_offset;
        %Store coordinates
        handles.ts( f_idx(fik) ).intl.coordsys.file_coords = file_coords;
        handles.ts( f_idx(fik) ).intl.coordsys.measpt_offset = measpt_offset;
        %Store measpt_coords only if non-nan
        pos_val = isnan( calc_measpt_coords )==0;
        handles.ts( f_idx(fik) ).intl.coordsys.measpt_coords(pos_val)=...
            calc_measpt_coords(pos_val);
        %Copy intl -> p
        handles.ts( f_idx(fik) ).p.coordsys.file_coords =...
            handles.ts( f_idx(fik) ).intl.coordsys.file_coords/...
            handles.gui.coordunitscale;
        handles.ts( f_idx(fik) ).p.coordsys.measpt_offset =...
            handles.ts( f_idx(fik) ).intl.coordsys.measpt_offset/...
            handles.gui.coordunitscale;
        handles.ts( f_idx(fik) ).p.coordsys.measpt_coords =...
            handles.ts( f_idx(fik) ).intl.coordsys.measpt_coords/...
            handles.gui.coordunitscale;
        handles.ts( f_idx(fik) ).enable_coords = {'On'};
        %Log
        handles.log{end+1,1} =...
            sprintf('    # done');
        if isvalid( handles.lbox )==1
            set(handles.lbox,'String',handles.log,'Value',length(handles.log));
        end
    else%if error occured -> do nothing just protocol
        %Log
        handles.log{end+1,1} =...
            sprintf('  ! skipped (ERROR)');
        if isvalid( handles.lbox )==1
            set(handles.lbox,'String',handles.log,'Value',length(handles.log));
        end
    end
end
%List errors - if any
if isempty(list_err_files)~=1
    warndlg(['Skipped reading coordinates of file(s) do to errors:',...
        list_err_files],'Skipped file(s)!');
end
%Closing waitbar:
close(wbar);
% Update handles structure
guidata(hObject, handles);
%Update GUI:
gui_update_ui_values(hObject, handles)









% ------------------------------------------------------------------------
% Export: parfile
% ------------------------------------------------------------------------
function ui_exp_pars_pushbutton_Callback(hObject, eventdata, handles)
%Ask if sure
sure = questdlg(sprintf('This exports the paramters of each selected %s',...
    'file to its own par-file! Overwrite check follows!'),...
    'Continue exporting?','Continue','Cancel','Cancel');
if strcmp(sure,'Cancel')==1 || strcmp(sure,'')
    return
end
%create dialog
fig_H = 120;
fig_W = 240;
fig_le = handles.gui.pos_screen(3)/2-fig_W/2;
fig_bo = handles.gui.pos_screen(4)/2-fig_H/2;
dlg_data.ui_fig = dialog('Name','Selecting parameters to export...',...
    'Position', [fig_le fig_bo fig_W fig_H]);
%Plot uicontrols
dlg_data.ui_tx = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 fig_H-25 fig_W-20 21], 'Style', 'text',...
    'String', 'Please select paramter groups to export');%,...
dlg_data.ui_dat_props = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 fig_H-45 fig_W-20 21], 'Style', 'checkbox',...
    'String', 'data properties', 'Value', 1, 'Enable', 'Off');%,...
dlg_data.ui_coordsys = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 fig_H-65 fig_W-20 21], 'Style', 'checkbox',...
    'String', 'coordinate transformation parameters', 'Value', 1);%,...
dlg_data.ui_errfilt = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 fig_H-85 fig_W-20 21], 'Style', 'checkbox',...
    'String', 'error filtering parameters', 'Value', 1);%,...
dlg_data.ui_ok = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 fig_H-110 (fig_W-20)/2 21], 'Style', 'pushbutton',...
    'String', 'OK', 'Value', 1);%,...
dlg_data.ui_cl = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10+(fig_W-20)/2 fig_H-110 (fig_W-20)/2 21],...
    'Style', 'pushbutton',...
    'String', 'Cancel', 'Value', 1);%,...
set(dlg_data.ui_fig,'CloseRequestFcn',...
    {@gui_exportdlg_pars_callbacks,handles.figure1,dlg_data,0});
set(dlg_data.ui_cl, 'Callback',...
    {@gui_exportdlg_pars_callbacks,handles.figure1,dlg_data,0});
set(dlg_data.ui_ok, 'Callback',...
    {@gui_exportdlg_pars_callbacks,handles.figure1,dlg_data,1});





% ------------------------------------------------------------------------
% Export: statistical results
% ------------------------------------------------------------------------
function ui_exp_stat_pushbutton_Callback(hObject, eventdata, handles)
%Ask if sure
sure = questdlg(sprintf('This exports the statistical results of all %s',...
    'selected files into one file! Overwrite check follows!'),...
    'Continue exporting?','Continue','Cancel','Cancel');
if strcmp(sure,'Cancel')==1 || strcmp(sure,'')
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Collect dataset 
%Get selected file index
f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
%Collect data: (keeping the same sorting as in f_idx)
valuarray = nan( length(f_idx), length(handles.ts_defaults.valu_data) );
ena_array = false( length(f_idx), length(handles.ts_defaults.valu_data) );
cds_array = nan( length(f_idx), 3 );
for fik=1:length(f_idx)
    %Enable state of valu-s
    ena_array(fik,:) =handles.ts( f_idx(fik) ).valu_enable;
    %Valu-s
    valuarray(fik,:) =handles.ts( f_idx(fik) ).valu_data;
    %Coordinates:
    cds_array(fik,:) =handles.ts( f_idx(fik) ).intl.coordsys.measpt_coords;
end
%Valu-s that are enabled for all files
posena_valu = min( ena_array,[],1);
%List of files:
fname_list = handles.fname( f_idx );
fpath_list = handles.fpath( f_idx );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check validity
% If no valu is present for all the files:
if sum( posena_valu )==0
    errordlg(sprintf('There are no commonly valid statistical %s',...
        'features for all selected files. Export cancelled!'),...
        '!!! ERROR !!!');
    return
end
% %valu-s that are enabled for all selected files:
% names_posena = handles.field_props.valu_names(posena_valu);
% sure=questdlg(...
%     ['Statistical features that are valid for all selected files are:',...
%     names_posena ],...
%     'Continue?','Continue','Cancel','Cancel');
% if strcmp(sure,'Cancel')==1 || strcmp(sure,'')
%     return
% end
%valu-s that are enabled AND are numeric (==non-nan) ==1
posena_valu_numeric = min( isnan(valuarray(:,posena_valu))==0,[],1);
%%% Later Check for nan's:
dat_typ_enabled = 'off';%not working so far
% dat_typ_enabled = 'on';%default == enable dat export
% If some of the enabled valu-s are NaN:
if isempty( find( posena_valu_numeric == 0, 1))==0
    sure=questdlg(sprintf(...
        'Some of the valid statistical features %s %s %s',...
        'are NaNs. For this reason, the results of this selection ',...
        'can only be exported as Excel sheet OR as Matlab dataset'),...
        'Continue?','Continue','Cancel','Cancel');
        %'unless you replace NaN''s by a numeric value at export'
    if strcmp(sure,'Cancel')==1 || strcmp(sure,'')
        return
    else
        dat_typ_enabled = 'off';
    end
end
%If some meas. point coordinate is NAN:
if isempty(find( isnan( cds_array ), 1))~=1
    sure=questdlg(sprintf(...
        'Some of the measurement point coordinates %s %s %s',...
        'are NaNs. For this reason, the results of this selection ',...
        'can only be exported as Excel sheet OR as Matlab dataset'),...
        'Continue?','Continue','Cancel','Cancel');
        %'unless you replace NaN''s by a numeric value at export'
    if strcmp(sure,'Cancel')==1 || strcmp(sure,'')
        return
    else
        dat_typ_enabled = 'off';
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create dialog
H_fig = 160;
W_fig = 320;
fig_le = handles.gui.pos_screen(3)/2-W_fig/2;
fig_bo = handles.gui.pos_screen(4)/2-H_fig/2;
dlg_data.ui_fig = dialog('Name','Selecting export format...',...
    'Position', [fig_le fig_bo W_fig H_fig]);
%Plot uicontrols
%Select export units
dlg_data.ui_unit_tx = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 H_fig-25 W_fig-20 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String', 'Please select output units:');%,...
dlg_data.ui_unit_m = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 30 H_fig-45 80 21], 'Style', 'radiobutton',...
    'String', 'm, m/s', 'Value', 1);%,...
dlg_data.ui_unit_cm = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 120 H_fig-45 80 21], 'Style', 'radiobutton',...
    'String', 'cm, cm/s', 'Value', 0);%,...
dlg_data.ui_unit_mm = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 210 H_fig-45 80 21], 'Style', 'radiobutton',...
    'String', 'mm, mm/s', 'Value', 0);%,...
H_unit = 50;%needed height for unit - including 5 pixel space on botton
%Select export formats:
dlg_data.ui_file_tx = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 H_fig-H_unit-25 W_fig-20 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String', 'Please select output file format:');
dlg_data.ui_fil_xls = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 30 H_fig-H_unit-45 W_fig-40 21], 'Style', 'radiobutton',...
    'String', '.xls: Excel sheet', 'Value', 1);
dlg_data.ui_fil_mat = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 30 H_fig-H_unit-65 W_fig-40 21], 'Style', 'radiobutton',...
    'String', '.mat: Matlab variable file', 'Value', 0);
H_file = H_unit+75;%needed height for unit - including 5 pixel space on botton
%Buttons:
dlg_data.ui_ok = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 H_fig-H_file-25 (W_fig-20)/2 21],...
    'Style', 'pushbutton',...
    'String', 'OK', 'Value', 0);%,...
dlg_data.ui_cl = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10+(W_fig-20)/2 H_fig-H_file-25 (W_fig-20)/2 21],...
    'Style', 'pushbutton',...
    'String', 'Cancel', 'Value', 1);%,...
%radiobutton callbacks
set(dlg_data.ui_unit_m, 'Callback',...
    {@gui_exportdlg_stat_ui_callbacks,dlg_data,dlg_data.ui_unit_m});
set(dlg_data.ui_unit_cm, 'Callback',...
    {@gui_exportdlg_stat_ui_callbacks,dlg_data,dlg_data.ui_unit_cm});
set(dlg_data.ui_unit_mm, 'Callback',...
    {@gui_exportdlg_stat_ui_callbacks,dlg_data,dlg_data.ui_unit_mm});
set(dlg_data.ui_fil_xls, 'Callback',...
    {@gui_exportdlg_stat_ui_callbacks,dlg_data,dlg_data.ui_fil_xls});
set(dlg_data.ui_fil_mat, 'Callback',...
    {@gui_exportdlg_stat_ui_callbacks,dlg_data,dlg_data.ui_fil_mat});
%Callbacks with dlg_data:
set(dlg_data.ui_ok, 'Callback',...
    {@gui_exportdlg_stat_callbacks,handles,dlg_data,1,...
    ena_array, valuarray, cds_array, posena_valu, fname_list, fpath_list});
set(dlg_data.ui_cl, 'Callback',...
    {@gui_exportdlg_stat_callbacks,handles,dlg_data,0,...
    [],[],[],[],[],[]});
set(dlg_data.ui_fig,'CloseRequestFcn',...
    {@gui_exportdlg_stat_callbacks,handles,dlg_data,0,...
    [],[],[],[],[],[]});





% ------------------------------------------------------------------------
% Export: plots as image file
% ------------------------------------------------------------------------
function ui_exp_fig_as_img_pushbutton_Callback(hObject, eventdata, handles)
% Get selected file index
f_lbx = get(handles.ui_listbox1,'Value');
f_idx = handles.sorted_index(  f_lbx  );
% Get active plots:
pos_active_figs = find( isnan(handles.plots.fig_h)==0 );
if isempty( pos_active_figs )==1
    errordlg('There is no plot active. Please check a plot first!',...
        '! ERROR !');
    return
end
%Num of figures:
num_figs = length( pos_active_figs )*length( f_idx );
%Ask if sure
sure = questdlg(sprintf('This exports an image file %s\n%s %d %s',...
    'of EACH active (checked) plot type for EACH selected file!',...
    'This will result in', num_figs, 'image files. Continue?'),...
    'Continue exporting?','Continue','Cancel','Cancel');
if strcmp(sure,'Continue')~=1
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create dialog
H_fig = 100;
W_fig = 320;
fig_le = handles.gui.pos_screen(3)/2-W_fig/2;
fig_bo = handles.gui.pos_screen(4)/2-H_fig/2;
dlg_data.ui_fig = dialog('Name','Selecting export format...',...
    'Position', [fig_le fig_bo W_fig H_fig]);
H_unit = 0;%needed height for unit - including 5 pixel space on botton
%Select export formats:
dlg_data.ui_file_tx = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 H_fig-H_unit-25 W_fig-20 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String', 'Please select output file format:');
dlg_data.ui_fil_tif = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 H_fig-H_unit-45 W_fig-20 21], 'Style', 'radiobutton',...
    'String', '.tiff image file', 'Value', 1);
dlg_data.ui_fil_fig = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 H_fig-H_unit-65 W_fig-20 21], 'Style', 'radiobutton',...
    'String', '.fig: Matlab figure file', 'Value', 0);
H_file = H_unit+70;%needed height for unit - including 5 pixel space on botton
%Buttons:
dlg_data.ui_ok = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 H_fig-H_file-25 (W_fig-20)/2 21],...
    'Style', 'pushbutton',...
    'String', 'OK', 'Value', 0);%,...
dlg_data.ui_cl = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10+(W_fig-20)/2 H_fig-H_file-25 (W_fig-20)/2 21],...
    'Style', 'pushbutton',...
    'String', 'Cancel', 'Value', 1);%,...
%radiobutton callbacks
set(dlg_data.ui_fil_tif, 'Callback',...
    {@gui_exportdlg_plot_ui_callbacks,dlg_data,dlg_data.ui_fil_tif});
set(dlg_data.ui_fil_fig, 'Callback',...
    {@gui_exportdlg_plot_ui_callbacks,dlg_data,dlg_data.ui_fil_fig});
%Callbacks with dlg_data:
set(dlg_data.ui_ok, 'Callback',...
    {@gui_exportdlg_plot_callbacks,handles,dlg_data,1,pos_active_figs});
set(dlg_data.ui_cl, 'Callback',...
    {@gui_exportdlg_plot_callbacks,handles,dlg_data,0,[]});
set(dlg_data.ui_fig,'CloseRequestFcn',...
    {@gui_exportdlg_plot_callbacks,handles,dlg_data,0,[]});


% ------------------------------------------------------------------------
% Export: variables as mat-file
% ------------------------------------------------------------------------
function ui_exp_ts_pushbutton_Callback(hObject,eventdata,handles)
%Ask if sure
sure = questdlg(sprintf('This exports the error filtered %s %s\n %s\n %s',...
    'time series of each selected files EITHER into one Matlab mat-file ',...
    'OR into ascii-file(s).',...
    '(Velocity units in the exported file: m/s)',...
    'Overwrite check follows!'),...
    'Continue exporting?','Continue','Cancel','Cancel');
if strcmp(sure,'Cancel')==1 || strcmp(sure,'')
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get selected file index
f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
%Check which enabled for export:
f_idx_enabled = [];
list_err_exporting  = {};
for fik=1:length(f_idx)
    %Add data - if fluct seri enabled
    if handles.ts( f_idx(fik) ).seri_enable(1)==1
        f_idx_enabled(end+1,1) = f_idx(fik);
    else
        list_err_exporting{1,end+1} = handles.fullname{ f_idx(fik)};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check validity - if fluct seri not is present for some files:
if isempty( list_err_exporting )==0
    sure = questdlg(['Following files cannot be exported '...
        '(time series not error filtered):',...
        list_err_exporting],...
        'Continue exporting?','Continue','Cancel','Cancel');
    if strcmp(sure,'Cancel')==1 || strcmp(sure,'')
        return
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create dialog
H_fig = 160;
W_fig = 320;
fig_le = handles.gui.pos_screen(3)/2-W_fig/2;
fig_bo = handles.gui.pos_screen(4)/2-H_fig/2;
dlg_data.ui_fig = dialog('Name','Selecting export format...',...
    'Position', [fig_le fig_bo W_fig H_fig]);
%Plot uicontrols
%ITT
H_ofs = 0;%needed height for unit - including 5 pixel space on botton
%Select export format:
dlg_data.ui_file_tx = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 H_fig-H_ofs-25 W_fig-20 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String', 'Please select output file format:');
%.asc
dlg_data.ui_fil_asc = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 30 H_fig-H_ofs-45 W_fig-40 21], 'Style', 'radiobutton',...
    'String', 'each time series exported to a separate ascii-file',...
    'Value', 1);
%.asc options - replace str
dlg_data.ui_repl_tx = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 45 H_fig-H_ofs-70 W_fig-20 21], 'Style', 'text',...
    'HorizontalAlignment','left',...
    'String', 'Replace intermittent invalid samples'' values by:');%,...
dlg_data.ui_repl_999 = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 45 H_fig-H_ofs-85 80 21], 'Style', 'radiobutton',...
    'String', '999', 'Value', 1);%,...
dlg_data.ui_repl_nan = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 105 H_fig-H_ofs-85 80 21], 'Style', 'radiobutton',...
    'String', 'NaN', 'Value', 0);%,...
dlg_data.ui_repl_int = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 165 H_fig-H_ofs-85 120 21], 'Style', 'radiobutton',...
    'String', 'linear interpolation', 'Value', 0);%,...
%.mat
dlg_data.ui_fil_mat = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 30 H_fig-H_ofs-115 W_fig-40 21], 'Style', 'radiobutton',...
    'String', 'one Matlab mat-file is exported with all time series',...
    'Value', 0);
H_ofs = H_ofs+120;%needed height for unit - including 5 pixel space on botton
%Buttons:
dlg_data.ui_ok = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10 H_fig-H_ofs-25 (W_fig-20)/2 21],...
    'Style', 'pushbutton',...
    'String', 'OK', 'Value', 0);%,...
dlg_data.ui_cl = uicontrol('Parent',dlg_data.ui_fig,...
    'HandleVisibility','callback','Units','pixels',...
    'Position',[ 10+(W_fig-20)/2 H_fig-H_ofs-25 (W_fig-20)/2 21],...
    'Style', 'pushbutton',...
    'String', 'Cancel', 'Value', 1);%,...
%radiobutton callbacks
set(dlg_data.ui_repl_999, 'Callback',...
    {@gui_exportdlg_ts_ui_callbacks,dlg_data,dlg_data.ui_repl_999});
set(dlg_data.ui_repl_nan, 'Callback',...
    {@gui_exportdlg_ts_ui_callbacks,dlg_data,dlg_data.ui_repl_nan});
set(dlg_data.ui_repl_int, 'Callback',...
    {@gui_exportdlg_ts_ui_callbacks,dlg_data,dlg_data.ui_repl_int});
set(dlg_data.ui_fil_mat, 'Callback',...
    {@gui_exportdlg_ts_ui_callbacks,dlg_data,dlg_data.ui_fil_mat});
set(dlg_data.ui_fil_asc, 'Callback',...
    {@gui_exportdlg_ts_ui_callbacks,dlg_data,dlg_data.ui_fil_asc});
%Callbacks with dlg_data:
set(dlg_data.ui_ok, 'Callback',...
    {@gui_exportdlg_ts_callbacks,handles,dlg_data,1,f_idx_enabled});
set(dlg_data.ui_cl, 'Callback',...
    {@gui_exportdlg_ts_callbacks,handles,dlg_data,0,[]});
set(dlg_data.ui_fig,'CloseRequestFcn',...
    {@gui_exportdlg_ts_callbacks,handles,dlg_data,0,[]});






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% NOT USED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ui_listbox1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_ts_fs_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_ts_fs_edit_Callback(hObject, eventdata, handles)

function ui_ts_lg_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_ts_lg_edit_Callback(hObject, eventdata, handles)

function ui_prob_type_edit_Callback(hObject, eventdata, handles)

function ui_prob_type_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_prob_nvr_edit_Callback(hObject, eventdata, handles)

function ui_prob_nvr_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_errf_cor_lambda_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_prob_config_edit_Callback(hObject, eventdata, handles)

function ui_prob_config_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_prob_hzvtvr_edit_Callback(hObject, eventdata, handles)

function ui_prob_hzvtvr_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_dure_mean_u_rel_edit_Callback(hObject, eventdata, handles)

function ui_dure_mean_u_rel_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_dure_mean_v_rel_edit_Callback(hObject, eventdata, handles)

function ui_dure_mean_v_rel_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_dure_mean_w_rel_edit_Callback(hObject, eventdata, handles)

function ui_dure_mean_w_rel_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_dure_var_u_rel_edit_Callback(hObject, eventdata, handles)

function ui_dure_var_u_rel_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_dure_var_v_rel_edit_Callback(hObject, eventdata, handles)

function ui_dure_var_v_rel_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_dure_var_w_rel_edit_Callback(hObject, eventdata, handles)

function ui_dure_var_w_rel_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_errf_valpct_edit_Callback(hObject, eventdata, handles)

function ui_errf_valpct_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_errf_valdur_edit_Callback(hObject, eventdata, handles)

function ui_errf_valdur_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_dure_mean_u_edit_Callback(hObject, eventdata, handles)

function ui_dure_mean_u_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_dure_mean_v_edit_Callback(hObject, eventdata, handles)

function ui_dure_mean_v_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_dure_mean_w_edit_Callback(hObject, eventdata, handles)

function ui_dure_mean_w_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_dure_var_u_edit_Callback(hObject, eventdata, handles)

function ui_dure_var_u_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_dure_var_v_edit_Callback(hObject, eventdata, handles)

function ui_dure_var_v_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_dure_var_w_edit_Callback(hObject, eventdata, handles)

function ui_dure_var_w_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_cds_fil_z_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_cds_fil_x_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_cds_fil_y_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_cds_ofs_z_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_cds_ofs_x_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_cds_ofs_y_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_cds_mpt_z_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_cds_mpt_x_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_cds_mpt_y_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_est_u_edit_Callback(hObject, eventdata, handles)

function ui_est_u_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_est_v_edit_Callback(hObject, eventdata, handles)

function ui_est_v_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_est_w_edit_Callback(hObject, eventdata, handles)

function ui_est_w_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_est_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_est_variable_popupmenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_errf_pst_hpfreq_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_errf_crop_1_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ui_errf_crop_2_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 
function chapter_02_main_gui_gets_and_sets
% 
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ------------------------------------------------------------
% Main-GUI/all: ui display values
% ------------------------------------------------------------
function gui_update_ui_values(hObject, handles)
if handles.gui.units_id == 1%mm
    fmt = '%.0f';
elseif handles.gui.units_id == 10%cm
    fmt = '%.1f';
elseif handles.gui.units_id == 1000%m
    fmt = '%.3f';
end
% Depending on number of files loaded:
if isfield(handles, 'sorted_index')==0 || ...
        length( get(handles.ui_listbox1,'Value') )~=1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% No loaded file OR Multiple-file-selection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------------------------------------------------------------
    % Time series file properties:
    % ------------------------------------------------------------
    %Show data_type:
    set(handles.ui_prob_config_edit,'String', 'n. a.' );
    %Show sampling frequency:
    set(handles.ui_ts_fs_edit,       'String', 'n. a.' );
    %Show time-series duration:
    set(handles.ui_ts_lg_edit,       'String', 'n. a.' );
    %Parameters from .hdr file:
    set(handles.ui_prob_type_edit,   'String', 'n. a.' );
    set(handles.ui_prob_nvr_edit,'String', 'n. a.' );
    set(handles.ui_prob_hzvtvr_edit,'String', 'n. a.' );
    %Required duration estimator
    set(handles.ui_est_u_edit,         'String','n. a.');
    set(handles.ui_est_v_edit,         'String','n. a.');
    set(handles.ui_est_w_edit,         'String','n. a.');
    % ------------------------------------------------------------
    % Error filtering - 'on' possible: depending on CONFIG. BUTTON Value
    % ------------------------------------------------------------    
    if get(handles.ui_errf_setdefaults_togglebutton, 'Value')==1
        % - Used components
        set( handles.ui_errf_cmp_u_checkbox,'Value',...
            handles.ts_defaults.p.errfilt.used_velcomp(1)   );
        set( handles.ui_errf_cmp_v_checkbox,'Value',...
            handles.ts_defaults.p.errfilt.used_velcomp(2)   );
        set( handles.ui_errf_cmp_w_checkbox,'Value',...
            handles.ts_defaults.p.errfilt.used_velcomp(3)   );
        % - Crop ts:
        ef_mrk = strcmpi(handles.ts_defaults.p.errfilt.used_types,'crop');
        if sum( ef_mrk )>0
            set( handles.ui_errf_crop_checkbox,           'Value',1);
            set( handles.ui_errf_crop_1_edit,'String', num2str(...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).lambda(1)));
            set( handles.ui_errf_crop_2_edit,'String', num2str(...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).lambda(2)));
        else
            set( handles.ui_errf_crop_checkbox,           'Value',0);
            set( handles.ui_errf_crop_1_edit,'String', 'NaN' );
            set( handles.ui_errf_crop_2_edit,'String', 'NaN' );
        end
        % - Spike filter == PST
        ef_mrk = strcmpi(handles.ts_defaults.p.errfilt.used_types,'pst');
        if sum( ef_mrk )>0
            set( handles.ui_errf_pst_checkbox,           'Value',1);
            set( handles.ui_errf_pst_cmp_u_checkbox,'Value',...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).used_comp(1) );
            set( handles.ui_errf_pst_cmp_v_checkbox,'Value',...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).used_comp(2) );
            set( handles.ui_errf_pst_cmp_w_checkbox,'Value',...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).used_comp(3) );
            if isempty( handles.ts_defaults.p.errfilt.typ( ef_mrk ).hipass_frq )
                set( handles.ui_errf_pst_usehpf_checkbox ,'Value',0);
                set( handles.ui_errf_pst_hpfreq_edit,'String', 'NaN');
            else
                set( handles.ui_errf_pst_usehpf_checkbox, 'Value',1);
                set( handles.ui_errf_pst_hpfreq_edit,'String', num2str(...
                    handles.ts_defaults.p.errfilt.typ( ef_mrk ).hipass_frq));
            end
        else
            set(handles.ui_errf_pst_checkbox,           'Value',0);
            set( handles.ui_errf_pst_cmp_u_checkbox,'Value', 1 );
            set( handles.ui_errf_pst_cmp_v_checkbox,'Value', 1 );
            set( handles.ui_errf_pst_cmp_w_checkbox,'Value', 1 );
            set( handles.ui_errf_pst_usehpf_checkbox,'Value', 0 );
            set( handles.ui_errf_pst_hpfreq_edit,'String', 1 );
        end
        % - Corr threshold init
        set(handles.ui_errf_cor_lambda_edit,'String','NaN');
        % - Filter out samples where any cor is lower than
        ef_mrk=strcmpi(handles.ts_defaults.p.errfilt.used_types,'cormin');
        if sum( ef_mrk )>0
            set(handles.ui_errf_cor_min_checkbox,        'Value',1);
            set(handles.ui_errf_cor_lambda_edit,'String', num2str(...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).lambda));
        else
            set(handles.ui_errf_cor_min_checkbox,        'Value',0)
        end
        % - Filter out samples where the average cor is lower than
        ef_mrk=strcmpi(handles.ts_defaults.p.errfilt.used_types,'coravg');
        if sum( ef_mrk )>0
            set(handles.ui_errf_cor_avg_checkbox,        'Value',1);
            set(handles.ui_errf_cor_lambda_edit,'String', num2str(...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).lambda));
        else
            set(handles.ui_errf_cor_avg_checkbox,        'Value',0);
        end
    else
        set( handles.ui_errf_cmp_u_checkbox,'Value', 0 );
        set( handles.ui_errf_cmp_v_checkbox,'Value', 0 );
        set( handles.ui_errf_cmp_w_checkbox,'Value', 0 );
        set( handles.ui_errf_crop_checkbox, 'Value', 0);
        set( handles.ui_errf_crop_1_edit,   'String', 'NaN' );
        set( handles.ui_errf_crop_2_edit,   'String', 'NaN' );
        set( handles.ui_errf_pst_checkbox,      'Value', 0);
        set( handles.ui_errf_pst_cmp_u_checkbox,'Value', 0 );
        set( handles.ui_errf_pst_cmp_v_checkbox,'Value', 0 );
        set( handles.ui_errf_pst_cmp_w_checkbox,'Value', 0 );
        set( handles.ui_errf_pst_usehpf_checkbox,'Value',0 );
        set( handles.ui_errf_pst_hpfreq_edit,   'String', 'n.a.' );
        set( handles.ui_errf_cor_lambda_edit,   'String', 'n. a.');
        set( handles.ui_errf_cor_min_checkbox,  'Value', 0)
        set( handles.ui_errf_cor_avg_checkbox,  'Value', 0);
    end
    % Validity after err.filt:
    set(handles.ui_errf_valpct_edit,'String',' n. a.' );    
    set(handles.ui_errf_valdur_edit,'String',' n. a.' );    
    % ------------------------------------------------------------
    % Coordinates:
    % ------------------------------------------------------------
    %%%% Show file coordinates:
    set(handles.ui_cds_fil_x_edit,'String','n. a.');
    set(handles.ui_cds_fil_y_edit,'String','n. a.');
    set(handles.ui_cds_fil_z_edit,'String','n. a.');
    %%%% Show measpt_ofs coordinates:
    set(handles.ui_cds_ofs_x_edit,'String','n. a.');
    set(handles.ui_cds_ofs_y_edit,'String','n. a.');
    set(handles.ui_cds_ofs_z_edit,'String','n. a.');
    %%%% Show measpt coordinates:
    set(handles.ui_cds_mpt_x_edit,'String','n. a.');
    set(handles.ui_cds_mpt_y_edit,'String','n. a.');
    set(handles.ui_cds_mpt_z_edit,'String','n. a.');
    % ------------------------------------------------------------
    % Duration error:
    % ------------------------------------------------------------
    set(handles.ui_dure_mean_u_edit,'String',' n. a.' );
    set(handles.ui_dure_mean_v_edit,'String',' n. a.' );
    set(handles.ui_dure_mean_w_edit,'String',' n. a.' );
    set(handles.ui_dure_var_u_edit, 'String',' n. a.' );
    set(handles.ui_dure_var_v_edit, 'String',' n. a.' );
    set(handles.ui_dure_var_w_edit, 'String',' n. a.' );
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Single-file-selection: length(f_idx)==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get selected file index
    f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
    % ------------------------------------------------------------
    % Time series file properties:
    % ------------------------------------------------------------
    %Show probe config:
    set(handles.ui_prob_config_edit,'String',...
        handles.ts( f_idx ).p.dat_props.probe_cfg );
    %Show sampling frequency:
    set(handles.ui_ts_fs_edit,'String',...
        sprintf('%d Hz',  handles.ts( f_idx ).p.dat_props.frq));
    %Show velocity range:
    if strcmp( handles.ts( f_idx ).p.dat_props.probe_nvr{1},'0')~=1
        %Nominal velocity range
        set(handles.ui_prob_nvr_edit,'String', sprintf( '%s %s',...
            handles.ts( f_idx ).p.dat_props.probe_nvr{1},...
            handles.ts( f_idx ).p.dat_props.probe_nvr{2} ));
        %Range for components parallel vs. normal to transmitter beam
        if strcmp( handles.ts( f_idx ).p.dat_props.type, 'vno')
            velrng_mrk = handles.viper_props.velrange_limits_vno(:,1)==...
                str2double( handles.ts( f_idx ).p.dat_props.probe_nvr{1} );
            set(handles.ui_prob_hzvtvr_edit,'String',...
                sprintf( '%1.2f || %1.2f',...
                handles.viper_props.velrange_limits_vno( velrng_mrk ,2),...
                handles.viper_props.velrange_limits_vno( velrng_mrk ,3) ));
        end
        if strcmp( handles.ts( f_idx ).p.dat_props.type, 'vec')
            velrng_mrk = handles.viper_props.velrange_limits_vec(:,1)==...
                str2double( handles.ts( f_idx ).p.dat_props.probe_nvr{1} );
            set(handles.ui_prob_hzvtvr_edit,'String',...
                sprintf( '%1.2f || %1.2f',...
                handles.viper_props.velrange_limits_vec( velrng_mrk ,2),...
                handles.viper_props.velrange_limits_vec( velrng_mrk ,3) ));
        end
        %Text:
        comp_names = {'u' 'v' 'w'};
        comp_n = ...
            comp_names(strcmp( handles.ts( f_idx ).intl.cmp_probe_np,'n'));
        comp_p = ...
            comp_names(strcmp( handles.ts( f_idx ).intl.cmp_probe_np,'p'));
        set(handles.ui_prob_text2,'String',...
            sprintf( 'ADV Vel. Ranges\n NVR &  %s,%s || %s',...
            comp_n{1}, comp_n{2}, comp_p{1} ));
    else
        set(handles.ui_prob_nvr_edit,'String', 'n. a.' );
        set(handles.ui_prob_hzvtvr_edit,'String','n. a.' );
        set(handles.ui_prob_text2,'String',...
            sprintf( 'ADV Vel. Ranges\n NVR &  u,v || w'));
    end
    %Show time-series duration:
    set(handles.ui_ts_lg_edit,'String',...
        sprintf('%5.0f s', length( handles.ts(f_idx).raw_t )/...
                handles.ts(f_idx).p.dat_props.frq) );%/...
    %data_type: %short: handles.ts( f_idx ).p.dat_props.type
    if strcmp( handles.ts( f_idx ).p.dat_props.type, 'vno')
        set(handles.ui_prob_type_edit,'String','Vectrino');
    elseif strcmp( handles.ts( f_idx ).p.dat_props.type, 'vec')
        set(handles.ui_prob_type_edit,'String','Vector');
    else
        set(handles.ui_prob_type_edit,'String','User def.');
    end
    %Required duration estimator
    gui_update_estimator_display( f_idx, handles )
    % ------------------------------------------------------------
    % Error filtering - depending on CONFIG. BUTTON Value
    % ------------------------------------------------------------
    if get(handles.ui_errf_setdefaults_togglebutton, 'Value')==1%on
        % - Used components
        set( handles.ui_errf_cmp_u_checkbox,'Value',...
            handles.ts_defaults.p.errfilt.used_velcomp(1)   );
        set( handles.ui_errf_cmp_v_checkbox,'Value',...
            handles.ts_defaults.p.errfilt.used_velcomp(2)   );
        set( handles.ui_errf_cmp_w_checkbox,'Value',...
            handles.ts_defaults.p.errfilt.used_velcomp(3)   );
        % - Crop ts:
        ef_mrk = strcmpi(handles.ts_defaults.p.errfilt.used_types,'crop');
        if sum( ef_mrk )>0
            set( handles.ui_errf_crop_checkbox,           'Value',1);
            set( handles.ui_errf_crop_1_edit,'String', num2str(...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).lambda(1)));
            set( handles.ui_errf_crop_2_edit,'String', num2str(...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).lambda(2)));
        else
            set( handles.ui_errf_crop_checkbox,           'Value',0);
            set( handles.ui_errf_crop_1_edit,'String', 'NaN' );
            set( handles.ui_errf_crop_2_edit,'String', 'NaN' );
        end
        % - Spike filter == PST
        ef_mrk = strcmpi(handles.ts_defaults.p.errfilt.used_types,'pst');
        if sum( ef_mrk )>0
            set( handles.ui_errf_pst_checkbox,           'Value',1);
            set( handles.ui_errf_pst_cmp_u_checkbox,'Value',...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).used_comp(1) );
            set( handles.ui_errf_pst_cmp_v_checkbox,'Value',...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).used_comp(2) );
            set( handles.ui_errf_pst_cmp_w_checkbox,'Value',...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).used_comp(3) );
            if isempty( handles.ts_defaults.p.errfilt.typ( ef_mrk ).hipass_frq )
                set( handles.ui_errf_pst_usehpf_checkbox ,'Value',0);
                set( handles.ui_errf_pst_hpfreq_edit,'String', 'NaN');
            else
                set( handles.ui_errf_pst_usehpf_checkbox, 'Value',1);
                set( handles.ui_errf_pst_hpfreq_edit,'String', num2str(...
                    handles.ts_defaults.p.errfilt.typ( ef_mrk ).hipass_frq));
            end
        else
            set(handles.ui_errf_pst_checkbox,           'Value',0);
            set( handles.ui_errf_pst_cmp_u_checkbox,'Value', 1 );
            set( handles.ui_errf_pst_cmp_v_checkbox,'Value', 1 );
            set( handles.ui_errf_pst_cmp_w_checkbox,'Value', 1 );
            set( handles.ui_errf_pst_usehpf_checkbox,'Value', 0 );
            set( handles.ui_errf_pst_hpfreq_edit,'String', 1 );
        end
        % - Corr threshold init
        set(handles.ui_errf_cor_lambda_edit,'String','NaN');
        % - Filter out samples where any cor is lower than
        ef_mrk=strcmpi(handles.ts_defaults.p.errfilt.used_types,'cormin');
        if sum( ef_mrk )>0
            set(handles.ui_errf_cor_min_checkbox,        'Value',1);
            set(handles.ui_errf_cor_lambda_edit,'String', num2str(...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).lambda));
        else
            set(handles.ui_errf_cor_min_checkbox,        'Value',0)
        end
        % - Filter out samples where the average cor is lower than
        ef_mrk=strcmpi(handles.ts_defaults.p.errfilt.used_types,'coravg');
        if sum( ef_mrk )>0
            set(handles.ui_errf_cor_avg_checkbox,        'Value',1);
            set(handles.ui_errf_cor_lambda_edit,'String', num2str(...
                handles.ts_defaults.p.errfilt.typ( ef_mrk ).lambda));
        else
            set(handles.ui_errf_cor_avg_checkbox,        'Value',0);
        end
    else%setdefaults==off -> show file specific values
        % - Used components
        set( handles.ui_errf_cmp_u_checkbox,'Value',...
            handles.ts( f_idx ).p.errfilt.used_velcomp(1)   );
        set( handles.ui_errf_cmp_v_checkbox,'Value',...
            handles.ts( f_idx ).p.errfilt.used_velcomp(2)   );
        set( handles.ui_errf_cmp_w_checkbox,'Value',...
            handles.ts( f_idx ).p.errfilt.used_velcomp(3)   );
        % - Crop ts:
        ef_mrk = strcmpi(handles.ts( f_idx ).p.errfilt.used_types,'crop');
        if sum( ef_mrk )>0
            set( handles.ui_errf_crop_checkbox,           'Value',1);
            set( handles.ui_errf_crop_1_edit,'String', num2str(...
                handles.ts( f_idx ).p.errfilt.typ( ef_mrk ).lambda(1)));
            set( handles.ui_errf_crop_2_edit,'String', num2str(...
                handles.ts( f_idx ).p.errfilt.typ( ef_mrk ).lambda(2)));
        else%file specific crop==off
            set( handles.ui_errf_crop_checkbox,           'Value',0);
            set( handles.ui_errf_crop_1_edit,'String', 'NaN' );
            set( handles.ui_errf_crop_2_edit,'String', 'NaN' );
        end
        % - Spike filter == PST
        ef_mrk = strcmpi(handles.ts( f_idx ).p.errfilt.used_types,'pst');
        if sum( ef_mrk )>0%ON for file
            set( handles.ui_errf_pst_checkbox,           'Value',1);
            set( handles.ui_errf_pst_cmp_u_checkbox,'Value',...
                handles.ts( f_idx ).p.errfilt.typ( ef_mrk ).used_comp(1)   );
            set( handles.ui_errf_pst_cmp_v_checkbox,'Value',...
                handles.ts( f_idx ).p.errfilt.typ( ef_mrk ).used_comp(2)   );
            set( handles.ui_errf_pst_cmp_w_checkbox,'Value',...
                handles.ts( f_idx ).p.errfilt.typ( ef_mrk ).used_comp(3)   );
            if isempty( handles.ts( f_idx ).p.errfilt.typ( ef_mrk ).hipass_frq )
                set( handles.ui_errf_pst_usehpf_checkbox ,'Value',0);
                set( handles.ui_errf_pst_hpfreq_edit,'String', 'NaN');
            else
                set( handles.ui_errf_pst_usehpf_checkbox, 'Value',1);
                set( handles.ui_errf_pst_hpfreq_edit,'String', num2str(...
                    handles.ts( f_idx ).p.errfilt.typ( ef_mrk ).hipass_frq));
            end
        else%file specific pst==off
            set(handles.ui_errf_pst_checkbox,           'Value',0);
            set( handles.ui_errf_pst_cmp_u_checkbox,'Value', 1 );
            set( handles.ui_errf_pst_cmp_v_checkbox,'Value', 1 );
            set( handles.ui_errf_pst_cmp_w_checkbox,'Value', 1 );
            set( handles.ui_errf_pst_usehpf_checkbox,'Value', 0 );
            set( handles.ui_errf_pst_hpfreq_edit,'String', 'NaN' );
        end
    end
    % - Corr threshold init
    set(handles.ui_errf_cor_lambda_edit,'String','NaN');
    % Cormin: Filter out samples where any cor is lower than
    ef_mrk = strcmpi(handles.ts( f_idx ).p.errfilt.used_types,'cormin');
    if sum( ef_mrk )>0
        set(handles.ui_errf_cor_min_checkbox,        'Value',1);
        set(handles.ui_errf_cor_lambda_edit,'String', num2str(...
            handles.ts( f_idx ).p.errfilt.typ( ef_mrk ).lambda));
    else
        set(handles.ui_errf_cor_min_checkbox,        'Value',0)
    end
    % Coravg: Filter out samples where the average cor is lower than
    ef_mrk = strcmpi(handles.ts( f_idx ).p.errfilt.used_types,'coravg');
    if sum( ef_mrk )>0
        set(handles.ui_errf_cor_avg_checkbox,        'Value',1);
        set(handles.ui_errf_cor_lambda_edit,'String', num2str(...
            handles.ts( f_idx ).p.errfilt.typ( ef_mrk ).lambda));
    else
        set(handles.ui_errf_cor_avg_checkbox,        'Value',0);
    end
    % Validity after err.filt:
    set(handles.ui_errf_valpct_edit,'String',...
        sprintf('%d',  round(handles.ts( f_idx ).valu_data( 16 ))) );
    set(handles.ui_errf_valdur_edit,'String',...
        sprintf('%0.1f',  round(handles.ts( f_idx ).valu_data( 24 ))) );
    % ------------------------------------------------------------
    % Coordinates:
    % ------------------------------------------------------------
    %%%% Show file coordinates:
    set(handles.ui_cds_fil_x_edit,'String',...
        sprintf( fmt, handles.ts( f_idx ).intl.coordsys.file_coords(1)));
    set(handles.ui_cds_fil_y_edit,'String',...
        sprintf( fmt, handles.ts( f_idx ).intl.coordsys.file_coords(2)));
    set(handles.ui_cds_fil_z_edit,'String',...
        sprintf( fmt, handles.ts( f_idx ).intl.coordsys.file_coords(3)));
    %%%% Show measpt_ofs coordinates:
    set(handles.ui_cds_ofs_x_edit,'String',...
        sprintf( fmt, handles.ts( f_idx ).intl.coordsys.measpt_offset(1)));
    set(handles.ui_cds_ofs_y_edit,'String',...
        sprintf( fmt, handles.ts( f_idx ).intl.coordsys.measpt_offset(2)));
    set(handles.ui_cds_ofs_z_edit,'String',...
        sprintf( fmt, handles.ts( f_idx ).intl.coordsys.measpt_offset(3)));
    %%%% Show measpt coordinates:
    set(handles.ui_cds_mpt_x_edit,'String',...
        sprintf( fmt, handles.ts( f_idx ).intl.coordsys.measpt_coords(1)));
    set(handles.ui_cds_mpt_y_edit,'String',...
        sprintf( fmt, handles.ts( f_idx ).intl.coordsys.measpt_coords(2)));
    set(handles.ui_cds_mpt_z_edit,'String',...
        sprintf( fmt, handles.ts( f_idx ).intl.coordsys.measpt_coords(3)));
    %%%% Show Direction of Probe X:
    set(handles.ui_cds_dir_posx_radiobutton,'Value',0);
    set(handles.ui_cds_dir_negx_radiobutton,'Value',0);
    set(handles.ui_cds_dir_posy_radiobutton,'Value',0);
    set(handles.ui_cds_dir_negy_radiobutton,'Value',0);
    if strcmp(handles.ts( f_idx ).intl.coordsys.dir_probe_u,'-x')
        set(handles.ui_cds_dir_negx_radiobutton,'Value',1);
    elseif strcmp(handles.ts( f_idx ).intl.coordsys.dir_probe_u,'+y')
        set(handles.ui_cds_dir_posy_radiobutton,'Value',1);
    elseif strcmp(handles.ts( f_idx ).intl.coordsys.dir_probe_u,'-y')
        set(handles.ui_cds_dir_negy_radiobutton,'Value',1);
    else%default: handles.ts( f_idx ).intl.coordsys.dir_probe_u=='+x'
        set(handles.ui_cds_dir_posx_radiobutton,'Value',1);
    end
    % ------------------------------------------------------------
    % Duration error:
    % ------------------------------------------------------------
    if get(handles.ui_dure_rel_radiobutton,'Value')==1%percent error
        %Editboxes - mean:
        set(handles.ui_dure_mean_u_edit,'String',...
            sprintf('%d',  round(handles.ts( f_idx ).valu_data( 34 ))) );
        set(handles.ui_dure_mean_v_edit,'String',...
            sprintf('%d',  round(handles.ts( f_idx ).valu_data( 35 ))) );
        set(handles.ui_dure_mean_w_edit,'String',...
            sprintf('%d',  round(handles.ts( f_idx ).valu_data( 36 ))) );
        set(handles.ui_dure_var_u_edit, 'String',...
            sprintf('%d',  round(handles.ts( f_idx ).valu_data( 40 ))) );
        set(handles.ui_dure_var_v_edit, 'String',...
            sprintf('%d',  round(handles.ts( f_idx ).valu_data( 41 ))) );
        set(handles.ui_dure_var_w_edit, 'String',...
            sprintf('%d',  round(handles.ts( f_idx ).valu_data( 42 ))) );
        %Unitboxes:
        set( handles.ui_dure_mean_unit_text, 'String', '%')
        set( handles.ui_dure_var_unit_text,  'String', '%')
    else%absolute error
        %Editboxes:
        set(handles.ui_dure_mean_u_edit,'String',...
            sprintf( fmt,      handles.ts( f_idx ).valu_data( 31 )) );
        set(handles.ui_dure_mean_v_edit,'String',...
            sprintf( fmt,      handles.ts( f_idx ).valu_data( 32 )) );
        set(handles.ui_dure_mean_w_edit,'String',...
            sprintf( fmt,      handles.ts( f_idx ).valu_data( 33 )) );
        set(handles.ui_dure_var_u_edit, 'String',...
            sprintf( fmt,      handles.ts( f_idx ).valu_data( 37 )) );
        set(handles.ui_dure_var_v_edit, 'String',...
            sprintf( fmt,      handles.ts( f_idx ).valu_data( 38 )) );
        set(handles.ui_dure_var_w_edit, 'String',...
            sprintf( fmt,      handles.ts( f_idx ).valu_data( 39 )) );
        %Unitboxes:
        if handles.gui.units_id == 10%== cm, cm/s
            set( handles.ui_dure_mean_unit_text, 'String', 'cm/s')
            set( handles.ui_dure_var_unit_text, 'String', 'cm2/s2')
        elseif handles.gui.units_id == 1000%== m, m/s
            set( handles.ui_dure_mean_unit_text, 'String', 'm/s')
            set( handles.ui_dure_var_unit_text, 'String', 'm2/s2')
        end
    end
end
%DVL:
if isfield( handles, 'dvl')==1
    [ ~, ~ ] = viper_dvl('dvl_update_ui_values', handles, { f_idx });
end






% ------------------------------------------------------------
% Main-GUI/all: ui validity (enable)
% ------------------------------------------------------------
function gui_update_ui_enables(hObject, handles)
%Depending on on|off state of Err.filt-EDIT & Err.filt-config & CDS-Edit
if get(handles.ui_errf_edit_n_cancel_toggle_b, 'Value')==0 &&...
        get(handles.ui_errf_setdefaults_togglebutton, 'Value')==0 && ...
        get(handles.ui_cds_edit_n_cancel_togglebutton, 'Value')==0
    % If Err.filt-EDIT and CDS-Edit off
    % ------------------------------------------------------------
    % Depending on number of files loaded:
    if isfield(handles, 'sorted_index')==0
        %No files loaded
        %%%% Main group:
        set(handles.ui__main_addfile_pushbutton,'Enable','On');
        set(handles.ui_listbox1,'Enable','Off');
        %%%% Plot group
        set(handles.ui_plot_fig1_checkbox,      'Enable','Off');
        set(handles.ui_plot_fig2_checkbox,      'Enable','Off');
        set(handles.ui_plot_fig3_checkbox,      'Enable','Off');
        set(handles.ui_plot_fig4_checkbox,      'Enable','Off');
        set(handles.ui_plot_fig5_checkbox,      'Enable','Off');
        set(handles.ui_plot_fig6_checkbox,      'Enable','Off');
        set(handles.ui_plot_fig9_checkbox,      'Enable','Off');
        set(handles.ui_plot_comp_u_checkbox,    'Enable','Off');
        set(handles.ui_plot_comp_v_checkbox,    'Enable','Off');
        set(handles.ui_plot_comp_w_checkbox,    'Enable','Off');
        set(handles.ui_plot_lims_set_togglebutton,'Enable','Off');
        set(handles.ui_plot_config_pushbutton,  'Enable','On');
        set(handles.ui_plot_add53tofig_togglebutton,'Enable','Off');
        set(handles.ui_plot_add3tofig_togglebutton,'Enable','Off');
        set(handles.ui_plot_add1tofig_togglebutton,'Enable','Off');
        %%%% Duration calculator
        set(handles.ui_est_variable_popupmenu,   'Enable','Off');
        set(handles.ui_est_edit,                 'Enable','Off');
        %%%% Err.filt group:
        set(handles.ui_errf_edit_n_cancel_toggle_b,'Enable','Off');
        set(handles.ui_errf_apply_pushbutton,      'Enable','Off');
        set(handles.ui_errf_reload_par_pushbutton, 'Enable','Off');
        set(handles.ui_errf_setdefaults_togglebutton, 'Enable','On');
        %
        set( handles.ui_errf_crop_checkbox,'Enable','Off');
        set( handles.ui_errf_crop_1_edit,  'Enable','Off');
        set( handles.ui_errf_crop_2_edit,  'Enable','Off');
        %
        set(handles.ui_errf_cmp_u_checkbox,     'Enable','Off');
        set(handles.ui_errf_cmp_v_checkbox,     'Enable','Off');
        set(handles.ui_errf_cmp_w_checkbox,     'Enable','Off');
        set(handles.ui_errf_pst_checkbox,       'Enable','Off');
        set(handles.ui_errf_pst_text,           'Enable','Off');
        set(handles.ui_errf_pst_cmp_u_checkbox, 'Enable','Off');
        set(handles.ui_errf_pst_cmp_v_checkbox, 'Enable','Off');
        set(handles.ui_errf_pst_cmp_w_checkbox, 'Enable','Off');
        set(handles.ui_errf_pst_usehpf_checkbox,'Enable','Off');
        set(handles.ui_errf_pst_hpfreq_edit,    'Enable','Off');
        set(handles.ui_errf_pst_hpfreq_txt,     'Enable','Off');
        set(handles.ui_errf_cor_text,           'Enable','Off');
        set(handles.ui_errf_cor_text_sign,      'Enable','Off');
        set(handles.ui_errf_cor_min_checkbox,   'Enable','Off');
        set(handles.ui_errf_cor_avg_checkbox,   'Enable','Off');
        set(handles.ui_errf_cor_lambda_edit,    'Enable','Off');
        %%%% Coords group
        set(handles.ui_cds_edit_n_cancel_togglebutton, 'Enable','Off');
        set(handles.ui_cds_apply_pushbutton,       'Enable','Off');
        set(handles.ui_cds_getfromfile_pushbutton, 'Enable','Off');
        %
        set(handles.ui_cds_fil_x_edit,'Enable','Off');
        set(handles.ui_cds_fil_y_edit,'Enable','Off');
        set(handles.ui_cds_fil_z_edit,'Enable','Off');
        set(handles.ui_cds_ofs_x_edit,'Enable','Off');
        set(handles.ui_cds_ofs_y_edit,'Enable','Off');
        set(handles.ui_cds_ofs_z_edit,'Enable','Off');
        set(handles.ui_cds_mpt_x_edit,'Enable','Off');
        set(handles.ui_cds_mpt_y_edit,'Enable','Off');
        set(handles.ui_cds_mpt_z_edit,'Enable','Off');
        set(handles.ui_cds_dir_posx_radiobutton,'Enable','Off');
        set(handles.ui_cds_dir_negx_radiobutton,'Enable','Off');
        set(handles.ui_cds_dir_posy_radiobutton,'Enable','Off');
        set(handles.ui_cds_dir_negy_radiobutton,'Enable','Off');
        %%%% Error estimator
        set(handles.ui_dure_rel_radiobutton,      'Enable','Off');
        set(handles.ui_dure_abs_radiobutton,      'Enable','Off');
        %%%% Others
        %Export group
        set(handles.ui_exp_pars_pushbutton, 'Enable','Off');
        set(handles.ui_exp_stat_pushbutton,     'Enable','Off');
        set(handles.ui_exp_fig_as_img_pushbutton, 'Enable','Off');
        set(handles.ui_exp_ts_pushbutton,  'Enable','Off');
    else%Files are loaded:
        %%%% Main group:
        set(handles.ui__main_addfile_pushbutton,'Enable','On');
        set(handles.ui_listbox1,'Enable','On');
        %%%% Plot group
        set(handles.ui_plot_fig1_checkbox,      'Enable','On');
        set(handles.ui_plot_fig2_checkbox,      'Enable','On');
        set(handles.ui_plot_fig3_checkbox,      'Enable','On');
        set(handles.ui_plot_fig4_checkbox,      'Enable','On');
        set(handles.ui_plot_fig5_checkbox,      'Enable','On');
        set(handles.ui_plot_fig6_checkbox,      'Enable','On');
        set(handles.ui_plot_fig9_checkbox,      'Enable','On');
        set(handles.ui_plot_comp_u_checkbox,'Enable','On');
        set(handles.ui_plot_comp_v_checkbox,'Enable','On');
        set(handles.ui_plot_comp_w_checkbox,'Enable','On');
        set(handles.ui_plot_lims_set_togglebutton,'Enable','On');
        set(handles.ui_plot_config_pushbutton,'Enable','On');
        %%%% Err.filt group:
        set(handles.ui_errf_apply_pushbutton,      'Enable','Off');
        set(handles.ui_errf_reload_par_pushbutton, 'Enable','Off');
        set(handles.ui_errf_setdefaults_togglebutton, 'Enable','On');
        %
        set( handles.ui_errf_crop_checkbox,'Enable','Off');
        set( handles.ui_errf_crop_1_edit,  'Enable','Off');
        set( handles.ui_errf_crop_2_edit,  'Enable','Off');
        %
        set(handles.ui_errf_cmp_u_checkbox,     'Enable','Off');
        set(handles.ui_errf_cmp_v_checkbox,     'Enable','Off');
        set(handles.ui_errf_cmp_w_checkbox,     'Enable','Off');
        set(handles.ui_errf_pst_checkbox,       'Enable','Off');
        set(handles.ui_errf_pst_text,           'Enable','Off');
        set(handles.ui_errf_pst_cmp_u_checkbox, 'Enable','Off');
        set(handles.ui_errf_pst_cmp_v_checkbox, 'Enable','Off');
        set(handles.ui_errf_pst_cmp_w_checkbox, 'Enable','Off');
        set(handles.ui_errf_pst_usehpf_checkbox,'Enable','Off');
        set(handles.ui_errf_pst_hpfreq_edit,    'Enable','Off');
        set(handles.ui_errf_pst_hpfreq_txt,     'Enable','Off');
        set(handles.ui_errf_cor_text,           'Enable','Off');
        set(handles.ui_errf_cor_text_sign,      'Enable','Off');
        set(handles.ui_errf_cor_min_checkbox,   'Enable','Off');
        set(handles.ui_errf_cor_avg_checkbox,   'Enable','Off');
        set(handles.ui_errf_cor_lambda_edit,    'Enable','Off');
        %%%% Coords group
        set(handles.ui_cds_edit_n_cancel_togglebutton, 'Enable','On');
        set(handles.ui_cds_apply_pushbutton,       'Enable','Off');
        set(handles.ui_cds_getfromfile_pushbutton, 'Enable','On');
        %
        set(handles.ui_cds_fil_x_edit,'Enable','Off');
        set(handles.ui_cds_fil_y_edit,'Enable','Off');
        set(handles.ui_cds_fil_z_edit,'Enable','Off');
        set(handles.ui_cds_ofs_x_edit,'Enable','Off');
        set(handles.ui_cds_ofs_y_edit,'Enable','Off');
        set(handles.ui_cds_ofs_z_edit,'Enable','Off');
        set(handles.ui_cds_mpt_x_edit,'Enable','Off');
        set(handles.ui_cds_mpt_y_edit,'Enable','Off');
        set(handles.ui_cds_mpt_z_edit,'Enable','Off');
        set(handles.ui_cds_dir_posx_radiobutton,'Enable','Off');
        set(handles.ui_cds_dir_negx_radiobutton,'Enable','Off');
        set(handles.ui_cds_dir_posy_radiobutton,'Enable','Off');
        set(handles.ui_cds_dir_negy_radiobutton,'Enable','Off');
        %%%% Error estimator
        set(handles.ui_dure_rel_radiobutton,      'Enable','On');
        set(handles.ui_dure_abs_radiobutton,      'Enable','On');
        %%%% Others
        %Export group
        set(handles.ui_exp_pars_pushbutton, 'Enable','On');
        set(handles.ui_exp_stat_pushbutton,     'Enable','On');
        set(handles.ui_exp_fig_as_img_pushbutton, 'Enable','On');
        set(handles.ui_exp_ts_pushbutton,  'Enable','On');
        %Lisbox:
        set(handles.ui__main_addfile_pushbutton,'Enable','On');
        set(handles.ui_listbox1,'Enable','On');
        %Depending on number of selected files:
        if numel(handles.sorted_index(get(handles.ui_listbox1,'Value')))==1
            %Files loaded & one file selected
            %%%% Plot group
            set(handles.ui_plot_add53tofig_togglebutton,    'Enable','On');
            set(handles.ui_plot_add3tofig_togglebutton,     'Enable','On');
            set(handles.ui_plot_add1tofig_togglebutton,     'Enable','On');
            %%%% Duration calculator
            set(handles.ui_est_variable_popupmenu,   'Enable','On');
            set(handles.ui_est_edit,                 'Enable','On');
            %%%% Err.filt group:
            set(handles.ui_errf_edit_n_cancel_toggle_b, 'Enable','On')
        else
            %Files loaded & multiple files selected
            %%%% Plot group
            set(handles.ui_plot_add53tofig_togglebutton,   'Enable','Off');
            set(handles.ui_plot_add3tofig_togglebutton,   'Enable','Off');
            set(handles.ui_plot_add1tofig_togglebutton,   'Enable','Off');
            %%%% Duration calculator
            set(handles.ui_est_variable_popupmenu,   'Enable','Off');
            set(handles.ui_est_edit,                 'Enable','Off');
            %%%% Err.filt group:
            set(handles.ui_errf_edit_n_cancel_toggle_b,    'Enable','Off');
        end
    end
else% If Err.filt-EDIT Err.filt-config or CDS-Edit on
    %(these can be on ONLY after at least loaded file!)
    %%%% Main group:
    set(handles.ui__main_addfile_pushbutton,'Enable','Off');
    set(handles.ui_listbox1,'Enable','Off');
    %%%% Plot group
    set(handles.ui_plot_fig1_checkbox,      'Enable','Off');
    set(handles.ui_plot_fig2_checkbox,      'Enable','Off');
    set(handles.ui_plot_fig3_checkbox,      'Enable','Off');
    set(handles.ui_plot_fig4_checkbox,      'Enable','Off');
    set(handles.ui_plot_fig5_checkbox,      'Enable','Off');
    set(handles.ui_plot_fig6_checkbox,      'Enable','Off');
    set(handles.ui_plot_fig9_checkbox,      'Enable','Off');
    set(handles.ui_plot_comp_u_checkbox,    'Enable','Off');
    set(handles.ui_plot_comp_v_checkbox,    'Enable','Off');
    set(handles.ui_plot_comp_w_checkbox,    'Enable','Off');
    set(handles.ui_plot_lims_set_togglebutton,'Enable','Off');
    set(handles.ui_plot_config_pushbutton,'Enable','Off');
    set(handles.ui_plot_add53tofig_togglebutton,'Enable','Off');
    set(handles.ui_plot_add3tofig_togglebutton,'Enable','Off');
    set(handles.ui_plot_add1tofig_togglebutton,'Enable','Off');
    %%%% Duration calculator
    set(handles.ui_est_variable_popupmenu,   'Enable','Off');
    set(handles.ui_est_edit,                 'Enable','Off');
    %%%% Err. filt.: most LATER
    set(handles.ui_errf_edit_n_cancel_toggle_b,'Enable','Off');
    set(handles.ui_errf_setdefaults_togglebutton,'Enable','Off');
    %%%% Coords group: most LATER
    set(handles.ui_cds_edit_n_cancel_togglebutton,'Enable','Off');
    set(handles.ui_cds_getfromfile_pushbutton, 'Enable','Off');
    %%%% Error estimator
    set(handles.ui_dure_rel_radiobutton,      'Enable','Off');
    set(handles.ui_dure_abs_radiobutton,      'Enable','Off');
    %%%% Others
    %Export group
    set(handles.ui_exp_pars_pushbutton, 'Enable','Off');
    set(handles.ui_exp_stat_pushbutton,     'Enable','Off');
    set(handles.ui_exp_fig_as_img_pushbutton, 'Enable','Off');
    set(handles.ui_exp_ts_pushbutton,  'Enable','Off');
    % ------------------------------------------------------------
    %Depending on which button on:
    %%%% Error filtering set defaults - no file needed:
    if get(handles.ui_errf_setdefaults_togglebutton, 'Value')==1
        set(handles.ui_errf_setdefaults_togglebutton,'Enable','On');
        set(handles.ui_errf_reload_par_pushbutton,   'Enable','Off');
        set(handles.ui_errf_apply_pushbutton,        'Enable','On');
        %
        set(handles.ui_errf_cmp_u_checkbox,     'Enable','On');
        set(handles.ui_errf_cmp_v_checkbox,     'Enable','On');
        set(handles.ui_errf_cmp_w_checkbox,     'Enable','On');
        set(handles.ui_errf_crop_checkbox,      'Enable','On');
        set(handles.ui_errf_pst_checkbox,       'Enable','On');
        %
        set(handles.ui_errf_cor_text,           'Enable','On');
        set(handles.ui_errf_cor_text_sign,      'Enable','On');
        set(handles.ui_errf_cor_min_checkbox,   'Enable','On');
        set(handles.ui_errf_cor_avg_checkbox,   'Enable','On');
        set(handles.ui_errf_cor_lambda_edit,    'Enable','On');
        %
        if get( handles.ui_errf_crop_checkbox, 'Value')==1
            set( handles.ui_errf_crop_1_edit,  'Enable','On');
            set( handles.ui_errf_crop_2_edit,  'Enable','On');
        else
            set( handles.ui_errf_crop_1_edit,  'Enable','Off');
            set( handles.ui_errf_crop_2_edit,  'Enable','Off');
        end
        %
        if get(handles.ui_errf_cmp_u_checkbox, 'Value')==0
            %%% pst
            set(handles.ui_errf_pst_cmp_u_checkbox, 'Value',0)
            set(handles.ui_errf_pst_cmp_u_checkbox, 'Enable','Off')
        else
            set(handles.ui_errf_pst_cmp_u_checkbox, 'Enable','On')
        end
        if get(handles.ui_errf_cmp_v_checkbox, 'Value')==0
            %%% pst
            set(handles.ui_errf_pst_cmp_v_checkbox, 'Value',0)
            set(handles.ui_errf_pst_cmp_v_checkbox, 'Enable','Off')
        else
            set(handles.ui_errf_pst_cmp_v_checkbox, 'Enable','On')
        end
        if get(handles.ui_errf_cmp_w_checkbox, 'Value')==0
            %%% pst
            set(handles.ui_errf_pst_cmp_w_checkbox, 'Value',0)
            set(handles.ui_errf_pst_cmp_w_checkbox, 'Enable','Off')
        else
            set(handles.ui_errf_pst_cmp_w_checkbox, 'Enable','On')
        end
        %
        if get(handles.ui_errf_pst_checkbox, 'Value')==1
            set(handles.ui_errf_pst_text,           'Enable','On')
            if get(handles.ui_errf_cmp_u_checkbox, 'Value')==0
                set(handles.ui_errf_pst_cmp_u_checkbox, 'Enable','Off')
            else
                set(handles.ui_errf_pst_cmp_u_checkbox, 'Enable','On')
            end
            if get(handles.ui_errf_cmp_v_checkbox, 'Value')==0
                set(handles.ui_errf_pst_cmp_v_checkbox, 'Enable','Off')
            else
                set(handles.ui_errf_pst_cmp_v_checkbox, 'Enable','On')
            end
            if get(handles.ui_errf_cmp_w_checkbox, 'Value')==0
                set(handles.ui_errf_pst_cmp_w_checkbox, 'Enable','Off')
            else
                set(handles.ui_errf_pst_cmp_w_checkbox, 'Enable','On')
            end
        else
            set(handles.ui_errf_pst_text,           'Enable','Off')
            set(handles.ui_errf_pst_cmp_u_checkbox, 'Enable','Off')
            set(handles.ui_errf_pst_cmp_v_checkbox, 'Enable','Off')
            set(handles.ui_errf_pst_cmp_w_checkbox, 'Enable','Off')
            set(handles.ui_errf_pst_usehpf_checkbox,'Enable','Off')
            set(handles.ui_errf_pst_usehpf_checkbox, 'Value',0)
        end
        %
        if get(handles.ui_errf_pst_checkbox, 'Value')==1
            set(handles.ui_errf_pst_usehpf_checkbox,'Enable','On');
        else
            set(handles.ui_errf_pst_usehpf_checkbox,'Enable','Off');
        end
        if get(handles.ui_errf_pst_usehpf_checkbox, 'Value')==1
            set(handles.ui_errf_pst_hpfreq_edit,    'Enable','On');
            set(handles.ui_errf_pst_hpfreq_txt,     'Enable','On');
        else
            set(handles.ui_errf_pst_hpfreq_edit,    'Enable','Off');
            set(handles.ui_errf_pst_hpfreq_txt,     'Enable','Off');
        end
    %%%% Error filtering edit - file needed:
    elseif get(handles.ui_errf_edit_n_cancel_toggle_b, 'Value')==1
        %Get errf_types to enable:
        f_idx  = handles.sorted_index(get(handles.ui_listbox1,'Value'));
        ef_vals = strcmpi(...
            handles.viper_props.errf_types_enable_matrix(:,1),...
            handles.ts( f_idx ).p.dat_props.type );
        errf_types_to_enable = ...
            handles.viper_props.errf_types_supported( ...
            handles.viper_props.errf_types_enable_matrix{ef_vals ,2}, 1);
        %Enabling buttons
        set(handles.ui_errf_edit_n_cancel_toggle_b,'Enable','On');
        set(handles.ui_errf_reload_par_pushbutton, 'Enable','On');
        set(handles.ui_errf_apply_pushbutton,      'Enable','On');
        %universal
        set(handles.ui_errf_cmp_u_checkbox,     'Enable','On');
        set(handles.ui_errf_cmp_v_checkbox,     'Enable','On');
        set(handles.ui_errf_cmp_w_checkbox,     'Enable','On');
        %cormin || coravg:
        if sum( strcmpi( errf_types_to_enable,'cormin') )>0 ||...
                sum( strcmpi( errf_types_to_enable,'coravg') )>0
            set(handles.ui_errf_cor_text,           'Enable','On');
            set(handles.ui_errf_cor_text_sign,      'Enable','On');
            set(handles.ui_errf_cor_min_checkbox,   'Enable','On');
            set(handles.ui_errf_cor_avg_checkbox,   'Enable','On');
            set(handles.ui_errf_cor_lambda_edit,    'Enable','On');
        end
        %crop
        if sum( strcmpi( errf_types_to_enable,'crop') )>0
            set(handles.ui_errf_crop_checkbox,      'Enable','On');
            if get( handles.ui_errf_crop_checkbox, 'Value')==1
                set( handles.ui_errf_crop_1_edit,  'Enable','On');
                set( handles.ui_errf_crop_2_edit,  'Enable','On');
            else
                set( handles.ui_errf_crop_1_edit,  'Enable','Off');
                set( handles.ui_errf_crop_2_edit,  'Enable','Off');
            end
        end
        %pst
        if sum( strcmpi( errf_types_to_enable,'pst') )>0
            set(handles.ui_errf_pst_checkbox,       'Enable','On');
            if get(handles.ui_errf_cmp_u_checkbox, 'Value')==0
                %%% pst
                set(handles.ui_errf_pst_cmp_u_checkbox, 'Value',0)
                set(handles.ui_errf_pst_cmp_u_checkbox, 'Enable','Off')
            else
                set(handles.ui_errf_pst_cmp_u_checkbox, 'Enable','On')
            end
            if get(handles.ui_errf_cmp_v_checkbox, 'Value')==0
                %%% pst
                set(handles.ui_errf_pst_cmp_v_checkbox, 'Value',0)
                set(handles.ui_errf_pst_cmp_v_checkbox, 'Enable','Off')
            else
                set(handles.ui_errf_pst_cmp_v_checkbox, 'Enable','On')
            end
            if get(handles.ui_errf_cmp_w_checkbox, 'Value')==0
                %%% pst
                set(handles.ui_errf_pst_cmp_w_checkbox, 'Value',0)
                set(handles.ui_errf_pst_cmp_w_checkbox, 'Enable','Off')
            else
                set(handles.ui_errf_pst_cmp_w_checkbox, 'Enable','On')
            end
            %
            if get(handles.ui_errf_pst_checkbox, 'Value')==1
                set(handles.ui_errf_pst_text,           'Enable','On')
                if get(handles.ui_errf_cmp_u_checkbox, 'Value')==0
                    set(handles.ui_errf_pst_cmp_u_checkbox, 'Enable','Off')
                else
                    set(handles.ui_errf_pst_cmp_u_checkbox, 'Enable','On')
                end
                if get(handles.ui_errf_cmp_v_checkbox, 'Value')==0
                    set(handles.ui_errf_pst_cmp_v_checkbox, 'Enable','Off')
                else
                    set(handles.ui_errf_pst_cmp_v_checkbox, 'Enable','On')
                end
                if get(handles.ui_errf_cmp_w_checkbox, 'Value')==0
                    set(handles.ui_errf_pst_cmp_w_checkbox, 'Enable','Off')
                else
                    set(handles.ui_errf_pst_cmp_w_checkbox, 'Enable','On')
                end
            else
                set(handles.ui_errf_pst_text,           'Enable','Off')
                set(handles.ui_errf_pst_cmp_u_checkbox, 'Enable','Off')
                set(handles.ui_errf_pst_cmp_v_checkbox, 'Enable','Off')
                set(handles.ui_errf_pst_cmp_w_checkbox, 'Enable','Off')
                set(handles.ui_errf_pst_usehpf_checkbox,'Enable','Off')
                set(handles.ui_errf_pst_usehpf_checkbox, 'Value',0)
            end
            %
            if get(handles.ui_errf_pst_checkbox, 'Value')==1
                set(handles.ui_errf_pst_usehpf_checkbox,'Enable','On');
            else
                set(handles.ui_errf_pst_usehpf_checkbox,'Enable','Off');
            end
            if get(handles.ui_errf_pst_usehpf_checkbox, 'Value')==1
                set(handles.ui_errf_pst_hpfreq_edit,    'Enable','On');
                set(handles.ui_errf_pst_hpfreq_txt,     'Enable','On');
            else
                set(handles.ui_errf_pst_hpfreq_edit,    'Enable','Off');
                set(handles.ui_errf_pst_hpfreq_txt,     'Enable','Off');
            end
        end
    else
        set(handles.ui_errf_apply_pushbutton,      'Enable','Off');
        set(handles.ui_errf_reload_par_pushbutton, 'Enable','Off');
        %
        set(handles.ui_errf_cmp_u_checkbox,     'Enable','Off');
        set(handles.ui_errf_cmp_v_checkbox,     'Enable','Off');
        set(handles.ui_errf_cmp_w_checkbox,     'Enable','Off');
        set(handles.ui_errf_crop_checkbox,      'Enable','Off');
        set(handles.ui_errf_crop_1_edit,        'Enable','Off');
        set(handles.ui_errf_crop_2_edit,        'Enable','Off');
        set(handles.ui_errf_pst_checkbox,       'Enable','Off');
        set(handles.ui_errf_pst_text,           'Enable','Off');
        set(handles.ui_errf_pst_cmp_u_checkbox, 'Enable','Off');
        set(handles.ui_errf_pst_cmp_v_checkbox, 'Enable','Off');
        set(handles.ui_errf_pst_cmp_w_checkbox, 'Enable','Off');
        set(handles.ui_errf_pst_usehpf_checkbox,'Enable','Off');
        set(handles.ui_errf_pst_hpfreq_edit,    'Enable','Off');
        set(handles.ui_errf_pst_hpfreq_txt,     'Enable','Off');
        set(handles.ui_errf_cor_text,           'Enable','Off');
        set(handles.ui_errf_cor_text_sign,      'Enable','Off');
        set(handles.ui_errf_cor_min_checkbox,   'Enable','Off');
        set(handles.ui_errf_cor_avg_checkbox,   'Enable','Off');
        set(handles.ui_errf_cor_lambda_edit,    'Enable','Off');
    end
    %%%% Coordinate group + Err.filt defaults
    if get(handles.ui_cds_edit_n_cancel_togglebutton, 'Value')==1
        set(handles.ui_cds_edit_n_cancel_togglebutton,'Enable','On');
        set(handles.ui_cds_apply_pushbutton,       'Enable','On');
        %
        set(handles.ui_cds_fil_x_edit,'Enable','On');
        set(handles.ui_cds_fil_y_edit,'Enable','On');
        set(handles.ui_cds_fil_z_edit,'Enable','On');
        set(handles.ui_cds_ofs_x_edit,'Enable','On');
        set(handles.ui_cds_ofs_y_edit,'Enable','On');
        set(handles.ui_cds_ofs_z_edit,'Enable','On');
        set(handles.ui_cds_mpt_x_edit,'Enable','On');
        set(handles.ui_cds_mpt_y_edit,'Enable','On');
        set(handles.ui_cds_mpt_z_edit,'Enable','On');
        set(handles.ui_cds_dir_posx_radiobutton,'Enable','On');
        set(handles.ui_cds_dir_negx_radiobutton,'Enable','On');
        set(handles.ui_cds_dir_posy_radiobutton,'Enable','On');
        set(handles.ui_cds_dir_negy_radiobutton,'Enable','On');
    else
        set(handles.ui_cds_apply_pushbutton,       'Enable','Off');
        %
        set(handles.ui_cds_fil_x_edit,'Enable','Off');
        set(handles.ui_cds_fil_y_edit,'Enable','Off');
        set(handles.ui_cds_fil_z_edit,'Enable','Off');
        set(handles.ui_cds_ofs_x_edit,'Enable','Off');
        set(handles.ui_cds_ofs_y_edit,'Enable','Off');
        set(handles.ui_cds_ofs_z_edit,'Enable','Off');
        set(handles.ui_cds_mpt_x_edit,'Enable','Off');
        set(handles.ui_cds_mpt_y_edit,'Enable','Off');
        set(handles.ui_cds_mpt_z_edit,'Enable','Off');
        set(handles.ui_cds_dir_posx_radiobutton,'Enable','Off');
        set(handles.ui_cds_dir_negx_radiobutton,'Enable','Off');
        set(handles.ui_cds_dir_posy_radiobutton,'Enable','Off');
        set(handles.ui_cds_dir_negy_radiobutton,'Enable','Off');
    end

end
%DVL - after loading first file
if isfield(handles, 'sorted_index')==1
    f_idx = handles.sorted_index(get(handles.ui_listbox1,'Value'));
    if isfield( handles, 'dvl')==1
        [ ~, ~ ] = viper_dvl('dvl_update_ui_enable', handles, { f_idx } );
    end
end



% ------------------------------------------------------------
% Main-GUI/Duration view: ui updates
% ------------------------------------------------------------
function gui_update_estimator_display( f_idx, handles )
%Get usr error limit
err_abs_usr = str2double( get(handles.ui_est_edit,'String') );
%Get usr variable:
variables = cellstr(get(handles.ui_est_variable_popupmenu,'String'));
% variables{ get(handles.ui_est_variable_popupmenu,'Value') }
if strcmp( variables{ get(handles.ui_est_variable_popupmenu,'Value') },...
        'mean' )==1
    if isnan( err_abs_usr )==1
        dur_u = NaN;
        dur_v = NaN;
        dur_w = NaN;
    else
        %Calculate durations: 2*T_i*var(u_i)/(err_abs/2)^2
        dur_u = 2 * handles.ts( f_idx ).valu_data( 25 ) *...
            handles.ts( f_idx ).valu_data( 4 ) /...
            (err_abs_usr/2)^2;
        dur_v = 2 * handles.ts( f_idx ).valu_data( 26 ) *...
            handles.ts( f_idx ).valu_data( 5 ) /...
            (err_abs_usr/2)^2;
        dur_w = 2 * handles.ts( f_idx ).valu_data( 27 ) *...
            handles.ts( f_idx ).valu_data( 6 ) /...
            (err_abs_usr/2)^2;
    end
    %Set unit:
    set( handles.ui_est_unit_text, 'String', 'cm/s')
end
if strcmp( variables{ get(handles.ui_est_variable_popupmenu,'Value') },...
        'variance' )==1
    if isnan( err_abs_usr )==1
        dur_u = NaN;
        dur_v = NaN;
        dur_w = NaN;
    else
        %Calculate durations: 2*T_i*var(u_i)/(err_abs/2)^2
        dur_u = 2 * handles.ts( f_idx ).valu_data( 25 ) *...
            (handles.ts( f_idx ).valu_data( 28 ) - ...
             handles.ts( f_idx ).valu_data( 4 )^2) /...
            (err_abs_usr/2)^2;
        dur_v = 2 * handles.ts( f_idx ).valu_data( 26 ) *...
            (handles.ts( f_idx ).valu_data( 29 ) - ...
             handles.ts( f_idx ).valu_data( 5 )^2) /...
            (err_abs_usr/2)^2;
        dur_w = 2 * handles.ts( f_idx ).valu_data( 27 ) *...
            (handles.ts( f_idx ).valu_data( 30 ) - ...
             handles.ts( f_idx ).valu_data( 6 )^2) /...
            (err_abs_usr/2)^2;
    end
    %Set unit:
    set( handles.ui_est_unit_text, 'String', 'cm2/s2')
end
%Display results:
set(handles.ui_est_u_edit,'String',...
    sprintf('%.0f',  dur_u ));
set(handles.ui_est_v_edit,'String',...
    sprintf('%.0f',  dur_v ));
set(handles.ui_est_w_edit,'String',...
    sprintf('%.0f',  dur_w ));










% ------------------------------------------------------------
% Main-GUI/Errorfiltering: Collect displayed err.filt parameters:
% ------------------------------------------------------------
function [curr_p_errfilt] = gui_errfiltpars_collect( handles )
% Case modifying this: a copy of this fcn is also in viper_dvl!
% f_idx can be empty if just modifying defaults parameters
%Collects err.filt parameters
curr_p_errfilt.used_types = {};
% - Used components
curr_p_errfilt.used_velcomp(1,1) =...
    get( handles.ui_errf_cmp_u_checkbox,'Value' );
curr_p_errfilt.used_velcomp(2,1) =...
    get( handles.ui_errf_cmp_v_checkbox,'Value' );
curr_p_errfilt.used_velcomp(3,1) =...
    get( handles.ui_errf_cmp_w_checkbox,'Value' );
%At lest one comp is on - this is checked at checkbox_function
% - Crop
if get(handles.ui_errf_crop_checkbox, 'Value')==1
    ef_num = length( curr_p_errfilt.used_types ) + 1;
    curr_p_errfilt.used_types{ ef_num }='crop';
    %Name
    ef_mrk = strcmpi(...
        handles.viper_props.errf_types_supported(:,1),...
        curr_p_errfilt.used_types{ ef_num } );
    curr_p_errfilt.typ( ef_num ).name =...
        handles.viper_props.errf_types_supported( ef_mrk ,2);
    %Pars
    curr_p_errfilt.typ( ef_num ).lambda(1) = str2double(...
        get(handles.ui_errf_crop_1_edit,'String') );
    curr_p_errfilt.typ( ef_num ).lambda(2) = str2double(...
        get(handles.ui_errf_crop_2_edit,'String') );
end
% - PST based spike filter
if get(handles.ui_errf_pst_checkbox, 'Value')==1
    ef_num = length( curr_p_errfilt.used_types ) + 1;
    curr_p_errfilt.used_types{ ef_num }='PST';
    %Name
    ef_mrk = strcmpi(...
        handles.viper_props.errf_types_supported(:,1),...
        curr_p_errfilt.used_types{ ef_num } );
    curr_p_errfilt.typ( ef_num ).name =...
        handles.viper_props.errf_types_supported( ef_mrk ,2);
    %Pars
    curr_p_errfilt.typ( ef_num ).used_comp(1,1) =...
        get( handles.ui_errf_pst_cmp_u_checkbox,'Value');
    curr_p_errfilt.typ( ef_num ).used_comp(2,1) =...
        get( handles.ui_errf_pst_cmp_v_checkbox,'Value');
    curr_p_errfilt.typ( ef_num ).used_comp(3,1) =...
        get( handles.ui_errf_pst_cmp_w_checkbox,'Value');
    curr_p_errfilt.typ( ef_num ).lambda = NaN;%Not stored
    %highpass:
    if get(handles.ui_errf_pst_usehpf_checkbox,'Value')==1 &&...
            isnan( str2double(...
            get(handles.ui_errf_pst_hpfreq_edit,'String') ))==0
        curr_p_errfilt.typ( ef_num ).hipass_frq = ...
            str2double( get(handles.ui_errf_pst_hpfreq_edit,'String') );
    else
        curr_p_errfilt.typ( ef_num ).hipass_frq = [];
    end
end
% - Filter out samples where any cor is lower than
if get(handles.ui_errf_cor_min_checkbox,        'Value')==1
    ef_num = length( curr_p_errfilt.used_types ) + 1;
    curr_p_errfilt.used_types{ ef_num }='cormin';
    %Name
    ef_mrk = strcmpi(...
        handles.viper_props.errf_types_supported(:,1),...
        curr_p_errfilt.used_types{ ef_num } );
    curr_p_errfilt.typ( ef_num ).name =...
        handles.viper_props.errf_types_supported( ef_mrk ,2);
    %Pars
    curr_p_errfilt.typ( ef_num ).lambda =...
        str2double( get(handles.ui_errf_cor_lambda_edit,'String') );
end
% - Filter out samples where the average cor is lower than
if get(handles.ui_errf_cor_avg_checkbox,        'Value')==1
    ef_num = length( curr_p_errfilt.used_types ) + 1;
    curr_p_errfilt.used_types{ ef_num }='coravg';
    %Name
    ef_mrk = strcmpi(...
        handles.viper_props.errf_types_supported(:,1),...
        curr_p_errfilt.used_types{ ef_num } );
    curr_p_errfilt.typ( ef_num ).name =...
        handles.viper_props.errf_types_supported( ef_mrk ,2);
    %Pars
    curr_p_errfilt.typ( ef_num ).lambda =...
        str2double( get(handles.ui_errf_cor_lambda_edit,'String') );
end
if isfield( handles, 'dvl')==1
    [ ~, curr_p_errfilt ] = viper_dvl('dvl_errfiltpars_guicollect',...
        handles, { curr_p_errfilt });
end




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 
function chapter_03_plot_figures_related_functions
% 
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ------------------------------------------------------------
% Plots: Initialize plot properties: axes limits, fig & axes sizes
% ------------------------------------------------------------
function [ plots_siz, new_figpos, new_axspos ] = gui_plots_init_sizes(...
    plots_siz, plot_to_init, handles )
%If plot_to_init==0: all plots
% 1. Create a temp fig
% 2. using temp_fig -> derive maximal possible fig sizes
% 3. create each plot into temp_fig and derive needed axes margins
% ------------------------------------------------------------
% 1. Create temp_fig of maximal possible size:
% ------------------------------------------------------------
temp_fig = figure('Units','pixels',...
    'PaperPositionMode','auto','Color',[1,1,1],'MenuBar','none',...
    'Toolbar','figure','InvertHardcopy','off','NumberTitle','Off',...
    'OuterPosition', [ handles.gui.lb_min,...
    handles.gui.wh_max.*[ handles.plots.siz.used_screen_wr 1] ],...
    'Visible','off');
% set( fig_tmp, 'WindowState', 'normal');%==non-maximized
% ------------------------------------------------------------
% 2. Derive available screen size for figs, fig_border_sizes, font size
% ------------------------------------------------------------
temp_fig_pos = get( temp_fig, 'Position');
dp_plt = get( temp_fig, 'OuterPosition') - temp_fig_pos;
siz_fig_plt_borders =...
    [-dp_plt(1) -dp_plt(2) dp_plt(3)-(-dp_plt(1)) dp_plt(4)-(-dp_plt(2))];
%space needed at top of figure window:
plots_siz.ontop = siz_fig_plt_borders(4) - siz_fig_plt_borders(2);
%Maximal plot-figure position and size:  left bottom width heihgt
plots_siz.lb_min = [ 1, 1 + handles.gui.siz_taskbar ];
fig1out = get( handles.figure1, 'OuterPosition');
plots_siz.wh_max = [...
    (handles.gui.pos_screen(3) - fig1out(3))*...
    handles.plots.siz.used_screen_wr,...
    handles.gui.pos_screen(4) - handles.gui.siz_taskbar - plots_siz.ontop];
%Font size for 10 pt fonts and max_width:
plots_siz.ft_px = handles.plots.siz.ft_pt *...
    plots_siz.wh_max(1)/max(handles.plots.siz.printwidth_mm)*(1/72*25);
%Proportion font pt->px
tmp_text = uicontrol('Parent',temp_fig, 'Style','text', 'Units','pixels',...
    'FontUnits','pixels','FontSize',plots_siz.ft_px,'String','AAA');
set( tmp_text,'FontUnits','points')
plots_siz.ft_px2pt = get( tmp_text,'FontSize')/plots_siz.ft_px;
% plots_siz.ft_px_pt = get( tmp_text,'FontSize')
%A defaults ax_mrgn for DVL:
plots_siz.dflt_axmrgns = zeros(1,4);
% ------------------------------------------------------------
% 3. Create each axes into fig_tmp and derive required ax_margins:
%    (only for real axes - for tables: see later)
% ------------------------------------------------------------
% copy the main .plots as basis for temporary adjustments:
temp_plots = handles.plots;
temp_plots.siz = plots_siz;
% set default offsets and the supplemental margins:
default_offset = 4*plots_siz.ft_px;
supl_margin = 0.5*plots_siz.ft_px;
%init variables
temp_plots.axspos = cell( size( handles.plots.axs_h ) );
temp_axs_margins  = cell( size( handles.plots.axs_h ) );
for p_idx=2:size( handles.plots.axs_h, 1)%Plot-types
if handles.plots.numaxs( p_idx )>0
    %a. adjust temp_plots for the temporary creation of each axes
    for s_i=1:handles.plots.numaxs( p_idx )%Subplots:
        temp_plots.axspos{p_idx,s_i} = ...
            [default_offset default_offset ...
            ( temp_fig_pos(3) - 1.5*default_offset )*...
                handles.plots.siz.axs_wr( p_idx, s_i ) ...
            ( temp_fig_pos(3) - 1.5*default_offset )*...
                handles.plots.siz.axs_wr( p_idx, s_i )*...
                handles.plots.siz.axs_ar( p_idx, s_i )];
    end
    %adjustments:
    if p_idx==2
        temp_plots.axspos{2,2}(4) = temp_plots.axspos{2,1}(4);
    end
    %c. create axes to temp_fig:
    [ temp_plots ] = gui_plots_create_figs( handles.figure1,...
        handles, temp_plots, p_idx, temp_fig  );
    %depending on numbner of axes: 
    if handles.plots.numaxs( p_idx )==1%one axes:
        %d. get axes margins (incl. labels)
        temp_axs_margins{ p_idx, 1:handles.plots.numaxs( p_idx ) } = get(...
            temp_plots.axs_h( p_idx, 1:handles.plots.numaxs( p_idx )), ...
            'TightInset')';
        %collect max values of first axes as default for quick fig creation
        plots_siz.dflt_axmrgns =...
            max( plots_siz.dflt_axmrgns, temp_axs_margins{ p_idx,1}' );
    else
        %d. get axes margins (incl. labels)
        temp_axs_margins( p_idx, 1:handles.plots.numaxs( p_idx )) = get(...
            temp_plots.axs_h( p_idx,1:handles.plots.numaxs( p_idx )),...
            'TightInset')';
        %collect max values of first axes as default for quick fig creation
        plots_siz.dflt_axmrgns =...
            max( plots_siz.dflt_axmrgns, temp_axs_margins{ p_idx,1} );
    end
    %adjustments:
    if p_idx==2
        temp_axs_margins{2,3}(1) = temp_axs_margins{2,1}(1);
    end
end
end
delete( temp_fig )
% ------------------------------------------------------------
% 4. Calculate required axes positions, then figure sizes
% ------------------------------------------------------------
new_axspos = cell( size( handles.plots.axs_h ) );
new_figpos = cell( size( handles.plots.fig_h ) );
%Go through axs - first exceptions:
%Plot==1: containing a table
if plot_to_init==0
    p_idx=1;
    h_fig_t = plots_siz.wh_max(2)/3;%fig height
    w_fig_t = plots_siz.wh_max(1);%fig width
    mrg_t = 5;%margins
    spc_t = 20;%space between tables
    h4tables = h_fig_t - 2*mrg_t;
    w4tables = w_fig_t - 2*mrg_t - spc_t;
    w_table1 = 0.33*w4tables;
    w_table2 = 0.66*w4tables;
    %Fill figpos variable:
    new_figpos{ p_idx } = ...
        [ plots_siz.lb_min w_fig_t h_fig_t ];
    %Fill axpos variable:
    new_axspos( p_idx, 1:2 ) = {...
        [ mrg_t  mrg_t  w_table1  h4tables],...
        [ mrg_t+w_table1+spc_t  mrg_t  w_table2  h4tables] };
end
%Plot==5: spectra - changes plot figure size at limit change
if plot_to_init==0 || plot_to_init==5
    p_idx=5; %spectra
    %If just modifying -> copy figpos & axpos
    if plot_to_init==5
        new_axspos = handles.plots.axspos;
        new_figpos = handles.plots.figpos;
    end
	%!!! Spectra dataaspectratio = 1, therefore first determine,
    %!!! whether required width or height are limiting
    %a. convert to numeric: (columns: le-bo-ri-to)
    curr_axs_margins = cell2mat(...
        temp_axs_margins( p_idx, 1:handles.plots.numaxs( p_idx )))';
    %WIDTH:
    %sum up margins within plot:
    sum_axs_w_margins = curr_axs_margins(1) + curr_axs_margins(3) +...
        supl_margin*(handles.plots.numaxs( p_idx )+1);
    %calc available width for axs_areas (using max fig width)
    w_for_axs_area =...
        plots_siz.wh_max(1)*...
        handles.plots.siz.printwidth_mm( p_idx )/...
        max( handles.plots.siz.printwidth_mm ) - ...
        sum_axs_w_margins;
    %HEIGHT:
    %sum up margins within plot:
    sum_axs_h_margins = curr_axs_margins(2) + curr_axs_margins(4) +...
        supl_margin*(handles.plots.numaxs( p_idx )+1);
    %calc available height for axs_areas (using max fig width)
    h_for_axs_area =  plots_siz.wh_max(2) - sum_axs_h_margins;
    %Compare h/w aspect ratio with limits h/w aspect ratio
    lims_ar = diff( handles.plots.axlims_usrdef{  p_idx ,1 }(2,:) )/...
        diff( handles.plots.axlims_usrdef{  p_idx ,1 }(1,:) );
    if w_for_axs_area * lims_ar > h_for_axs_area%h is limiting:
        %axes height:
        new_axs_h = h_for_axs_area;
        %axes width:
        new_axs_w = h_for_axs_area/lims_ar;
        %fig_h
        new_fig_h = plots_siz.wh_max(2);
        %fig_w:
        new_fig_w = new_axs_w + sum_axs_w_margins;
    else
        %axes height:
        new_axs_w = w_for_axs_area;
        %axes width:
        new_axs_h = w_for_axs_area*lims_ar;
        %fig_h
        new_fig_w = new_axs_w + sum_axs_w_margins;
        %fig_w:
        new_fig_h = new_axs_h + sum_axs_h_margins;
    end
    %ax_le:
    new_axs_le = curr_axs_margins(1) + supl_margin;
    %ax_bo:
    new_axs_bo = curr_axs_margins(2) + supl_margin;
    %Concatenate axs pos
    new_axspos_mat =...
        vertcat(new_axs_le,new_axs_bo,new_axs_w,new_axs_h);
    %Fill axpos variable:
    new_axspos( p_idx, 1:handles.plots.numaxs( p_idx ) ) = mat2cell(...
        reshape( new_axspos_mat, [1, handles.plots.numaxs( p_idx )*4 ]),...
        [1], repmat( 4, [1,handles.plots.numaxs( p_idx )]) );
    %Fill figpos variable:
    new_figpos{ p_idx } = ...
        [ plots_siz.lb_min new_fig_w new_fig_h];
end
%All other plots:
if plot_to_init==0
for p_idx=1:size( handles.plots.axs_h, 1)%Plot-types
%Skip exceptions (above)
if p_idx==1 || p_idx==5
    continue
end
if handles.plots.numaxs( p_idx )>0
    %subplot arrangement and plot type dependent:
    if handles.plots.numaxs( p_idx )==1%one axes:
        %a. convert to numeric: (columns: le-bo-ri-to)
        curr_axs_margins = cell2mat(...
            temp_axs_margins( p_idx, 1:handles.plots.numaxs( p_idx )))';
        %WIDTH:
        %b. use full available width
        new_fig_w =...
            plots_siz.wh_max(1)*...
            handles.plots.siz.printwidth_mm( p_idx )/...
            max( handles.plots.siz.printwidth_mm );
        %c. sum up margins within plot:
        sum_axs_w_margins = curr_axs_margins(1) + curr_axs_margins(3) +...
            supl_margin*(handles.plots.numaxs( p_idx )+1);
        %d. calc available width for axs_areas (using max fig width)
        w_for_axs_area = new_fig_w - sum_axs_w_margins;
        %e. calc axes widths:
        new_axs_w = w_for_axs_area *...
            handles.plots.siz.axs_wr(p_idx,1:handles.plots.numaxs(p_idx));      
        %HEIGHT:
        %f. calc axes heights:
        new_axs_h = w_for_axs_area *...
            handles.plots.siz.axs_wr(p_idx,1:handles.plots.numaxs(p_idx)).*...
            handles.plots.siz.axs_ar(p_idx,1:handles.plots.numaxs(p_idx));
        %LE, BO:
        %g. ax_le:
        new_axs_le = curr_axs_margins(:,1) + supl_margin;
        %h. determine common axes bo:
        new_axs_bo = curr_axs_margins(:,2) + supl_margin;
        %i. Concatenate axs pos
        new_axspos_mat =...
            vertcat(new_axs_le,new_axs_bo,new_axs_w,new_axs_h);
        %j. determine required figure height
        new_fig_h = curr_axs_margins(2) + curr_axs_margins(4) + ...
            new_axs_h + supl_margin*(handles.plots.numaxs( p_idx )+1);
    elseif p_idx==3 %subplots over each other:
        %a. convert to numeric: (columns: le-bo-ri-to)
        curr_axs_margins = cell2mat(...
            temp_axs_margins( p_idx, 1:handles.plots.numaxs( p_idx ))');
        %WIDTH:
        %b. use full available width
        new_fig_w =...
            plots_siz.wh_max(1)*...
            handles.plots.siz.printwidth_mm( p_idx )/...
            max( handles.plots.siz.printwidth_mm );
        %c. sum up margins horizontally within fig & add suppl margins:
        sum_max_axes_w_margins =...
            max( curr_axs_margins(:,1) ) +...
            max( curr_axs_margins(:,3) ) + ...
            supl_margin + supl_margin;
        %d. calc available width for axs_areas (using max fig width)
        w_for_axs_area = new_fig_w - sum_max_axes_w_margins;
        %e. calc axes widths:
        new_axs_w = w_for_axs_area *...
            handles.plots.siz.axs_wr(p_idx,1:handles.plots.numaxs(p_idx));
        %HEIGHT:
        %f. calc axes heights:
        new_axs_h = w_for_axs_area *...
            handles.plots.siz.axs_wr(p_idx,1:handles.plots.numaxs(p_idx)).*...
            handles.plots.siz.axs_ar(p_idx,1:handles.plots.numaxs(p_idx));
        %LE, BO:
        %g. ax_le:
        new_axs_le = repmat( max( curr_axs_margins(:,1) ),...
            [1, handles.plots.numaxs( p_idx ) ]) + ...
            repmat( supl_margin, [1 handles.plots.numaxs( p_idx )]);
        %h. determine common axes bo:
        new_axs_bo = cumsum( curr_axs_margins(:,2)' + ...
            repmat( supl_margin, [1 handles.plots.numaxs( p_idx )]) +...
            [0 new_axs_h(1:handles.plots.numaxs( p_idx )-1)] );
        %i. Concatenate axs pos with reverse order, so that uppest == 1:
        new_axspos_mat =...
            vertcat(new_axs_le,new_axs_bo,new_axs_w,new_axs_h);
        new_axspos_mat =...
            new_axspos_mat(:,(handles.plots.numaxs( p_idx ):-1:1));
        %j. determine required figure height
        new_fig_h = sum(...
            curr_axs_margins(:,2)' + ...
            repmat( supl_margin, [1 handles.plots.numaxs( p_idx )]) + ...
            new_axs_h + ...
            curr_axs_margins(:,4)') + supl_margin;
    else%subplots beside each other:
        %a. convert to numeric: (columns: le-bo-ri-to)
        curr_axs_margins = cell2mat(...
            temp_axs_margins( p_idx, 1:handles.plots.numaxs( p_idx ))');
        %WIDTH:
        %b. use full available width
        new_fig_w =...
            plots_siz.wh_max(1)*...
            handles.plots.siz.printwidth_mm( p_idx )/...
            max( handles.plots.siz.printwidth_mm );
        %c. sum up margins horizontally within fig & add suppl margins:
        sum_axs_margins = sum( curr_axs_margins, 1);
        sum_axs_w_margins = sum_axs_margins(1) + sum_axs_margins(3)+...
            supl_margin*(handles.plots.numaxs( p_idx )+1);
        %d. calc available width for axs_areas (using max fig width)
        w_for_axs_area = new_fig_w - sum_axs_w_margins;
        %e. calc axes widths:
        new_axs_w = w_for_axs_area *...
            handles.plots.siz.axs_wr(p_idx,1:handles.plots.numaxs(p_idx));
        %HEIGHT:
        %f. calc axes heights:
        new_axs_h = w_for_axs_area *...
            handles.plots.siz.axs_wr(p_idx,1:handles.plots.numaxs(p_idx)).*...
            handles.plots.siz.axs_ar(p_idx,1:handles.plots.numaxs(p_idx));
        %LE, BO:
        %g. ax_le:
        new_axs_le = cumsum(curr_axs_margins(:,1)' +...
            repmat( supl_margin, [1 handles.plots.numaxs( p_idx )]) +...
            [ 0 new_axs_w(1:handles.plots.numaxs( p_idx )-1) ]);
        %h. determine common axes bo:
        new_axs_bo = repmat( max( curr_axs_margins(:,2) ), ...
            [1 handles.plots.numaxs( p_idx )]) + supl_margin;
        %i. Concatenate axs pos
        new_axspos_mat =...
            vertcat(new_axs_le,new_axs_bo,new_axs_w,new_axs_h);
        %j. determine required figure height
        new_fig_h = max(...
            new_axs_bo + new_axs_h + curr_axs_margins(:,4)', [], 2) +...
            supl_margin;
    end
    %Fill axpos variable:
    new_axspos( p_idx, 1:handles.plots.numaxs( p_idx ) ) = mat2cell(...
        reshape( new_axspos_mat, [1, handles.plots.numaxs( p_idx )*4 ]),...
        [1], repmat( 4, [1,handles.plots.numaxs( p_idx )]) );
    %Fill figpos variable:
    new_figpos{ p_idx } = ...
        [ plots_siz.lb_min new_fig_w new_fig_h];
    %adjustments:
    if p_idx==2
        new_axspos{2,2}(4) = new_axspos{2,1}(4);
    end
end
end
end





% ------------------------------------------------------------
% Plots: Create plot_figures:
% ------------------------------------------------------------
function [ curr_plots ] = gui_plots_create_figs( main_fig_h,...
    handles, curr_plots, fig2create, handle_test_fig)
% fig_idx:  create figure with given fig_idx
% handle_test_fig:
%       - for plot-figure creation: empty
%       - for testing -> dest_fig == handle_test_fig
% -------------------------------------------------------------------------
%%% Statistics box (pos 1 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if fig2create==1 || (curr_plots.fixed_figs==0 && ...
        fig2create==0 && get(handles.ui_plot_fig1_checkbox, 'Value')==1)
    plot_idx = 1;
    curr_plots.fig_h( plot_idx ) = figure(...
        'Position', curr_plots.figpos{ plot_idx },...
        'MenuBar','none',...
        'NumberTitle','Off', 'Name','Newfig' );
    %figure('Units','pixels',...
    %'PaperPositionMode','auto','Color',[1,1,1],...
    %'Toolbar','figure','InvertHardcopy','off',...
    set( curr_plots.fig_h( plot_idx ), 'CloseRequestFcn',...
        {@guiexec_plots_close,main_fig_h,...
        curr_plots.fig_h( plot_idx ), plot_idx });
    % Table 2:
    curr_plots.axs_h( plot_idx ,1) = uitable(...
        'Parent', curr_plots.fig_h( plot_idx ),...
        'Position', curr_plots.axspos{ plot_idx,1} );
    pos_no_comp = find( handles.field_props.valu_uvw == 0 );
    set(curr_plots.axs_h( plot_idx ,1),...
        'ColumnName',[],...
        'RowName',handles.field_props.valu_names( pos_no_comp ));
    % Table 1:
    curr_plots.axs_h( plot_idx ,2) = uitable(...
        'Parent', curr_plots.fig_h( plot_idx ),...
        'Position', curr_plots.axspos{ plot_idx ,2});
    set(curr_plots.axs_h( plot_idx ,2),...
        'ColumnName',{'u','v','w'},...
        'RowName',handles.field_props.valu_groups);
end
% -------------------------------------------------------------------------
% %%% QUALITY PLOT (pos 2 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if fig2create==2 || (curr_plots.fixed_figs==0 && ...
        fig2create==0 && get(handles.ui_plot_fig2_checkbox, 'Value')==1)
    plot_idx = 2;
    %Define destination figure:
    if isempty( handle_test_fig )==1
        curr_plots.fig_h( plot_idx ) = figure('Units','pixels',...
            'PaperPositionMode','auto','Color',[1,1,1],'MenuBar','none',...
            'Toolbar','figure','InvertHardcopy','off','NumberTitle','Off',...
            'Position',curr_plots.figpos{ plot_idx }, 'Name','Newfig' );
        set( curr_plots.fig_h( plot_idx ), 'CloseRequestFcn',...
            {@guiexec_plots_close,main_fig_h,...
            curr_plots.fig_h( plot_idx ), plot_idx});
    else
        curr_plots.fig_h( plot_idx ) = handle_test_fig;
    end
    %Create and set axes:
    if isempty(find(isnan( curr_plots.axlims_minmax{ plot_idx ,1} ),1))==1
        [ currlims ] = fcn_sokoray_fit_lims2data(...
            curr_plots.axlims_minmax{ plot_idx ,1}(:,1),...
            curr_plots.axlims_minmax{ plot_idx ,1}(:,2));
    else
        currlims = curr_plots.axlims_usrdef{ plot_idx ,1};
    end
    curr_plots.axs_h( plot_idx ,1) = axes(...
        'Parent',curr_plots.fig_h( plot_idx ),...
        'Units','pixels', 'Box', 'on',...
        'Position',curr_plots.axspos{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px,...
        'YLim',currlims(2,:),'YTick', [0 50 60 70 80 90 100],...
        'XLim',currlims(1,:),'YGrid','on','XGrid','on');
    xlabel( curr_plots.axs_h( plot_idx ,1),...
        curr_plots.axs_label_x{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    ylabel( curr_plots.axs_h( plot_idx ,1),...
        curr_plots.axs_label_y{ plot_idx,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    if isempty(find(isnan( curr_plots.axlims_usrdef{ plot_idx ,2} ),1))==1
        currlims = curr_plots.axlims_minmax{ plot_idx ,2};
    else
        currlims = curr_plots.axlims_usrdef{ plot_idx ,2};
    end
    
    curr_plots.axs_h( plot_idx ,2) = axes(...
        'Parent',curr_plots.fig_h( plot_idx ),...
        'Units','pixels', 'Box', 'on',...
        'Position',curr_plots.axspos{ plot_idx ,2},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels',...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px,...
        'YLim',currlims(2,:),'YTick', [0 50 60 70 80 90 100],...
        'YTickLabel',{''},...
        'XLim',currlims(1,:),'XTick', [0 100],'XTickLabel',{''},...
        'YGrid','off','XGrid','off');
    if isempty(find(isnan( curr_plots.axlims_minmax{ plot_idx ,3} ),1))==1
        currlims = curr_plots.axlims_minmax{ plot_idx ,3};
    else
        currlims = curr_plots.axlims_usrdef{ plot_idx ,3};
    end
    curr_plots.axs_h( plot_idx ,3) = axes(...
        'Parent',curr_plots.fig_h( plot_idx ),...
        'Units','pixels', 'Box', 'on',...
        'Position',curr_plots.axspos{ plot_idx ,3},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px,...
        'Color','none','YLim',currlims(2,:) ,'YTick', [0 5 10 15 20 25],...
        'XLim',currlims(1,:),'YGrid','on','XGrid','on');
    xlabel( curr_plots.axs_h( plot_idx ,3),...
        curr_plots.axs_label_x{ plot_idx ,3},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    ylabel( curr_plots.axs_h( plot_idx ,3),...
        curr_plots.axs_label_y{ plot_idx ,3},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
end
% -------------------------------------------------------------------------
% %%% TS Raw vs. (err.filtd & intpoltd) PLOT (pos 3 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if fig2create==3 || (curr_plots.fixed_figs==0 && ...
        fig2create==0 && get(handles.ui_plot_fig3_checkbox,'Value')==1)
    plot_idx = 3;
%     curr_plots.figpos{ plot_idx }
    %Define destination figure:
    if isempty( handle_test_fig )==1
        curr_plots.fig_h( plot_idx ) = figure('Units','pixels',...
            'PaperPositionMode','auto','Color',[1,1,1],'MenuBar','none',...
            'Toolbar','figure','InvertHardcopy','off','NumberTitle','Off',...
            'Position',curr_plots.figpos{ plot_idx }, 'Name','Newfig' );
        set( curr_plots.fig_h( plot_idx ), 'CloseRequestFcn',...
            {@guiexec_plots_close,main_fig_h,...
            curr_plots.fig_h( plot_idx ), plot_idx });
    else
        curr_plots.fig_h( plot_idx ) = handle_test_fig;
    end
    %Create and set axes:
    if isempty(find(isnan( curr_plots.axlims_minmax{ plot_idx ,1} ),1))==1
        %currlims = curr_plots.axlims_minmax{3,1}
        [ currlims ] = fcn_sokoray_fit_lims2data(...
            curr_plots.axlims_minmax{ plot_idx ,1}(:,1),...
            curr_plots.axlims_minmax{ plot_idx ,1}(:,2));
    else
        currlims = curr_plots.axlims_usrdef{ plot_idx ,1};
    end
    curr_plots.axs_h( plot_idx ,1) = axes(...
        'Parent',curr_plots.fig_h( plot_idx ),...
        'Units','pixels', 'Box','on', 'GridLineStyle',':',...
        'XGrid','on','YGrid','on','XMinorTick','on','YMinorTick','on',...
        'Position',curr_plots.axspos{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px,...
        'XLim',currlims(1,:),'YLim',currlims(2,:));
    xlabel( curr_plots.axs_h( plot_idx ,1),...
        curr_plots.axs_label_x{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    ylabel( curr_plots.axs_h( plot_idx ,1),...
        curr_plots.axs_label_y{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    curr_plots.axs_h( plot_idx ,2) = axes(...
        'Parent',curr_plots.fig_h( plot_idx ),...
        'Units','pixels', 'Box','on', 'GridLineStyle',':',...
        'XGrid','on','YGrid','on','XMinorTick','on','YMinorTick','on',...
        'Position',curr_plots.axspos{ plot_idx ,2},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px,...
        'XLim',currlims(1,:),'YLim',currlims(2,:));
    xlabel( curr_plots.axs_h( plot_idx ,2),...
        curr_plots.axs_label_x{ plot_idx ,2},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    ylabel( curr_plots.axs_h( plot_idx ,2),...
        curr_plots.axs_label_y{ plot_idx ,2},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    % ---------------------------------------------------------------------
    % Link axes
    % ---------------------------------------------------------------------
    linkaxes([curr_plots.axs_h( plot_idx ,1);...
        curr_plots.axs_h( plot_idx ,2)],  'xy');
end
% -------------------------------------------------------------------------
% %%% Velocity histograms (pos4 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if fig2create==4 || (curr_plots.fixed_figs==0 && ...
        fig2create==0 && get(handles.ui_plot_fig4_checkbox,'Value')==1)
    plot_idx = 4;
    %Define destination figure:
    if isempty( handle_test_fig )==1
        curr_plots.fig_h( plot_idx ) = figure('Units','pixels',...
            'PaperPositionMode','auto','Color',[1,1,1],'MenuBar','none',...
            'Toolbar','figure','InvertHardcopy','off','NumberTitle','Off',...
            'Position',curr_plots.figpos{ plot_idx }, 'Name','Newfig' );
        set( curr_plots.fig_h( plot_idx ), 'CloseRequestFcn',...
            {@guiexec_plots_close,main_fig_h,...
            curr_plots.fig_h( plot_idx ), plot_idx });
    else
        curr_plots.fig_h( plot_idx ) = handle_test_fig;
    end
    %Create and set axes:
    if isempty(find(isnan( curr_plots.axlims_minmax{ plot_idx ,1} ),1))==1
        %currlims = curr_plots.axlims_minmax{4,1};
        [ currlims ] = fcn_sokoray_fit_lims2data(...
            curr_plots.axlims_minmax{ plot_idx ,1}(:,1),...
            curr_plots.axlims_minmax{ plot_idx ,1}(:,2));
    else
        currlims = curr_plots.axlims_usrdef{ plot_idx ,1};
    end
    curr_plots.axs_h( plot_idx ,1) = axes(...
        'Parent',curr_plots.fig_h( plot_idx ),...
        'Units','pixels', 'Box', 'on',...
        'Position',curr_plots.axspos{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px,...
        'XLim',currlims(1,:),'YLim',currlims(2,:));
    xlabel( curr_plots.axs_h( plot_idx ,1),...
        curr_plots.axs_label_x{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    ylabel( curr_plots.axs_h( plot_idx ,1),...
        curr_plots.axs_label_y{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    if isempty(find(isnan( curr_plots.axlims_minmax{ plot_idx ,2} ),1))==1
        %currlims = curr_plots.axlims_minmax{4,2};
        [ currlims ] = fcn_sokoray_fit_lims2data(...
            curr_plots.axlims_minmax{ plot_idx ,2}(:,1),...
            curr_plots.axlims_minmax{ plot_idx ,2}(:,2));
    else
        currlims = curr_plots.axlims_usrdef{ plot_idx ,2};
    end
    curr_plots.axs_h( plot_idx ,2) = axes(...
        'Parent',curr_plots.fig_h( plot_idx ),...
        'Units','pixels', 'Box', 'on',...
        'Position',curr_plots.axspos{ plot_idx ,2},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px,...
        'XLim',currlims(1,:),'YLim',currlims(2,:));
    xlabel( curr_plots.axs_h( plot_idx ,2),...
        curr_plots.axs_label_x{ plot_idx ,2},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    ylabel( curr_plots.axs_h( plot_idx ,2),...
        curr_plots.axs_label_y{ plot_idx ,2},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    if isempty(find(isnan( curr_plots.axlims_minmax{ plot_idx ,3} ),1))==1
        %currlims = curr_plots.axlims_minmax{4,3};
        [ currlims ] = fcn_sokoray_fit_lims2data(...
            curr_plots.axlims_minmax{ plot_idx ,3}(:,1),...
            curr_plots.axlims_minmax{ plot_idx ,3}(:,2));
    else
        currlims = curr_plots.axlims_usrdef{ plot_idx ,3};
    end
    curr_plots.axs_h( plot_idx ,3) = axes(...
        'Parent',curr_plots.fig_h( plot_idx ),...
        'Units','pixels', 'Box', 'on',...
        'Position',curr_plots.axspos{ plot_idx ,3},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px,...
        'XLim',currlims(1,:),'YLim',currlims(2,:));
    xlabel( curr_plots.axs_h( plot_idx ,3),...
        curr_plots.axs_label_x{ plot_idx ,3},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    ylabel( curr_plots.axs_h( plot_idx ,3),...
        curr_plots.axs_label_y{ plot_idx ,3},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
end
% -------------------------------------------------------------------------
% %%% Spectrum PLOT (pos 5 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if fig2create==5 || (curr_plots.fixed_figs==0 && ...
        fig2create==0 && get(handles.ui_plot_fig5_checkbox,'Value')==1)
    plot_idx = 5;
    %Define destination figure:
    if isempty( handle_test_fig )==1
        % CREATE NEW:
        curr_plots.fig_h( plot_idx ) = figure('Units','pixels',...
            'PaperPositionMode','auto','Color',[1,1,1],'MenuBar','none',...
            'Toolbar','figure','InvertHardcopy','off','NumberTitle','Off',...
            'Position',curr_plots.figpos{ plot_idx }, 'Name','Newfig' );
        set( curr_plots.fig_h( plot_idx ), 'CloseRequestFcn',...
            {@guiexec_plots_close,main_fig_h,...
            curr_plots.fig_h( plot_idx ), plot_idx });
    else
        curr_plots.fig_h( plot_idx ) = handle_test_fig;
    end
    %Create and set axes:
    %Always using user-defined limits, otherwise figure recreated at update 
    powlims(1,:) = curr_plots.axlims_usrdef{  plot_idx ,1 }(1,:);
    powlims(2,:) = curr_plots.axlims_usrdef{  plot_idx ,1 }(2,:);
    currlims = 10.^powlims;
    curr_plots.axs_h( plot_idx ,1) = axes(...
        'Parent',curr_plots.fig_h( plot_idx ),...
        'Units','pixels', 'Box','on', 'GridLineStyle',':',...
        'XGrid','on','YGrid','on','XMinorTick','on','YMinorTick','on',...
        'Position',curr_plots.axspos{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px,...
        'XScale','log','YScale','log',...
        'XLim',currlims(1,:),'YLim',currlims(2,:));
    xlabel( curr_plots.axs_h( plot_idx ,1),...
        curr_plots.axs_label_x{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    ylabel( curr_plots.axs_h( plot_idx ,1),...
        curr_plots.axs_label_y{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    %Correcting axes ticklabels (assumption: ticklabel width == height)
    x_range = floor(log10(currlims(1,2))) - ceil(log10(currlims(1,1)));
    y_range = floor(log10(currlims(2,2))) - ceil(log10(currlims(2,1)));
    %x_siz_px = curr_plots.axspos{ plot_idx ,1}(3);
    %y_siz_px = curr_plots.axspos{ plot_idx ,1}(4);
    x_space4numtick =...
        curr_plots.axspos{ plot_idx ,1}(3)/(2*curr_plots.siz.ft_px);
    y_space4numtick =...
        curr_plots.axspos{ plot_idx ,1}(4)/(2*curr_plots.siz.ft_px);
    new_xticks = ...
        ( ceil(log10(currlims(1,1))) : ...
        ceil( x_range/x_space4numtick ) : floor(log10(currlims(1,2))) );
    new_yticks = ...
        ( ceil(log10(currlims(2,1))) : ...
        ceil( y_range/y_space4numtick ) : floor(log10(currlims(2,2))) );
    set( curr_plots.axs_h( plot_idx ,1), 'XTick', 10.^new_xticks );
    set( curr_plots.axs_h( plot_idx ,1), 'YTick', 10.^new_yticks );

end
% -------------------------------------------------------------------------
% %%% Cumulative mean PLOT (pos 6 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if fig2create==6 || (curr_plots.fixed_figs==0 && ...
        fig2create==0 && get(handles.ui_plot_fig6_checkbox,'Value')==1)
    plot_idx = 6;
    %Define destination figure:
    if isempty( handle_test_fig )==1
        curr_plots.fig_h( plot_idx ) = figure('Units','pixels',...
            'PaperPositionMode','auto','Color',[1,1,1],'MenuBar','none',...
            'Toolbar','figure','InvertHardcopy','off','NumberTitle','Off',...
            'Position',curr_plots.figpos{ plot_idx }, 'Name','Newfig' );
        set( curr_plots.fig_h( plot_idx ), 'CloseRequestFcn',...
            {@guiexec_plots_close,main_fig_h,...
            curr_plots.fig_h( plot_idx ), plot_idx });
    else
        curr_plots.fig_h( plot_idx ) = handle_test_fig;
    end
    %Create and set axes:
    if isempty(find(isnan( curr_plots.axlims_minmax{ plot_idx ,1} ),1))==1
        %currlims = curr_plots.axlims_minmax{6,1};
        [ currlims ] = fcn_sokoray_fit_lims2data(...
            curr_plots.axlims_minmax{ plot_idx ,1}(:,1),...
            curr_plots.axlims_minmax{ plot_idx ,1}(:,2));
    else
        currlims = curr_plots.axlims_usrdef{ plot_idx ,1};
    end
    curr_plots.axs_h( plot_idx ,1) = axes(...
        'Parent',curr_plots.fig_h( plot_idx ),...
        'Units','pixels', 'Box','on', 'GridLineStyle',':',...
        'XGrid','on','YGrid','on','XMinorTick','on','YMinorTick','on',...
        'Position',curr_plots.axspos{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px,...
        'XLim',currlims(1,:),'YLim',currlims(2,:));
    xlabel( curr_plots.axs_h( plot_idx ,1),...
        curr_plots.axs_label_x{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    ylabel( curr_plots.axs_h( plot_idx ,1),...
        curr_plots.axs_label_y{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
end
% -------------------------------------------------------------------------
% %%% Autocorr function (pos 9 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if fig2create==9 || (curr_plots.fixed_figs==0 && ...
        fig2create==0 && get(handles.ui_plot_fig9_checkbox,'Value')==1)
    plot_idx = 9;
    %Define destination figure:
    if isempty( handle_test_fig )==1
        curr_plots.fig_h( plot_idx ) = figure('Units','pixels',...
            'PaperPositionMode','auto','Color',[1,1,1],'MenuBar','none',...
            'Toolbar','figure','InvertHardcopy','off','NumberTitle','Off',...
            'Position',curr_plots.figpos{ plot_idx }, 'Name','Newfig' );
        set( curr_plots.fig_h( plot_idx ), 'CloseRequestFcn',...
            {@guiexec_plots_close,main_fig_h,...
            curr_plots.fig_h( plot_idx ), plot_idx });
    else
        curr_plots.fig_h( plot_idx ) = handle_test_fig;
    end
    %Create and set axes:
    if isempty(find(isnan( curr_plots.axlims_minmax{ plot_idx ,1} ),1))==1
        %currlims = curr_plots.axlims_minmax{9,1};
        [ currlims ] = fcn_sokoray_fit_lims2data(...
            curr_plots.axlims_minmax{ plot_idx ,1}(:,1),...
            curr_plots.axlims_minmax{ plot_idx ,1}(:,2));
    else
        currlims = curr_plots.axlims_usrdef{ plot_idx ,1};
    end
    curr_plots.axs_h( plot_idx ,1) = axes(...
        'Parent',curr_plots.fig_h( plot_idx ),...
        'Units','pixels', 'Box','on', 'GridLineStyle',':',...
        'XGrid','on','YGrid','on','XMinorTick','on','YMinorTick','on',...
        'Position',curr_plots.axspos{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px,...
        'XLim',currlims(1,:),'YLim',currlims(2,:));
    xlabel( curr_plots.axs_h( plot_idx ,1),...
        curr_plots.axs_label_x{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    ylabel( curr_plots.axs_h( plot_idx ,1),...
        curr_plots.axs_label_y{ plot_idx ,1},...
        'FontName',curr_plots.ft_name,...
        'Fontunits','pixels','FontSize',curr_plots.siz.ft_px);
    set(curr_plots.axs_h( plot_idx ,1),...
        'Xlim',curr_plots.axlims_usrdef{  plot_idx ,1 }(1,:));
end


% ------------------------------------------------------------------------
% Plots: close plot figures
% WARNING: also called at figure1_CloseRequestFcn !!!
% ------------------------------------------------------------------------
function guiexec_plots_close( src, eventdata,...
    h_main_fig, handle_fig, plot_idx)
% ------------------------------------------------------------
% Executes, if a figureplot is closed by [x] or by unchecking checkbox
% ------------------------------------------------------------
% If uncheck checkbox
%     - input: ( [] , [] , handle_main_GUI, [], plot_idx )
% If close windows with [x] == closerequestFcn
%     - input: ( src, eventdata, handle_main_GUI, handle_fig, fig_idx )
%%%% Info:
% get first the current handles of the main gui
% ------------------------------------------------------------------------
% WARNING:update_fig_closing also called at figure1_CloseRequestFcn !!!
% ------------------------------------------------------------------------

% ------------------------------------------------------------
% %%%% Get handles of main GUI:
% ------------------------------------------------------------
handles = guidata( h_main_fig );
% ------------------------------------------------------------
% %%%% Save current figure le-bo-positions as default:
% ------------------------------------------------------------
figpos = get( handles.plots.fig_h( plot_idx ), 'Position' );
handles.plots.figpos{ plot_idx }(1:2) = figpos(1:2);

% ------------------------------------------------------------
% %%%% The checkbox has been unchecked:
% ------------------------------------------------------------
if isempty( handle_fig )
    % > Uncheck checkbox
    % # done
    if handles.plots.fixed_figs==1
        % > Close the figure window:
        delete( handles.plots.fig_h( plot_idx ) );
    end
    % > Remove from database handles.plots.fig_h(fig_idx)-ot nan-ra
    handles.plots.fig_h( plot_idx )=nan;%Update handles.plots.fig_h
    %Update handles structure
    guidata(h_main_fig, handles);
% ------------------------------------------------------------
% %%%% The figure window was closed by [x]
% ------------------------------------------------------------
else%figure [x] -> CloseRequestFcn
    % > Close the figure
    figlist = findobj(0,'Type','figure');%Find out which > figure list
    delete( figlist(1) )%last selected=first on the figure list
    % > Remove from database handles.plots.fig_h
    %if figure handles is in handles.plots.fig_h:
    if isempty( find( handles.plots.fig_h == handle_fig, 1) )~=1
        handles.plots.fig_h( plot_idx )=nan;%Update handles.plots.fig_h
        %Update handles structure
        guidata(h_main_fig, handles);
        % > Uncheck checkbox
        if handles.plots.fixed_figs==1
            if plot_idx==1
                set(handles.ui_plot_fig1_checkbox,         'Value',0);
            end
            if plot_idx==2
                set(handles.ui_plot_fig2_checkbox,          'Value',0);
            end
            if plot_idx==3
                set(handles.ui_plot_fig3_checkbox,     'Value',0);
            end
            if plot_idx==4
                set(handles.ui_plot_fig4_checkbox,       'Value',0);
            end
            if plot_idx==5
                set(handles.ui_plot_fig5_checkbox,     'Value',0);
            end
            if plot_idx==6
                set(handles.ui_plot_fig6_checkbox,      'Value',0);
            end
            if plot_idx==9
                set(handles.ui_plot_fig9_checkbox,     'Value',0);
            end
        end
    end
end
% ------------------------------------------------------------------------
% WARNING: update_fig_closing also called at figure1_CloseRequestFcn !!!
% ------------------------------------------------------------------------



% ------------------------------------------------------------
% Plots: Update plot_figures content for selected file and components:
% ------------------------------------------------------------
function gui_plots_update( hObject, handles, fig_to_plot )
%%% input argument: fig_to_plot
% ==0:   update all (e.g. selction change)
% ==num: update only one!
%DURING FILL PLOT:
% if single file selection && needed seri enabled -> fill
% else (multiple file selection or needed seri not enabled)->fill with n.a.
% -------------------------------------------------------------------------
%%% Get file and component indexes:
% -------------------------------------------------------------------------
%Get selected file index
f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );

%Get selected components:
used_velcomp = [0 0 0];
used_velcomp_names = cell(0);
if get(handles.ui_plot_comp_u_checkbox,'Value')==1
    used_velcomp(1)=1;
    used_velcomp_names{end+1} = 'u';
end
if get(handles.ui_plot_comp_v_checkbox,'Value')==1
    used_velcomp(2)=1;
    used_velcomp_names{end+1} = 'v';
end
if get(handles.ui_plot_comp_w_checkbox,'Value')==1
    used_velcomp(3)=1;
    used_velcomp_names{end+1} = 'w';
end
% -------------------------------------------------------------------------
%%% Statistics box (pos 1 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if (fig_to_plot==0 && get(handles.ui_plot_fig1_checkbox, 'Value')==1) ||...
        fig_to_plot==1
    plot_idx = 1;
    % ---------------------------------------------------------------------
    % FILL TABLE with actual selections content:
    % ---------------------------------------------------------------------
    pos_wi_comp = find( handles.field_props.valu_uvw );%non-zero
    pos_no_comp = find( handles.field_props.valu_uvw == 0 );%zeros
    %%% Case Single file selection
    if length(f_idx)==1
        %Get valid data:
        curr_data = handles.ts( f_idx ).valu_data;
        curr_data( handles.ts( f_idx ).valu_enable==0 ) = nan;
        %separate to comp and non-comp
        data_for_2 = curr_data(pos_wi_comp);
        data2show_1 = curr_data(pos_no_comp)';
    %%% Case Multiple file selection
    else
        data_for_2 =...
            nan(size( handles.field_props.valu_names(pos_wi_comp) ));
        data2show_1 =...
            nan(size( handles.field_props.valu_names(pos_no_comp) ))';
    end
    %%% Reshape data2show_1 to needed sizes:
    data2show_2 = permute(...
        reshape( data_for_2, [ 3,numel(data_for_2)/3 ]), [2 1]);
    set( handles.plots.axs_h(plot_idx,2) , 'Data',     data2show_2);
    set( handles.plots.axs_h(plot_idx,1) , 'Data',     data2show_1);
    % ---------------------------------------------------------------------
    % Modify figure name
    % ---------------------------------------------------------------------
    if length(f_idx)==1
        set( handles.plots.fig_h( plot_idx ),'Name',...
            [handles.fname{f_idx} ' - ' handles.plots.names{ plot_idx }]);
    else
        set( handles.plots.fig_h( plot_idx ),'Name',...
            handles.plots.names{ plot_idx});
    end
end
% -------------------------------------------------------------------------
% %%% QUALITY PLOT (pos 2 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if (fig_to_plot==0 && get(handles.ui_plot_fig2_checkbox, 'Value')==1) ||...
        fig_to_plot==2
    plot_idx = 2;
    % ---------------------------------------------------------------------
    % Clear axes 
    % ---------------------------------------------------------------------
    cla( handles.plots.axs_h( plot_idx,1) )
    cla( handles.plots.axs_h( plot_idx,2) )
    cla( handles.plots.axs_h( plot_idx,3) )
    % ---------------------------------------------------------------------
    % FILL PLOT with actual selections content:
    % ---------------------------------------------------------------------
    axes( handles.plots.axs_h( plot_idx,1) );
    %%% Case Single file selection
    if length(f_idx)==1 &&...
            strcmp( handles.ts( f_idx ).p.dat_props.type, 'usr' )~=1
        hold on
        barh( handles.ts(f_idx).corsnr_dist{1,2},...
            handles.ts(f_idx).corsnr_dist{1,3}(2:end) ,0.75,'c');
        line(handles.ts(f_idx).corsnr_dist{1,4}(2:end),...
            handles.ts(f_idx).corsnr_dist{1,2} ,...
            'Color','b','Marker','.');
        %line( [0 100], [70 70] ,'Color','k','LineStyle',':');
    %%% Case Multiple file selection
    else
        oldlims(1,:) = get( handles.plots.axs_h( plot_idx,1) , 'Xlim');
        oldlims(2,:) = get( handles.plots.axs_h( plot_idx,1) , 'Ylim');
        text( oldlims(1,2)/2,oldlims(2,2)/2,{'not';'available'},...
            'FontSize',18,'HorizontalAlignment','center')
    end
    axes( handles.plots.axs_h( plot_idx,2) );
    %%% Case Single file selection
    % & if err.filt done and at least 2 valid values
    if length(f_idx)==1 && sum( handles.ts( f_idx ).vals )>1 &&...
            strcmp( handles.ts( f_idx ).p.dat_props.type, 'usr' )~=1
        hold on
        %Boxplot:
        pmi = handles.ts( f_idx ).valu_data( 19 );
        p25 = handles.ts( f_idx ).valu_data( 20 );
        pme = handles.ts( f_idx ).valu_data( 21 );
        p75 = handles.ts( f_idx ).valu_data( 22 );
        pma = handles.ts( f_idx ).valu_data( 23 );
        rectangle( 'Position',[65, p25, 20, p75-p25] )
        line([70 80],[pmi pmi],'Color','k','LineStyle','-');
        line([65 85],[pme pme],'Color','k','LineStyle','-','LineWidth',2);
        line([70 80],[pma pma],'Color','k','LineStyle','-');
        line([75 75],[pmi p25],'Color','k','LineStyle','-','LineWidth',1);
        line([75 75],[p75 pma],'Color','k','LineStyle','-','LineWidth',1);
        %70%
        line( [60 90], [70 70] ,'Color','k','LineStyle',':');
    %%% Case Multiple file selection or not enough valid values:
    % - nothing to plot
    end
    axes( handles.plots.axs_h( plot_idx,3) );
    %%% Case Single file selection
    if length(f_idx)==1 &&...
            strcmp( handles.ts( f_idx ).p.dat_props.type, 'usr' )~=1
        hold on
        barh( handles.ts(f_idx).corsnr_dist{3,2},...
            handles.ts(f_idx).corsnr_dist{3,3}(2:end),0.75,'c');
        line(handles.ts(f_idx).corsnr_dist{3,4}(2:end),...
            handles.ts(f_idx).corsnr_dist{3,2} ,...
            'Color','b','Marker','.');
        %line( [0 100], [15 15] ,'Color','k','LineStyle',':');
    %%% Case Multiple file selection
    else
        oldlims(1,:) = get( handles.plots.axs_h( plot_idx,3) , 'Xlim');
        oldlims(2,:) = get( handles.plots.axs_h( plot_idx,3) , 'Ylim');
        text( oldlims(1,2)/2,oldlims(2,2)/2,{'not';'available'},...
            'FontSize',18,'HorizontalAlignment','center')
    end
    % ---------------------------------------------------------------------
    % Modify figure name
    % ---------------------------------------------------------------------
    if length(f_idx)==1
        set( handles.plots.fig_h( plot_idx ),'Name',...
            [handles.fname{f_idx} ' - ' handles.plots.names{ plot_idx,1}]);
    else
        set( handles.plots.fig_h( plot_idx ),'Name',...
            handles.plots.names{ plot_idx,1});
    end
end
% -------------------------------------------------------------------------
% %%% TS Raw vs. (err.filtd & intpoltd) PLOT (pos 3 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if (fig_to_plot==0 && get(handles.ui_plot_fig3_checkbox,'Value')==1) ||...
        fig_to_plot==3
    plot_idx = 3;
    % ---------------------------------------------------------------------
    % Clear axes
    % ---------------------------------------------------------------------
    cla( handles.plots.axs_h( plot_idx,1) )
    cla( handles.plots.axs_h( plot_idx,2) )
    % ---------------------------------------------------------------------
    % FILL PLOT with actual selections content:
    % ---------------------------------------------------------------------
    %%% Raw time-series (where nans -> no values) + bad pos of err.filt
    axes( handles.plots.axs_h( plot_idx,1) );
    %%% Case Single file selection
    if length(f_idx)==1
        hold on
        for cmp_idx=3:-1:1
        if used_velcomp(cmp_idx)==1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %plot ts
            plot(handles.ts( f_idx ).raw_t,...
                 handles.ts( f_idx ).raw_veldata(:,cmp_idx),...
                 handles.plots.cmp_color{cmp_idx})%,'Parent',axes4ts);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %mark invalid posititons:
%             %Spike filter markers
%             errpos = find(handles.ts( f_idx ).efilt_res.errmarker(:,...
%                 strcmpi(handles.ts( f_idx ).efilt_res.eftyp,'VCT'))==1);
%             h_mrk_vct = plot( handles.ts( f_idx ).raw_t( errpos ),...
%                 handles.ts( f_idx ).raw_veldata( errpos ,cmp_idx),...
%                 'ok', 'Parent',handles.plots.axs_h( plot_idx,1));
%             errpos = find(handles.ts( f_idx ).efilt_res.errmarker(:,...
%                 strcmpi(handles.ts( f_idx ).efilt_res.eftyp,'PST'))==1);
%             h_mrk_pst = plot( handles.ts( f_idx ).raw_t( errpos ),...
%                 handles.ts( f_idx ).raw_veldata( errpos ,cmp_idx),...
%                 'ok', 'Parent',handles.plots.axs_h( plot_idx,1));
%             %Cor-filter markers
%             errpos = find(handles.ts( f_idx ).efilt_res.errmarker(:,...
%                 strcmpi(handles.ts( f_idx ).efilt_res.eftyp,'cormin'))==1);
%             h_mrk_cmi = plot( handles.ts( f_idx ).raw_t( errpos ),...
%                 handles.ts( f_idx ).raw_veldata( errpos ,cmp_idx),...
%                 'xk', 'Parent',handles.plots.axs_h( plot_idx,1));
%             errpos = find(handles.ts( f_idx ).efilt_res.errmarker(:,...
%                 strcmpi(handles.ts( f_idx ).efilt_res.eftyp,'coravg'))==1);
%             h_mrk_cma = plot( handles.ts( f_idx ).raw_t( errpos ),...
%                 handles.ts( f_idx ).raw_veldata( errpos ,cmp_idx),...
%                 'xk', 'Parent',handles.plots.axs_h( plot_idx,1));
%             if handles.plots.axs_legend_on(plot_idx)==1
%                 set(get(get(h_mrk_vct,'Annotation'),'LegendInformation'),...
%                     'IconDisplayStyle','off');
%                 set(get(get(h_mrk_pst,'Annotation'),'LegendInformation'),...
%                     'IconDisplayStyle','off');
%                 set(get(get(h_mrk_cmi,'Annotation'),'LegendInformation'),...
%                     'IconDisplayStyle','off');
%                 set(get(get(h_mrk_cma,'Annotation'),'LegendInformation'),...
%                     'IconDisplayStyle','off');
%             end
        end
        end
        if handles.plots.axlims_to_use( plot_idx,1)==1
            set(handles.plots.axs_h( plot_idx,1),...
                'Xlim',handles.plots.axlims_usrdef{ plot_idx,1 }(1,:),...
                'Ylim',handles.plots.axlims_usrdef{ plot_idx,1 }(2,:));
        end
        %Legend
        if handles.plots.axs_legend_on(plot_idx)==1
            legend( handles.plots.axs_h( plot_idx,1),...
                used_velcomp_names,...
                'FontName',handles.plots.ft_name,...
                'FontSize',handles.plots.siz.ft_px*handles.plots.siz.ft_px2pt);%Fontunits: points
            %'Location','Best')
        end
    %%% Case Multiple file selection
    else
        if handles.plots.axlims_to_use( plot_idx,1)==1
            currlims = handles.plots.axlims_usrdef{ plot_idx,1 };
        else
            [ currlims ] = fcn_sokoray_fit_lims2data(...
                handles.plots.axlims_minmax{ plot_idx, 1}(:,1),...
                handles.plots.axlims_minmax{ plot_idx, 1}(:,2));
        end
        text( (currlims(1,2)+currlims(1,1))/2,...
            (currlims(2,2)+currlims(2,1))/2,...
            {'not';'available'},...
            'FontSize',18,'HorizontalAlignment','center')
        set(handles.plots.axs_h( plot_idx ,1),...
            'Xlim',currlims(1,:), 'Ylim',currlims(2,:));
    end
    % ---------------------------------------------------------------------
    %%% Err-filtered (and interpolated) time series:
    axes( handles.plots.axs_h( plot_idx,2) );
    %%% Case Single file selection && if needed seri enabled:
    if length(f_idx)==1 && handles.ts( f_idx ).seri_enable(1)==1
        hold on
        for cmp_idx=3:-1:1%go through the components (1=u, 2=v, 3=w)
        if used_velcomp( cmp_idx )==1
            %Showing original and interpolated data
%             plot(handles.ts( f_idx ).seri_x.t,...
%                 handles.ts( f_idx ).seri_data{1}(:,cmp_idx) + ...
%                 handles.ts( f_idx ).intl.mean_correctn( cmp_idx ) + ...
%                 handles.ts( f_idx ).valu_data( cmp_idx ),...
%                 handles.plots.cmp_color{ cmp_idx });%,'Parent',axes4tsok);
            %Showing original data with gaps: nans at gaps needed!
%             tic
            temp_uvw = nan(length( handles.ts( f_idx ).raw_t ),1);
            temp_uvw( handles.ts( f_idx ).vals ) =...
                handles.ts( f_idx ).raw_veldata( handles.ts( f_idx ).vals, cmp_idx );
%             toc
            plot( handles.ts( f_idx ).raw_t(handles.ts( f_idx ).intl.used_ts_seg),...
                temp_uvw(handles.ts( f_idx ).intl.used_ts_seg),...
                handles.plots.cmp_color{ cmp_idx } );%,'Parent',axes4tsok);
        end       
        end
        if handles.plots.axlims_to_use( plot_idx,2)==1
            set(handles.plots.axs_h( plot_idx,2),...
                'Xlim',handles.plots.axlims_usrdef{ plot_idx,2 }(1,:),...
                'Ylim',handles.plots.axlims_usrdef{ plot_idx,2 }(2,:));
        end
        %Legend
        if handles.plots.axs_legend_on(plot_idx)==1
            legend( handles.plots.axs_h( plot_idx,2),...
                used_velcomp_names,...
                'FontName',handles.plots.ft_name,...
                'FontSize',handles.plots.siz.ft_px*handles.plots.siz.ft_px2pt);%Fontunits: points
            %'Location','Best')
        end
    %%% Case Multiple file selection or if needed seri NOT enabled
    else
        if handles.plots.axlims_to_use( plot_idx,2)==1
            currlims = handles.plots.axlims_usrdef{ plot_idx,2 };
        else
            [ currlims ] = fcn_sokoray_fit_lims2data(...
                handles.plots.axlims_minmax{ plot_idx, 2}(:,1),...
                handles.plots.axlims_minmax{ plot_idx, 2}(:,2));
        end
        text( (currlims(1,2)+currlims(1,1))/2,...
            (currlims(2,2)+currlims(2,1))/2,...
            {'not';'available'},...
            'FontSize',18,'HorizontalAlignment','center')
        set(handles.plots.axs_h( plot_idx ,2),...
            'Xlim',currlims(1,:), 'Ylim',currlims(2,:));
    end
    % ---------------------------------------------------------------------
    % Link axes
    % ---------------------------------------------------------------------
    linkaxes([handles.plots.axs_h( plot_idx,1);...
        handles.plots.axs_h( plot_idx,2)],  'xy');
    % ---------------------------------------------------------------------
    % Modify figure name
    % ---------------------------------------------------------------------
    if length(f_idx)==1
        set( handles.plots.fig_h( plot_idx ),'Name',...
            [handles.fname{f_idx} ' - ' handles.plots.names{ plot_idx,1}]);
    else
        set( handles.plots.fig_h( plot_idx ),'Name',...
            handles.plots.names{ plot_idx,1});
    end
end
% -------------------------------------------------------------------------
% %%% Velocity histograms (pos4 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if (fig_to_plot==0 && get(handles.ui_plot_fig4_checkbox,'Value')==1) ||...
        fig_to_plot==4
    plot_idx = 4;
    % ---------------------------------------------------------------------
    % Clear axes
    % ---------------------------------------------------------------------
    cla( handles.plots.axs_h( plot_idx,1) )
    cla( handles.plots.axs_h( plot_idx,2) )
    cla( handles.plots.axs_h( plot_idx,3) )
    % ---------------------------------------------------------------------
    % FILL PLOT with actual selections content:
    % ---------------------------------------------------------------------
    for cmp_idx=3:-1:1
        axes( handles.plots.axs_h( plot_idx,cmp_idx) );
        %%% Case Single file selection
        if length(f_idx)==1 && used_velcomp(cmp_idx)==1
            hold on
            bar( handles.ts( f_idx ).hist.vel_of_bin(:,cmp_idx),...
                handles.ts( f_idx ).hist.p_per_bin(:,cmp_idx) ,1,...
                handles.plots.cmp_color{cmp_idx});
            if handles.plots.axlims_to_use( plot_idx,cmp_idx )==1
                set(handles.plots.axs_h( plot_idx,cmp_idx),...
                    'Xlim',handles.plots.axlims_usrdef{ plot_idx,cmp_idx }(1,:),...
                    'Ylim',handles.plots.axlims_usrdef{ plot_idx,cmp_idx }(2,:));
            end
        %%% Case Multiple file selection
        else
            if handles.plots.axlims_to_use( plot_idx,cmp_idx )==1
                currlims = handles.plots.axlims_usrdef{ plot_idx,cmp_idx };
            else
                [ currlims ] = fcn_sokoray_fit_lims2data(...
                    handles.plots.axlims_minmax{ plot_idx, cmp_idx}(:,1),...
                    handles.plots.axlims_minmax{ plot_idx, cmp_idx}(:,2));
            end
            text( (currlims(1,2)+currlims(1,1))/2,...
                (currlims(2,2)+currlims(2,1))/2,...
                {'not';'available'},...
                'FontSize',18,'HorizontalAlignment','center')
            set(handles.plots.axs_h( plot_idx ,cmp_idx),...
                'Xlim',currlims(1,:), 'Ylim',currlims(2,:));
        end
    end
    % ---------------------------------------------------------------------
    % Modify figure name
    % ---------------------------------------------------------------------
    if length(f_idx)==1
        set( handles.plots.fig_h( plot_idx ),'Name',...
            [handles.fname{f_idx} ' - ' handles.plots.names{ plot_idx,1}]);
    else
        set( handles.plots.fig_h( plot_idx ),'Name',...
            handles.plots.names{ plot_idx,1});
    end
end
% -------------------------------------------------------------------------
% %%% Spectrum PLOT (pos 5 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if (fig_to_plot==0 && get(handles.ui_plot_fig5_checkbox,'Value')==1) ||...
        fig_to_plot==5
    plot_idx = 5;
    % ---------------------------------------------------------------------
    % Clear axes 
    % ---------------------------------------------------------------------
    %Get old position - NOT for Spectra
    %Clear axes:
    cla( handles.plots.axs_h( plot_idx,1) )
    %Get limits:
    powlims(1,:) = handles.plots.axlims_usrdef{ plot_idx,1 }(1,:);
    powlims(2,:) = handles.plots.axlims_usrdef{ plot_idx,1 }(2,:);
    currlims=10.^powlims;
    % ---------------------------------------------------------------------
    % FILL PLOT with actual selections content:
    % ---------------------------------------------------------------------
    % Plot spectra
    axes( handles.plots.axs_h( plot_idx,1) );%no hold on, otherwise linear
    %%% Case Single file selection && if needed seri enabled:
    if length(f_idx)==1 &&...
            handles.ts( f_idx ).seri_enable(8)==1
        hold on
        for cmp_idx=3:-1:1%go through the components (1=u, 2=v, 3=w=
        if used_velcomp( cmp_idx )==1 &&...
                handles.ts( f_idx ).p.errfilt.used_velcomp( cmp_idx )==1
            loglog(handles.ts( f_idx ).seri_x.f,...
                handles.ts( f_idx ).seri_data{8}(:,cmp_idx),...
                handles.plots.cmp_color{ cmp_idx },...
                'Parent',handles.plots.axs_h( plot_idx,1));
            %hold on
        end
        end
        if handles.viper_props.userdef.plot_line53 == 1
            loglog( handles.viper_props.userdef.plot_line53_x,...
                    handles.viper_props.userdef.plot_line53_y ,...
                    '--k');%'--','Color',[0.5,0.5,0.5]);
            %End of line
            text( 10^(log10(handles.viper_props.userdef.plot_line53_x(2))+0),...
                  10^(log10(handles.viper_props.userdef.plot_line53_y(2))+0.25),...
                {'-5/3'},'Color',[0.5,0.5,0.5],...
                'FontName',handles.plots.ft_name,...
                'FontSize',handles.plots.siz.ft_px*handles.plots.siz.ft_px2pt)
        end
        if handles.viper_props.userdef.plot_line3 == 1
            loglog( handles.viper_props.userdef.plot_line3_x,...
                    handles.viper_props.userdef.plot_line3_y ,...
                    '--k');%'--','Color',[0.5,0.5,0.5]);
            %End of line
            text( 10^(log10(handles.viper_props.userdef.plot_line3_x(1))+0.25),...
                  10^(log10(handles.viper_props.userdef.plot_line3_y(1))+0),...
                {'-3'},'Color',[0.5,0.5,0.5],...
                'FontName',handles.plots.ft_name,...
                'FontSize',handles.plots.siz.ft_px*handles.plots.siz.ft_px2pt)
        end
        if handles.viper_props.userdef.plot_line1 == 1
            loglog( handles.viper_props.userdef.plot_line1_x,...
                    handles.viper_props.userdef.plot_line1_y ,...
                    '--k');%'--','Color',[0.5,0.5,0.5]);
            %End of line
            text( 10^(log10(handles.viper_props.userdef.plot_line1_x(2))+0),...
                  10^(log10(handles.viper_props.userdef.plot_line1_y(2))+0.25),...
                {'-1'},'Color',[0.5,0.5,0.5],...
                'FontName',handles.plots.ft_name,...
                'FontSize',handles.plots.siz.ft_px*handles.plots.siz.ft_px2pt)
        end
        set(handles.plots.axs_h( plot_idx,1),...
            'GridLineStyle',':','XGrid','on','YGrid','on',...
            'XMinorTick','on','YMinorTick','on',...
            'XLim',currlims(1,:),'YLim',currlims(2,:),...
            'XMinorGrid', 'off', 'YMinorGrid', 'off');
        %Legend
        if handles.plots.axs_legend_on(plot_idx)==1
            legend( handles.plots.axs_h( plot_idx,1),...
                used_velcomp_names,...
                'FontName',handles.plots.ft_name,...
                'FontSize',handles.plots.siz.ft_px*handles.plots.siz.ft_px2pt);%Fontunits: points
            %'Location','Best')
        end
%         xlabel( handles.plots.axs_h( plot_idx,1) ,...
%             handles.plots.axs_label_x{ plot_idx,1},...
%             'FontName',handles.plots.ft_name,...
%             'FontSize',handles.plots.siz.ft_px);
%         ylabel( handles.plots.axs_h( plot_idx,1) ,...
%             handles.plots.axs_label_y{ plot_idx,1},...
%             'FontName',handles.plots.ft_name,...
%             'FontSize',handles.plots.siz.ft_px);
    %%% Case Multiple file selection or if needed seri NOT enabled
    else
        currlims(1,:) = get( handles.plots.axs_h(5,1) , 'Xlim');
        currlims(2,:) = get( handles.plots.axs_h(5,1) , 'Ylim');
        text( (currlims(1,2)+currlims(1,1))/2,...
            (currlims(2,2)+currlims(2,1))/2,...
            {'not';'available'},...
            'FontSize',18,'HorizontalAlignment','center')
        set(handles.plots.axs_h( plot_idx ,1),...
            'Xlim',currlims(1,:), 'Ylim',currlims(2,:));
    end
    % ---------------------------------------------------------------------
    % Modify figure name
    % ---------------------------------------------------------------------
    if length(f_idx)==1
        set( handles.plots.fig_h( plot_idx ),'Name',...
            [handles.fname{f_idx} ' - ' handles.plots.names{ plot_idx,1}]);
    else
        set( handles.plots.fig_h( plot_idx ),'Name',...
            handles.plots.names{ plot_idx,1});
    end
end
% -------------------------------------------------------------------------
% %%% Cumulative mean PLOT (pos 6 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if (fig_to_plot==0 && get(handles.ui_plot_fig6_checkbox,'Value')==1) ||...
        fig_to_plot==6
    plot_idx = 6;
    % ---------------------------------------------------------------------
    % Clear axes
    % ---------------------------------------------------------------------
    cla( handles.plots.axs_h( plot_idx,1) )
    % ---------------------------------------------------------------------
    % FILL PLOT with actual selections content:
    % ---------------------------------------------------------------------
    % Plot cummean
    axes( handles.plots.axs_h( plot_idx,1) );
    %%% Case Single file selection && if needed seri enabled:
    if length(f_idx)==1 && handles.ts( f_idx ).seri_enable(4)==1
        hold on
        for cmp_idx=3:-1:1%go through the components (1=u, 2=v, 3=w)
        if used_velcomp( cmp_idx )==1
            plot(handles.ts( f_idx ).seri_x.t,...
                handles.ts( f_idx ).seri_data{4}(:,cmp_idx),...
                handles.plots.cmp_color{ cmp_idx });%,'Parent',axes4tsok);
        end
        end
        if handles.plots.axlims_to_use( plot_idx,1)==1
            set(handles.plots.axs_h( plot_idx,1),...
                'Xlim',handles.plots.axlims_usrdef{ plot_idx,1 }(1,:),...
                'Ylim',handles.plots.axlims_usrdef{ plot_idx,1 }(2,:));
        end
    %%% Case Multiple file selection or if needed seri NOT enabled
    else
        if handles.plots.axlims_to_use( plot_idx,1)==1
            currlims = handles.plots.axlims_usrdef{ plot_idx,1 };
        else
            [ currlims ] = fcn_sokoray_fit_lims2data(...
                handles.plots.axlims_minmax{ plot_idx, 1}(:,1),...
                handles.plots.axlims_minmax{ plot_idx, 1}(:,2));
        end
        text( (currlims(1,2)+currlims(1,1))/2,...
            (currlims(2,2)+currlims(2,1))/2,...
            {'not';'available'},...
            'FontSize',18,'HorizontalAlignment','center')
        set(handles.plots.axs_h( plot_idx ,1),...
            'Xlim',currlims(1,:), 'Ylim',currlims(2,:));
    end
    %Legend
    if handles.plots.axs_legend_on(plot_idx)==1
        legend( handles.plots.axs_h( plot_idx,1),...
            used_velcomp_names,...
            'FontName',handles.plots.ft_name,...
            'FontSize',handles.plots.siz.ft_px*handles.plots.siz.ft_px2pt);%Fontunits: points
        %'Location','Best')
    end
    % ---------------------------------------------------------------------
    % Modify figure name
    % ---------------------------------------------------------------------
    if length(f_idx)==1
        set( handles.plots.fig_h( plot_idx ),'Name',...
            [handles.fname{f_idx} ' - ' handles.plots.names{ plot_idx,1}]);
    else
        set( handles.plots.fig_h( plot_idx ),'Name',...
            handles.plots.names{ plot_idx,1});
    end
end
% -------------------------------------------------------------------------
% %%% Autocorr function (pos 9 in handles.plots.fig_h)
% -------------------------------------------------------------------------
if (fig_to_plot==0 && get(handles.ui_plot_fig9_checkbox,'Value')==1) ||...
        fig_to_plot==9
    plot_idx = 9;
    % ---------------------------------------------------------------------
    % Clear axes
    % ---------------------------------------------------------------------
    cla( handles.plots.axs_h( plot_idx,1) )
    % ---------------------------------------------------------------------
    % FILL PLOT with actual selections content:
    % ---------------------------------------------------------------------
    % Plot autocorr
    axes( handles.plots.axs_h( plot_idx,1) );
    %%% Case Single file selection && if needed seri enabled:
    if length(f_idx)==1 &&...
            handles.ts( f_idx ).seri_enable(3)==1
        hold on
        for cmp_idx=3:-1:1%go through the components (1=u, 2=v, 3=w)
        if used_velcomp( cmp_idx )==1
            plot(handles.ts( f_idx ).seri_x.t,...
                handles.ts( f_idx ).seri_data{3}(:,cmp_idx),...
                handles.plots.cmp_color{ cmp_idx });%,'Parent',axes4tsok);
        end
        end
        if handles.plots.axlims_to_use( plot_idx,1)==1
            set(handles.plots.axs_h( plot_idx,1),...
                'Xlim',handles.plots.axlims_usrdef{ plot_idx,1 }(1,:),...
                'Ylim',handles.plots.axlims_usrdef{ plot_idx,1 }(2,:));
        end
    %%% Case Multiple file selection or if needed seri NOT enabled
    else
        if handles.plots.axlims_to_use( plot_idx,1)==1
            currlims = handles.plots.axlims_usrdef{ plot_idx,1 };
        else
            [ currlims ] = fcn_sokoray_fit_lims2data(...
                handles.plots.axlims_minmax{ plot_idx, 1}(:,1),...
                handles.plots.axlims_minmax{ plot_idx, 1}(:,2));
        end
        text( (currlims(1,2)+currlims(1,1))/2,...
            (currlims(2,2)+currlims(2,1))/2,...
            {'not';'available'},...
            'FontSize',18,'HorizontalAlignment','center')
        set(handles.plots.axs_h( plot_idx ,1),...
            'Xlim',currlims(1,:), 'Ylim',currlims(2,:));
    end
    %Legend
    if handles.plots.axs_legend_on(plot_idx)==1
        legend( handles.plots.axs_h( plot_idx,1),...
            used_velcomp_names,...
            'FontName',handles.plots.ft_name,...
            'FontSize',handles.plots.siz.ft_px*handles.plots.siz.ft_px2pt);%Fontunits: points
        %'Location','Best')
    end
    % ---------------------------------------------------------------------
    % Modify figure name
    % ---------------------------------------------------------------------
    if length(f_idx)==1
        set( handles.plots.fig_h( plot_idx ),'Name',...
            [handles.fname{f_idx} ' - ' handles.plots.names{ plot_idx,1}]);
    else
        set( handles.plots.fig_h( plot_idx ),'Name',...
            handles.plots.names{ plot_idx,1});
    end
end
% -------------------------------------------------------------------------
% %%% DVL
% -------------------------------------------------------------------------
if fig_to_plot==0 && isfield( handles, 'dvl')==1
    [ ~, ~ ] = viper_dvl('dvl_update_plots', handles, {[]} );
end




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 
function chapter_04_dialogs_and_functions_for_plot_setup
% 
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------
% Plot properties: Calculates min/max values of plot axes :
% ------------------------------------------------------------
function [ data_info,axlims_minmax ]= gui_plots_calc_axlims_minmax(handles)
data_info.raw_mins = nan(1,3);
data_info.raw_maxs = nan(1,3);
data_info.seri_mins = ...
    nan(size( handles.ts_defaults.seri_data, 1),3);
data_info.seri_maxs = ...
    nan(size( handles.ts_defaults.seri_data, 1),3);
data_info.t_minmax = nan(1,2);
data_info.d_minmax = nan(1,2);
data_info.f_minmax = nan(1,2);
data_info.hisp_max = nan(1,3);
data_info.hisv_mins = nan(1,3);
data_info.hisv_maxs = nan(1,3);
% ------------------------------------------------------------
% Determine data minima and maxima over all files for "used_velcomp"
% for unused velcomp -> NaN
for f_idx = 1:length(handles.ts)
    cmp_idx = find( handles.ts( f_idx ).p.errfilt.used_velcomp==1 )';
    %Raw u,v,w-minima compared to other files
    data_info.raw_mins(cmp_idx) = min(...
        data_info.raw_mins(cmp_idx),...
        min(handles.ts( f_idx ).raw_veldata(:,cmp_idx),[],1));
    %Raw u,v,w-maxima compared to other files
    data_info.raw_maxs(cmp_idx) = max(...
        data_info.raw_maxs(cmp_idx),...
        max(handles.ts( f_idx ).raw_veldata(:,cmp_idx),[],1));
    %Seri
    for seri_idx = 1:size( handles.ts( f_idx ).seri_data, 1)
    if handles.ts( f_idx ).seri_enable( seri_idx )==1
        if seri_idx==1
            %Seri u,v,w-minima compared to other files
            data_info.seri_mins( seri_idx, cmp_idx ) = min(...
                data_info.seri_mins( seri_idx, cmp_idx ), ...
                min(...
                handles.ts( f_idx ).seri_data{seri_idx}(:,cmp_idx),[],1) +...
                handles.ts( f_idx ).valu_data(cmp_idx));
            %data_info.seri_mins( seri_idx, : )
            %Seri u,v,w-maxima compared to other files
            data_info.seri_maxs( seri_idx, cmp_idx ) = max(...
                data_info.seri_maxs( seri_idx, cmp_idx ), ...
                max(...
                handles.ts( f_idx ).seri_data{seri_idx}(:,cmp_idx),[],1) +...
                handles.ts( f_idx ).valu_data(cmp_idx));
            %data_info.seri_maxs( seri_idx, : )
        else
            %Seri u,v,w-minima compared to other files
            data_info.seri_mins( seri_idx, cmp_idx ) = min(...
                data_info.seri_mins( seri_idx, cmp_idx ), ...
                min(...
                handles.ts( f_idx ).seri_data{seri_idx}(:,cmp_idx),[],1));
            %data_info.seri_mins( seri_idx, : )
            %Seri u,v,w-maxima compared to other files
            data_info.seri_maxs( seri_idx, cmp_idx ) = max(...
                data_info.seri_maxs( seri_idx, cmp_idx ), ...
                max(...
                handles.ts( f_idx ).seri_data{seri_idx}(:,cmp_idx),[],1));
            %data_info.seri_maxs( seri_idx, : )
        end
    end
    end
    %non-value-variables (mostly x-variables):
    data_info.t_minmax = [...
        min( data_info.t_minmax(1), min(handles.ts( f_idx ).raw_t) ),...
        max( data_info.t_minmax(2), max(handles.ts( f_idx ).raw_t) )];
    data_info.d_minmax = max( data_info.d_minmax,...
        [ 0 max(handles.ts( f_idx ).seri_x.tlag) ] );
    data_info.f_minmax(1) = min(...
        data_info.f_minmax(1),...
        handles.ts( f_idx ).seri_x.f(1) );
    %pwelch produces lower freqs than this: (Why? TODO)
%         1/(length( handles.ts(f_idx).seri_x.t )/...
%         handles.ts(f_idx).p.dat_props.frq) );%min freq
    data_info.f_minmax(2) = max(...
        data_info.f_minmax(2),...
        max(handles.ts( f_idx ).seri_x.f) );%max freq
    %Histogram u.v.w-p-minima/maxima:
    data_info.hisp_max(cmp_idx) = max( data_info.hisp_max(cmp_idx),...
        max(handles.ts( f_idx ).hist.p_per_bin(:,cmp_idx), [], 1) );
    %Hist u,v,w-minima/maxima:
    data_info.hisv_mins(cmp_idx) = min( data_info.hisv_mins(cmp_idx),...
        min( handles.ts( f_idx ).hist.vel_of_bin(:,cmp_idx), [], 1) );
    data_info.hisv_maxs(cmp_idx) = max( data_info.hisv_maxs(cmp_idx),...
        max( handles.ts( f_idx ).hist.vel_of_bin(:,cmp_idx), [], 1) );
end
% ------------------------------------------------------------
% Update plot axlims_minmax:%use old as init -> no nans added
axlims_minmax = handles.plots.axlims_usrdef;
for p_i=1:size( handles.plots.axs_h, 1)%Plot-types
%if notadjustable: p_i==1, p_i==2, take old value:
if handles.plots.adjustable_axlims(p_i)==0
    axlims_minmax(p_i,:) = handles.plots.axlims_usrdef(p_i,:);
else%adjustables:
    if p_i==3%ts: limits of axs1==axs2
        %Get idx of seri used as basis:
        seri_idx = handles.plots.reference_seri( p_i );
        axlims_minmax{p_i,1} = [...
            data_info.t_minmax;...
            min([data_info.raw_mins, data_info.seri_mins( seri_idx,: )]),...
            max([data_info.raw_maxs, data_info.seri_maxs( seri_idx,: )]) ];
        axlims_minmax{p_i,2} =...
            axlims_minmax{p_i,1};
        %if minmax of a velcomp is nan -> disappears through using min/max
    end
    if p_i==4%hist
        %if minmax of a velcomp is nan -> it results in a nan min/max value
        %therefore: use next non-nan for such values!
        %1. calc non-nans:
        cmp_val = find(isnan( data_info.hisp_max.*...
            data_info.hisv_mins.*data_info.hisv_maxs)==0);
        for i_c=1:length( cmp_val )
            axlims_minmax{p_i,cmp_val(i_c)} = [...
                data_info.hisv_mins(cmp_val(i_c)),...
                data_info.hisv_maxs(cmp_val(i_c));...
                0, data_info.hisp_max(cmp_val(i_c)) ];
        end
        %2. for nans use the first non-nan:
        cmp_nan = find(isnan( data_info.hisp_max.*...
            data_info.hisv_mins.*data_info.hisv_maxs)==1);
        for i_c=1:length( cmp_nan )
            axlims_minmax{p_i,cmp_nan(i_c)} = ...
                axlims_minmax{p_i,cmp_val(1)};
        end
    end
    if p_i==5%spectra
        %Get idx of seri used as basis:
        %data_info.f_minmax
        seri_idx = handles.plots.reference_seri( p_i );
        axlims_minmax{p_i,1} = [...
            floor(log10(data_info.f_minmax(1))),...
            ceil( log10(data_info.f_minmax(2)));...
            floor(log10(min(data_info.seri_mins( seri_idx, : )))),...
            ceil( log10(max(data_info.seri_maxs( seri_idx, : )))) ];
        %if minmax of a velcomp is nan -> disappears through using min/max
    end
    if p_i==6%cummean
        %Get idx of seri used as basis:
        seri_idx = handles.plots.reference_seri( p_i );
        axlims_minmax{p_i,1} = [...
            data_info.t_minmax;...
            min(data_info.seri_mins( seri_idx, : )),...
            max(data_info.seri_maxs( seri_idx, : )) ];
        %if minmax of a velcomp is nan -> disappears through using min/max
    end
    if p_i==9%cummean
        %Get idx of seri used as basis:
        seri_idx = handles.plots.reference_seri( p_i );
        axlims_minmax{p_i,1} = [...
            data_info.d_minmax;...
            min(data_info.seri_mins( seri_idx, : )),...
            max(data_info.seri_maxs( seri_idx, : )) ];
        %if minmax of a velcomp is nan -> disappears through using min/max
    end
end
end










% ------------------------------------------------------------
% Plot-limits-box: callbacks
% ------------------------------------------------------------
function guiexec_plotlimbox_callbacks( src, eventdata, ...
    h_main_fig, curr_uicontrol)
handles = guidata(h_main_fig);
p_val = find( handles.plots.adjustable_axlims );
%--------------------------------------------------------------------------
%radiobutton change -> update editboxes:
update_needed = 0;
p_sel = find( handles.plotlimbox.p_rb == curr_uicontrol );
if isempty( p_sel )==0
    %set other rb-s to off
    set( handles.plotlimbox.p_rb( p_val ), 'Value',0);
    set( handles.plotlimbox.p_rb( p_sel ), 'Value',1);
    %Change axes selection to 1:
    set( handles.plotlimbox.a_rb( : ), 'Value',0);
    set( handles.plotlimbox.a_rb( 1 ), 'Value',1);
    %Set apply button enable:
    set( handles.plotlimbox.apply, 'Enable','Off');
    %Set axes radiobutton selectability:
    for a_i=1:size( handles.plots.axs_h, 2)
        set( handles.plotlimbox.a_rb( a_i ),...
            'String', handles.plots.axs_names{ p_sel, a_i } )
        if isempty( handles.plots.axlims_usrdef{ p_sel,a_i } )==0
            set( handles.plotlimbox.a_rb( a_i ), 'Enable','On');
        else
            set( handles.plotlimbox.a_rb( a_i ), 'Enable','Off');
        end
    end
    %Set checkbox selectability:
    if p_sel==5
        set( handles.plotlimbox.check, 'Enable','Off');
    elseif handles.plots.adjustable_axlims( p_sel )==1
        set( handles.plotlimbox.check, 'Enable','On');
    else
        set( handles.plotlimbox.check, 'Enable','Off');
    end
    
    update_needed = 1;
end
a_sel = find( handles.plotlimbox.a_rb == curr_uicontrol );
if isempty( a_sel )==0
    %set other rb-s to off
    set( handles.plotlimbox.a_rb( : ), 'Value',0);
    set( handles.plotlimbox.a_rb( a_sel ), 'Value',1);
    %Set apply button enable:
    set( handles.plotlimbox.apply, 'Enable','Off');
    update_needed = 1;
end
%--------------------------------------------------------------------------
%Final selection:
%do not use logicals!
curr_p_rb = find(cell2mat( get(handles.plotlimbox.p_rb( p_val ),'Value')));
p_idx = find(...
    handles.plotlimbox.p_rb == handles.plotlimbox.p_rb(p_val(curr_p_rb)));
curr_a_rb = find(cell2mat( get(handles.plotlimbox.a_rb( : ),'Value')));
a_idx = find(...
    handles.plotlimbox.a_rb == handles.plotlimbox.a_rb(curr_a_rb));
if update_needed == 1
    %load limits:
    set( handles.plotlimbox.al_ed(1),...
        'String', num2str(handles.plots.axlims_usrdef{ p_idx,a_idx }(1,1)),...
        'BackgroundColor',[0.941,0.941,0.941]);
    set( handles.plotlimbox.al_ed(2),...
        'String', num2str(handles.plots.axlims_usrdef{ p_idx,a_idx }(1,2)),...
        'BackgroundColor',[0.941,0.941,0.941]);
    set( handles.plotlimbox.al_ed(3),...
        'String', num2str(handles.plots.axlims_usrdef{ p_idx,a_idx }(2,1)),...
        'BackgroundColor',[0.941,0.941,0.941]);
    set( handles.plotlimbox.al_ed(4),...
        'String', num2str(handles.plots.axlims_usrdef{ p_idx,a_idx }(2,2)),...
        'BackgroundColor',[0.941,0.941,0.941]);
    %load limits_to_use value and set checkbox bgcolor:
    set( handles.plotlimbox.check,...
        'Value', handles.plots.axlims_to_use( p_idx,a_idx ),...
        'BackgroundColor',[0.941,0.941,0.941]);
    %Disable applying
    set( handles.plotlimbox.apply,'Enable','Off',...
        'BackgroundColor',[0.941,0.941,0.941]);
end
%--------------------------------------------------------------------------
%editbox changes:
if isempty(find(handles.plotlimbox.al_ed == curr_uicontrol,1))==0
    %Set editbox bgcolor:
    set( curr_uicontrol,...
        'BackgroundColor',[1,0.6,0.6]);
    %Enable applying
    set( handles.plotlimbox.apply,'Enable','On','BackgroundColor',[1,0.6,0.6])
end
%--------------------------------------------------------------------------
%pushbuttons
%Calculates common limits of all loaded files using all components
if handles.plotlimbox.collect == curr_uicontrol
    %load values to editboxes:
    set( handles.plotlimbox.al_ed(1),'String',...
        num2str(handles.plots.axlims_minmax{p_idx,a_idx}(1,1)),...
        'BackgroundColor',[1,0.6,0.6]);
    set( handles.plotlimbox.al_ed(2),'String',...
        num2str(handles.plots.axlims_minmax{p_idx,a_idx}(1,2)),...
        'BackgroundColor',[1,0.6,0.6]);
    set( handles.plotlimbox.al_ed(3),'String',...
        num2str(handles.plots.axlims_minmax{p_idx,a_idx}(2,1)),...
        'BackgroundColor',[1,0.6,0.6]);
    set( handles.plotlimbox.al_ed(4),'String',...
        num2str(handles.plots.axlims_minmax{p_idx,a_idx}(2,2)),...
        'BackgroundColor',[1,0.6,0.6]);
    %Enable applying
    set( handles.plotlimbox.apply,'Enable','On','BackgroundColor',[1,0.6,0.6])
end
%--------------------------------------------------------------------------
%checkbox changes:
if isempty(find(handles.plotlimbox.check == curr_uicontrol,1))==0
    %Set editbox bgcolor:
    set( curr_uicontrol,...
        'BackgroundColor',[1,0.6,0.6]);
    %Enable applying
    set( handles.plotlimbox.apply,'Enable','On','BackgroundColor',[1,0.6,0.6])
    %Display warning:
    if get( curr_uicontrol, 'Value')==0
        warndlg(...
        '''Using automatic limits'' is applied at next plot creation!',...
        'Warning','modal')
    end
end
%--------------------------------------------------------------------------
%Store values from editboxes and checkbox setting
if handles.plotlimbox.apply == curr_uicontrol
    input = str2double(...
        regexprep( get( handles.plotlimbox.al_ed,'String'), ',', '.') )';
    %Check inputs:
    if sum( isnan( input ) )~=0
        errordlg('At least one input invalid!','Error');
        return
    end
    lims = reshape(input,[2,2])';
    if lims(1,1)>lims(1,2)
        errordlg('X_max has to be larger than X_min!','Error');
        return
    end
    if lims(2,1)>lims(2,2)
        errordlg('Y_max has to be larger than Y_min!','Error');
        return
    end
    %Store inputs:
    handles.plots.axlims_usrdef{ p_idx,a_idx } = reshape(input,[2,2])';
    handles.plots.axlims_to_use( p_idx,a_idx ) =...
        get( handles.plotlimbox.check, 'Value' );
    %Plot-specific issues:
    if p_idx==3%if 1 changed -> adjust 2 and the other way round
        if a_idx==1
            handles.plots.axlims_usrdef{ p_idx,2 } =...
                handles.plots.axlims_usrdef{ p_idx,a_idx };
            handles.plots.axlims_to_use( p_idx,2 ) =...
                handles.plots.axlims_to_use( p_idx,a_idx );
        end
        if a_idx==2
            handles.plots.axlims_usrdef{ p_idx,1 } =...
                handles.plots.axlims_usrdef{ p_idx,a_idx };
            handles.plots.axlims_to_use( p_idx,1 ) =...
                handles.plots.axlims_to_use( p_idx,a_idx );
        end
    end
    % Update handles structure
    guidata( h_main_fig, handles);
    % Update plots:
    % If plot 5 changed -> recalc sizes. If exist: recreate!
    if p_idx==5
        if isnan( handles.plots.fig_h( p_idx ) )==0
            recreate_spectra = 1;
        else
            recreate_spectra = 0;
        end
        if recreate_spectra == 1
            % close plot
            guiexec_plots_close( [], [], handles.figure1, [], 5)
            % Get handles of main GUI:
            handles = guidata( h_main_fig );
        end
        %reinit sizes:
        [ handles.plots.siz, handles.plots.figpos, handles.plots.axspos ]=...
            gui_plots_init_sizes( handles.plots.siz, 5, handles );
        if recreate_spectra == 1
            %create plot
            [ handles.plots ] = gui_plots_create_figs(handles.figure1,...
                handles, handles.plots, 5, []  );
        end
        % Update handles structure
        guidata( h_main_fig, handles);
        % Fill plot:
        gui_plots_update( h_main_fig, handles, 0 );
    elseif handles.plots.axlims_to_use( p_idx,a_idx )==1
        gui_plots_update( h_main_fig, handles, 0 );%Update using 0
    end
    %Disable applying
    set( handles.plotlimbox.apply,'Enable','Off',...
        'BackgroundColor',[0.941,0.941,0.941]);
    %Set editbox bgcolor:
    set( handles.plotlimbox.al_ed(1),...
        'BackgroundColor',[0.941,0.941,0.941]);
    set( handles.plotlimbox.al_ed(2),...
        'BackgroundColor',[0.941,0.941,0.941]);
    set( handles.plotlimbox.al_ed(3),...
        'BackgroundColor',[0.941,0.941,0.941]);
    set( handles.plotlimbox.al_ed(4),...
        'BackgroundColor',[0.941,0.941,0.941]);
    %Set checkbox bgcolor:
    set( handles.plotlimbox.check,...
        'BackgroundColor',[0.941,0.941,0.941]);
end





% ------------------------------------------------------------
% Plot-limits-box: close
% ------------------------------------------------------------
function guiexec_plotlimbox_closerequest( src, eventdata, h_fig1 )
sure = questdlg(...
    'Settings not saving at close automatically! Close anyway?',...
    'Sure?',...
    'Close','Cancel','Cancel');
%Case Cancel
if strcmp(sure,'Cancel')==1 || strcmp(sure,'');
    return
end
handles = guidata(h_fig1);
%Close window:
delete( handles.plotlimbox.fig )
%delete data:
handles.plotlimbox = [];
% Update handles structure
guidata( h_fig1, handles);
%Toggle off button:
set( handles.ui_plot_lims_set_togglebutton, 'Value', 0)









% ------------------------------------------------------------
% Plot-config-box: callbacks
% ------------------------------------------------------------
function guiexec_plotsetupbox_callbacks( src, eventdata, ...
    h_main_fig, curr_uicontrol)
handles = guidata(h_main_fig);
%--------------------------------------------------------------------------
%Common inputs
if handles.plotsetupbox.ft_pm == curr_uicontrol
    %Set editbox bgcolor:
    set( curr_uicontrol,...
        'BackgroundColor',[1,0.6,0.6]);
    %Enable applying
    set( handles.plotsetupbox.app_co,'Enable','On',...
        'BackgroundColor',[1,0.6,0.6])
end
if handles.plotsetupbox.w_ed == curr_uicontrol
    %Set editbox bgcolor:
    set( curr_uicontrol,...
        'BackgroundColor',[1,0.6,0.6]);
    %Enable applying
    set( handles.plotsetupbox.app_co,'Enable','On',...
        'BackgroundColor',[1,0.6,0.6])
end
if handles.plotsetupbox.fs_ed == curr_uicontrol
    %Set editbox bgcolor:
    set( curr_uicontrol,...
        'BackgroundColor',[1,0.6,0.6]);
    %Enable applying
    set( handles.plotsetupbox.app_co,'Enable','On',...
        'BackgroundColor',[1,0.6,0.6])
end
if handles.plotsetupbox.dw_ed == curr_uicontrol
    %Set editbox bgcolor:
    set( curr_uicontrol,...
        'BackgroundColor',[1,0.6,0.6]);
    %Enable applying
    set( handles.plotsetupbox.app_co,'Enable','On',...
        'BackgroundColor',[1,0.6,0.6])
end
%--------------------------------------------------------------------------
%Common: Store values from editboxes 
if handles.plotsetupbox.app_co == curr_uicontrol
    pm_cell = get( handles.plotsetupbox.ft_pm,'String');
    input_ft_name = pm_cell(get( handles.plotsetupbox.ft_pm,'Value'));
    input_w_mm   = str2double(...
        regexprep( get( handles.plotsetupbox.w_ed,'String'), ',', '.'));
    input_ft_pt  = str2double(...
        regexprep( get( handles.plotsetupbox.fs_ed,'String'), ',', '.'));
    input_w_disp = str2double(...
        regexprep( get( handles.plotsetupbox.dw_ed,'String'), ',', '.'));
    %Check inputs:
    tmp_siz_ft_mm = input_ft_pt*(1/72*25);
    if isnan( input_w_mm )~=0 || isnan( input_ft_pt )~=0 ||...
            isnan( input_w_disp )~=0
        errordlg('At least one input invalid!','Error');
        return
    elseif input_ft_pt < 4
        errordlg('Font size too small!','Error');
        return
    elseif input_w_mm < 8*tmp_siz_ft_mm
        errordlg('Image width too small!','Error');
        return
    elseif input_w_disp <=10
        errordlg('Plot on-screen width too small!','Error');
        return
    elseif input_w_disp > 100
        errordlg('Plot on-screen width too large!','Error');
        return
    end
    %Store inputs:
    handles.plots.ft_name = input_ft_name{1};
    handles.plots.siz.printwidth_mm(:) = input_w_mm;
    handles.plots.siz.printwidth_mm(5) =...
        handles.plots.siz.printwidth_mm(5)/2;
    handles.plots.siz.ft_pt = input_ft_pt;
    handles.plots.siz.used_screen_wr = input_w_disp/100;
    % ReInitialize plot properties
    [ handles.plots.siz, handles.plots.figpos, handles.plots.axspos ] =...
        gui_plots_init_sizes( handles.plots.siz, 0, handles );
    % Update handles structure
    guidata( h_main_fig, handles);
    %Warning: takes effect at reopen
    warndlg('New setting will take effect only after reopening plots',...
        'Warning','modal')
    %Set uicontrols:
    %Disable applying
    set( handles.plotsetupbox.app_co,'Enable','Off',...
        'BackgroundColor',[0.941,0.941,0.941]);
    %Set editbox bgcolor:
    set( handles.plotsetupbox.ft_pm,...
        'BackgroundColor',[0.941,0.941,0.941]);
    set( handles.plotsetupbox.w_ed,...
        'BackgroundColor',[0.941,0.941,0.941]);
    set( handles.plotsetupbox.fs_ed,...
        'BackgroundColor',[0.941,0.941,0.941]);
    set( handles.plotsetupbox.dw_ed,...
        'BackgroundColor',[0.941,0.941,0.941]);
end
%--------------------------------------------------------------------------
%Plot-specific: radiobutton change -> update editboxes:
p_val = find( handles.plots.adjustable_axlabels );
update_needed = 0;
p_sel = find( handles.plotsetupbox.p_rb == curr_uicontrol );
if isempty( p_sel )==0
    %set other rb-s to off
    set( handles.plotsetupbox.p_rb( p_val ), 'Value',0);
    set( handles.plotsetupbox.p_rb( p_sel ), 'Value',1);
    %Change axes selection to 1:
    set( handles.plotsetupbox.a_rb( : ), 'Value',0);
    set( handles.plotsetupbox.a_rb( 1 ), 'Value',1);
    %Set legend checkbox:
    set( handles.plotsetupbox.ag_cb,...
        'Value',handles.plots.axs_legend_on(p_sel));
    if handles.plots.adjustable_axlegend(p_sel)==1
        set( handles.plotsetupbox.ag_cb, 'Enable','On');
    else
        set( handles.plotsetupbox.ag_cb, 'Enable','Off');
    end
    %Set apply button enable:
    set( handles.plotsetupbox.app_lb, 'Enable','Off');
    %Set axes radiobutton selectability:
    for a_i=1:size( handles.plots.axs_h, 2)
        set( handles.plotsetupbox.a_rb( a_i ),...
            'String', handles.plots.axs_names{ p_sel, a_i } );
        if isempty( handles.plots.axlims_usrdef{ p_sel,a_i } )==0
            set( handles.plotsetupbox.a_rb( a_i ), 'Enable','On');
        else
            set( handles.plotsetupbox.a_rb( a_i ), 'Enable','Off');
        end
    end
    update_needed = 1;
end
a_sel = find( handles.plotsetupbox.a_rb == curr_uicontrol );
if isempty( a_sel )==0
    %set other rb-s to off
    set( handles.plotsetupbox.a_rb( : ), 'Value',0);
    set( handles.plotsetupbox.a_rb( a_sel ), 'Value',1);
    %Set apply button enable:
    set( handles.plotsetupbox.app_lb, 'Enable','Off');
    update_needed = 1;
end
%--------------------------------------------------------------------------
%Plot-specific: Final selection:
%do not use logicals!
pb_idx = find(cell2mat( get(handles.plotsetupbox.p_rb( p_val ),'Value')));
p_idx = find(...
    handles.plotsetupbox.p_rb == handles.plotsetupbox.p_rb(p_val(pb_idx)));
ab_idx = find(cell2mat( get(handles.plotsetupbox.a_rb( : ),'Value')));
a_idx = find(...
    handles.plotsetupbox.a_rb == handles.plotsetupbox.a_rb(ab_idx));
if update_needed == 1
    %load limits:
    set( handles.plotsetupbox.al_ed(1),...
        'String', handles.plots.axs_label_x_roots{ p_idx,a_idx },...
        'BackgroundColor',[0.941,0.941,0.941]);
    set( handles.plotsetupbox.al_ed(2),...
        'String', handles.plots.axs_label_y_roots{ p_idx,a_idx },...
        'BackgroundColor',[0.941,0.941,0.941]);
    %Disable applying
    set( handles.plotsetupbox.app_lb,'Enable','Off',...
        'BackgroundColor',[0.941,0.941,0.941]);
end
%--------------------------------------------------------------------------
%Plot-specific: editbox changes:
if isempty(find(handles.plotsetupbox.al_ed == curr_uicontrol,1))==0
    %Set editbox bgcolor:
    set( curr_uicontrol,...
        'BackgroundColor',[1,0.6,0.6]);
    %Enable applying
    set( handles.plotsetupbox.app_lb,'Enable','On',...
        'BackgroundColor',[1,0.6,0.6])
end
%Plot-specific: checkbox changes:
if isempty(find(handles.plotsetupbox.ag_cb == curr_uicontrol,1))==0
    %Set editbox bgcolor:
    set( curr_uicontrol,...
        'BackgroundColor',[1,0.6,0.6]);
    %Enable applying
    set( handles.plotsetupbox.app_lb,'Enable','On',...
        'BackgroundColor',[1,0.6,0.6])
end
%--------------------------------------------------------------------------
%Plot-specific: Store values from editboxes and checkbox setting
if handles.plotsetupbox.app_lb == curr_uicontrol
    %Check inputs: not needed
    %Store inputs:
    handles.plots.axs_label_x_roots{ p_idx,a_idx } =...
        get( handles.plotsetupbox.al_ed(1),'String');
    handles.plots.axs_label_y_roots{ p_idx,a_idx } = ...
        get( handles.plotsetupbox.al_ed(2),'String');
    handles.plots.axs_legend_on(p_idx) =...
        get( handles.plotsetupbox.ag_cb, 'Value');
    %Regenerate axis labels:
    if handles.gui.units_id == 10%== cm, cm/s
        %Plot labels:
        handles.plots.axs_label_x = strcat(...
            handles.plots.axs_label_x_roots,...
            handles.plots.axs_label_x_units_cm);
        handles.plots.axs_label_y = strcat(...
            handles.plots.axs_label_y_roots,...
            handles.plots.axs_label_y_units_cm);
    elseif handles.gui.units_id == 1000%== m, m/s
        %Plot labels:
        handles.plots.axs_label_x = strcat(...
            handles.plots.axs_label_x_roots,...
            handles.plots.axs_label_x_units_m);
        handles.plots.axs_label_y = strcat(...
            handles.plots.axs_label_y_roots,...
            handles.plots.axs_label_y_units_m);
    end
    % ReInitialize plot properties
    [ handles.plots.siz, handles.plots.figpos, handles.plots.axspos ] =...
        gui_plots_init_sizes( handles.plots.siz, 0, handles );
    % Update handles structure
    guidata( h_main_fig, handles);
    %Warning: takes effect at reopen
    warndlg('New setting will take effect only after reopening plots',...
        'Warning','modal')
    %Set uicontrols:
    %Disable applying
    set( handles.plotsetupbox.app_lb,'Enable','Off',...
        'BackgroundColor',[0.941,0.941,0.941]);
    %Set editbox bgcolor:
    set( handles.plotsetupbox.al_ed(1),...
        'BackgroundColor',[0.941,0.941,0.941]);
    set( handles.plotsetupbox.al_ed(2),...
        'BackgroundColor',[0.941,0.941,0.941]);
    %Set checkbox color:
    set( handles.plotsetupbox.ag_cb,...
        'BackgroundColor',[0.941,0.941,0.941]);
end





% ------------------------------------------------------------
% Plot-config-box: close
% ------------------------------------------------------------
function guiexec_plotsetupbox_closerequest( src, eventdata, h_fig1 )
sure = questdlg(...
    'Settings not saving at close automatically! Close anyway?',...
    'Sure?',...
    'Close','Cancel','Cancel');
%Case Cancel
if strcmp(sure,'Cancel')==1 || strcmp(sure,'')
    return
end
handles = guidata(h_fig1);
%Close window:
delete( handles.plotsetupbox.fig )
%delete data:
handles.plotsetupbox = [];
% Update handles structure
guidata( h_fig1, handles);




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 
function chapter_05_dialogs_and_functions_for_exporting
% 
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Exporting paramters dialog: dialog callbacks
% -------------------------------------------------------------------------
function gui_exportdlg_pars_callbacks(src,eventdata,h_fig1, dlg_data, b_cl)
% -------------------------------------------------------------------------
% If Cancel -> close dlg
% -------------------------------------------------------------------------
if b_cl==0
    delete( dlg_data.ui_fig )
    return
end
% -------------------------------------------------------------------------
% Else get user selection and start start exporting...
% -------------------------------------------------------------------------
pargroups2export = {};
if get( dlg_data.ui_dat_props, 'Value')==1
    pargroups2export{end+1} = 'dat_props';
end
if get( dlg_data.ui_coordsys, 'Value')==1
    pargroups2export{end+1} = 'coordsys';
end
if get( dlg_data.ui_errfilt, 'Value')==1
    pargroups2export{end+1} = 'errfilt';
end
delete( dlg_data.ui_fig )
wbar = waitbar(0.01,'Exporting parameters to .par file, please wait...',...
    'WindowStyle','modal');%waitbar( 1/10, wbar)
% Get handles
handles = guidata( h_fig1 );% No handles input -> Get handles of main GUI:
% Get selected file index
f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
%Update Log:
handles.log{end+1,1} = sprintf('>>> Exporting parameters to parfile(s)');
if isvalid( handles.lbox )==1
    set(handles.lbox,'String',handles.log,'Value',length(handles.log));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Go through files, go through pargroups
% > check whether selected pargroup exportable
% > check whether selected pargroups are overwritten in parfile
gui_pargroups = fieldnames( handles.ts_defaults.p );%Supported pargroups
file_err_marker = false(length( f_idx ),1);
file_wrn_marker = false(length( f_idx ),1);
for fik = 1:length( f_idx )
    %Log
    handles.log{end+1,1} =...
        sprintf('  > %s', handles.fname{ f_idx(fik) } );
    if isvalid( handles.lbox )==1
        set(handles.lbox,'String',handles.log,'Value',length(handles.log));
    end
    %Init msg variables
    curr_log = {};%log shown at the end
    curr_cancel_export = NaN;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check pargroup validity in gui (Export only valid pars)
    curr_pargroups = {};
    % > dat_props -> always valid
    if isempty( find(strcmp('dat_props',pargroups2export),1) )~=1
        curr_pargroups{end+1} = 'dat_props';
    end
    % > coordsys:
    if isempty( find(strcmp('coordsys',pargroups2export),1) )~=1
        %if only nans - gives just a warning
        if sum((isnan(handles.ts( f_idx(fik) ).p.coordsys.measpt_coords)+...
                isnan(handles.ts( f_idx(fik) ).p.coordsys.measpt_offset)+...
                isnan(handles.ts( f_idx(fik) ).p.coordsys.file_coords)))==9
            file_wrn_marker(fik) = 1;
            curr_log{end+1,1} =...
                '    * Coordinate(s) contain only NaNs ';
        end
        curr_pargroups{end+1} = 'coordsys';
    end
    % > errfilt
    % TODO: maybe unneded, as if EDIT==off: always valid)
    if isempty( find(strcmp('errfilt',pargroups2export),1) )~=1
        %Validity check:
        [ err_msg ] = fcn_errfiltpars_check_validity(...
            handles.ts(f_idx(fik)).p.errfilt, [], isfield(handles,'dvl') );
        if isempty( err_msg )~=1
            file_err_marker(fik) = 1;
            curr_log{end+1,1} =...
                '    ! Error: invalid errfilt parameters';
        else
            curr_pargroups{end+1} = 'errfilt';
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Load parset from parfile + check overwriting (if contained) & ask
    if exist([handles.fpath{f_idx(fik)},handles.fname{f_idx(fik)},...
            '.par'],'file')==2
        
        % 1. Load all supported pars from parfile
        [ parfile_pars, parfile_pargroups, unknown_pars, log_msg_c ] =...
            fcn_parfile_load(...
            [ handles.fpath{f_idx(fik)},handles.fname{f_idx(fik)},'.par'],...
            gui_pargroups,...
            handles.ts_defaults.p,...
            handles.viper_props.errf_types_supported(:,1));
        % If parfile_pars empty == even parset not existing
        % -> OVERWRITE FILE with curr_pargroups
        % If parfile newer than gui_parset_ver:
        % -> ABORT == No overwrite possible
        if isempty( parfile_pars ) ~=1 && parfile_pars.parset_props.ver >...
                handles.ts_defaults.p.parset_props.ver
            file_err_marker(fik) = 1;
            curr_log{end+1,1} =...
                '    ! Error: not overwriting parfile of newer version!';
            handles.log(end+1:end+length(curr_log),1) = curr_log;
            continue
        end
        % If parfile of older than gui_parset_version or of current version
        % -> parfile_pars and parfile_pargroups can be used
        
        % 2. Overwrite check for parfile of current file
        % if curr_pargroups{i_g} is in parfile -> ask whether overwrite:
        for i_g=1:length(curr_pargroups)
            if isnan(curr_cancel_export) && isempty(find(strcmpi(...
                    curr_pargroups{i_g}, parfile_pargroups ),1))~=1
                sure = questdlg(sprintf('Selected parameters %s %s%s %s',...
                    'already exist in parfile',...
                    handles.fname{f_idx(fik)},'.par?','Overwrite them?'),...
                    'File exists!','Overwrite','Cancel','Cancel');
                if strcmp(sure,'Overwrite')~=1
                    curr_cancel_export = 1;
                else
                    curr_cancel_export = 0;
                end
            end
        end
        %If cancelling overwrite:
        if curr_cancel_export == 1
            file_err_marker(fik) = 1;
            curr_log{1,end+1} =...
                '    ! Error: not overwriting parfile (not permitted)';
            handles.log(end+1:end+length(curr_log),1) = curr_log;
            continue
        end
    else
        parfile_pars = [];
        parfile_pargroups = {};
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create parset to export
    % parset_props !always!
    export_pars.parset_props = handles.ts_defaults.p.parset_props;
    export_pars.parset_props.id = now;
    % others: depending on selection and source (handles/parfile)
    % Go through pargroups existing in gui:
    % Info: Since dat_props checkbox disabled == always exported
    for i_g=2:length(gui_pargroups)
        %if gui_pargroups{i_g} is among selected -> add from handles
        %else (== gui_pargroups{i_g} is not among selected)
        %     if gui_pargroups{i_g} is in parfile -> add from parfile
        %     else don't add
        if isempty(find(strcmpi(...
                gui_pargroups{i_g},curr_pargroups),1))==0
            export_pars.(gui_pargroups{i_g}) = ...
                handles.ts( f_idx(fik) ).p.(gui_pargroups{i_g});
        elseif isempty(find(strcmpi(...
                gui_pargroups{i_g},parfile_pargroups),1))==0
            export_pars.(gui_pargroups{i_g}) = ...
                parfile_pars.(gui_pargroups{i_g});
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Export: (overwrite already confirmed!)
    % ------------------------------------------------------------
    % write to file:
    err_msg_write = fcn_parfile_write_ascii(...
        [ handles.fpath{ f_idx(fik) },handles.fname{ f_idx(fik) },'.par'],...
        export_pars);
    % -----------------------------------------------------------------
    % Write log:
    if isempty( err_msg_write )~=1%error occured at export
        file_err_marker(fik) = 1;
        curr_log{1,end+1} =...
            '    ! Error: parfile export failed';
    end
    handles.log(end+1:end+length(curr_log),1) = curr_log;
    if file_err_marker(fik) == 0
        handles.log{end+1,1} = '    # done';
    end
    if isvalid( handles.lbox )==1
        set(handles.lbox,'String',handles.log,'Value',length(handles.log));
    end
    
    waitbar(fik/length( f_idx ))
    
end
%Update handles structure
guidata(handles.figure1, handles);
% -------------------------------------------------------------------------
% List errors - if any - and messages
if sum( file_err_marker )~=0%isempty(list_err_files)~=1
    list_err_files = handles.fname( f_idx(file_err_marker) );
    warndlg([{'Skipped exporting some parameters of file(s):'};...
        list_err_files],'Skipped file(s)!','modal');
end
% -------------------------------------------------------------------------
% List warning - if any - and messages
if sum( file_wrn_marker )~=0%isempty(list_err_files)~=1
    list_wrn_files = handles.fname( f_idx(file_wrn_marker) );
    warndlg([{'Parameter exporting with warning for file(s):'};...
        list_wrn_files],'File(s) with warning!','modal');
end
close(wbar);




% -------------------------------------------------------------------------
% Exporting statistical results dialog: ui callbacks
% -------------------------------------------------------------------------
function gui_exportdlg_stat_ui_callbacks(src,eventdata, dlg_data, curr_ui)
% -------------------------------------------------------------------------
% Radiobutton switches
% -------------------------------------------------------------------------
ui_new_unit = zeros(3,1);
if curr_ui == dlg_data.ui_unit_m
    ui_new_unit( 1 ) = 1;
elseif curr_ui == dlg_data.ui_unit_cm
    ui_new_unit( 2 ) = 1;
elseif curr_ui == dlg_data.ui_unit_mm
    ui_new_unit( 3 ) = 1;
end
if isempty( find( ui_new_unit, 1 ))==0
    set(dlg_data.ui_unit_m, 'Value', ui_new_unit(1) )
    set(dlg_data.ui_unit_cm, 'Value', ui_new_unit(2) )
    set(dlg_data.ui_unit_mm, 'Value', ui_new_unit(3) )
end
ui_new_fil = zeros(5,1);
if curr_ui == dlg_data.ui_fil_xls
    ui_new_fil( 1 ) = 1;
elseif curr_ui == dlg_data.ui_fil_mat
    ui_new_fil( 2 ) = 1;
end
if isempty( find( ui_new_fil, 1 ))==0
    set(dlg_data.ui_fil_xls, 'Value', ui_new_fil(1) )
    set(dlg_data.ui_fil_mat, 'Value', ui_new_fil(2) )
end

% -------------------------------------------------------------------------
% Exporting statistical results dialog: dialog callbacks
% -------------------------------------------------------------------------
function gui_exportdlg_stat_callbacks(src,eventdata,handles,dlg_data,b_cl,...
    ena_array, valuarray, cds_array, posena_valu, fname_list, fpath_list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case Cancel -> close dlg
if b_cl==0
    delete( dlg_data.ui_fig )
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate valu_names with appropriate units
% Get user unit selection
if get( dlg_data.ui_unit_m, 'Value')==1
    scaling = handles.gui.units_id/1000;%cm->m: 1/100
    %export_unit_str = 'm/s, m';
    valu_names = strcat(...
        handles.field_props.valu_name_roots,...
        repmat({' '},[1,size( handles.field_props.valu_name_roots, 2)]),...
        handles.field_props.valu_units_m);
    cds_names = { 'X [m]', 'Y [m]', 'Z [m]' };
end
if get( dlg_data.ui_unit_cm, 'Value')==1
    scaling = handles.gui.units_id/10;%cm->cm:1
    %export_unit_str = 'cm/s, cm';
    valu_names = strcat(...
        handles.field_props.valu_name_roots,...
        repmat({' '},[1,size( handles.field_props.valu_name_roots, 2)]),...
        handles.field_props.valu_units_cm);
    cds_names = { 'X [cm]', 'Y [cm]', 'Z [cm]' };
end
if get( dlg_data.ui_unit_mm, 'Value')==1
    scaling = handles.gui.units_id;%cm->mm:10
    %export_unit_str = 'mm/s, mm';
    valu_names = strcat(...
        handles.field_props.valu_name_roots,...
        repmat({' '},[1,size( handles.field_props.valu_name_roots, 2)]),...
        handles.field_props.valu_units_mm);
    cds_names = { 'X [mm]', 'Y [mm]', 'Z [mm]' };
end
%scale factors
valu_scaling = scaling.^handles.field_props.valu_unitscalingpower;
cds_scaling = [scaling, scaling, scaling];
%scaled valu-s
valuarray = valuarray.*repmat( valu_scaling, [size( valuarray, 1),1] );
cds_array = cds_array.*repmat( cds_scaling,  [size( cds_array, 1),1] );
% Get user file type selection
if get( dlg_data.ui_fil_xls, 'Value')==1
    export_typ = 'xls-pts';
end
if get( dlg_data.ui_fil_mat, 'Value')==1
    export_typ = 'mat-pts';
end
%Close dialog
delete( dlg_data.ui_fig )
%Update Log:
handles.log{end+1,1} = sprintf('>>> Exporting statistical results...');
if isvalid( handles.lbox )==1
    set(handles.lbox,'String',handles.log,'Value',length(handles.log));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exporting: transform data to needed dataset AND EXPORT
if strcmp( export_typ, 'xls-pts')==1
    %Ask a file name for saving:
    [Filename,Pathname,Ftypeidx]=uiputfile(...
        {'*.xls','Save to Excel file (*.xls)'},...
        'Please select file name for exporting','adv_results');
    if isequal(Filename,0)
        errordlg('No file selected, export cancelled!',' Error! ');
        return
    end
    [path,Fname,ext] = fileparts(Filename);
    if strcmp(ext,'.xls')==0
        errordlg('Wrong file extension, export cancelled!',' Error! ');
        return
    end
    if exist([Pathname,Filename],'file')==2
        sure = questdlg('This overwrites/replaces existing file!',...
            'Overwrite file?','Continue','Cancel','Cancel');
        if strcmp(sure,'Continue')~=1 %strcmp(sure,'Cancel')==1 || strcmp(sure,'');
            return
        end
        
    end
    wbar = waitbar(0.01,'Exporting data to an Excel sheet, please wait...',...
        'WindowStyle','modal');
    %Build dataset:
    export_vars = cell(0);
    %Filenames:
    header = {'Path','Filename'};
    %Sorting based on coordinates: (if a coordinate is nan -> no sorting)
    if isempty(find( isnan( cds_array ), 1))==1%no nans
        % Sort files by coordinates:
        [sorted_cds, sorted ] = sortrows( cds_array );
    else
        sorted = (1:size( cds_array, 1))';
    end
    header = cat(2,header,cds_names);
    filenames = cat(2,fpath_list(sorted),fname_list(sorted));
    export_vars = cat(2,export_vars,num2cell( cds_array( sorted, :) ));
    %exportable (valid) valu-s and valu names:
    %Info: MATLAB converts NaN values to empty Excel cells -> coord nans OK
    header =...
        cat(2,header,valu_names(posena_valu));
    export_vars =...
        cat(2,export_vars,num2cell( valuarray(sorted,posena_valu) ));
    %Delete old excel, if existing:
    if exist([Pathname,Filename],'file')==2
        delete([Pathname,Filename])
    end
    %Write to xls:
    xlswrite([Pathname,Filename], header,...
        'statistics', 'A1')
    xlswrite([Pathname,Filename], filenames,...
        'statistics', 'A2')
    xlswrite([Pathname,Filename], export_vars,...
        'statistics', 'C2')
end
if strcmp( export_typ, 'mat-pts')==1
    %Ask a file name for saving:
    [Filename,Pathname,Ftypeidx]=uiputfile(...
        {'*.mat','MAT file (*.mat)'},...
        'Please select file name for exporting','adv_results');
    if isequal(Filename,0)
        errordlg('No file selected, export cancelled!',' Error! ');
        return
    end
    [path,Fname,ext] = fileparts(Filename);
    if strcmp(ext,'.mat')==0
        errordlg('Wrong file extension, export cancelled!',' Error! ');
        return
    end
    wbar = waitbar(0.01,'Exporting data to a MAT-file, please wait...',...
        'WindowStyle','modal');
    %Build dataset:
    variable_names_in_columns = cell(0);
    %coordinates: Sort files by coordinates. CASE nan present -> no sorting
    if isempty(find( isnan( cds_array ), 1))==1%no nans
        [sorted_cds, sorted ] = sortrows( cds_array );
    else
        sorted = (1:size( cds_array, 1))';
    end
    filenames_in_rows = strcat(fpath_list(sorted),fname_list(sorted));
    %exportable (valid) valu-s and valu names:
    variable_names_in_columns =...
        cat(2,cds_names,valu_names(posena_valu));
    variables =...
        cat(2,cds_array( sorted, :), valuarray(sorted,posena_valu) );
    %Write to mat:
    save([Pathname,Filename],...
        'filenames_in_rows','variable_names_in_columns','variables',...
        '-v7.3');
end

handles.log{end+1,1} =...
	sprintf('  > Finished exporting data to:');
handles.log{end+1,1} =...
	sprintf('    %s%s',...
    Pathname,Filename );
%Update handles structure
guidata(handles.figure1, handles);
%Update Log:
if isvalid( handles.lbox )==1
    set(handles.lbox,'String',handles.log,'Value',length(handles.log));
end
%Closing waitbar:
close(wbar);




% -------------------------------------------------------------------------
% Exporting statistical results dialog: ui callbacks
% -------------------------------------------------------------------------
function gui_exportdlg_plot_ui_callbacks(src,eventdata, dlg_data, curr_ui)
% -------------------------------------------------------------------------
% Radiobutton switches
% -------------------------------------------------------------------------
ui_new_fil = zeros(5,1);
if curr_ui == dlg_data.ui_fil_tif
    ui_new_fil( 1 ) = 1;
elseif curr_ui == dlg_data.ui_fil_fig
    ui_new_fil( 2 ) = 1;
end
if isempty( find( ui_new_fil, 1 ))==0
    set(dlg_data.ui_fil_tif, 'Value', ui_new_fil(1) )
    set(dlg_data.ui_fil_fig, 'Value', ui_new_fil(2) )
end

% -------------------------------------------------------------------------
% Exporting statistical results dialog: dialog callbacks
% -------------------------------------------------------------------------
function gui_exportdlg_plot_callbacks(src,eventdata,handles,dlg_data,b_cl,...
    pos_active_figs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case Cancel -> close dlg
if b_cl==0
    delete( dlg_data.ui_fig )
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get user file type selection
if get( dlg_data.ui_fil_tif, 'Value')==1
    export_typ = 'tif-file';
    fil_ext = '.tif';
end
if get( dlg_data.ui_fil_fig, 'Value')==1
    export_typ = 'fig-file';
    fil_ext = '.fig';
    warning off
end
%Close dialog
delete( dlg_data.ui_fig )
%Update Log:
handles.log{end+1,1} =...
	sprintf('>>> Exporting plot(s) to images file(s)...');
if isvalid( handles.lbox )==1
    set(handles.lbox,'String',handles.log,'Value',length(handles.log));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get selected file index
f_lbx = get(handles.ui_listbox1,'Value');
f_idx = handles.sorted_index(  f_lbx  );
% f_idx = handles.sorted_index(  get(handles.ui_listbox1,'Value')  );
% ------------------------------------------------------------
% Files one-by-one
% ------------------------------------------------------------
for fik=1:length(f_idx)
    %Display log:
    handles.log{end+1,1} =...
    	sprintf('  > %s... ',handles.fname{f_idx(fik)} );
    if isvalid( handles.lbox )==1
        set(handles.lbox,'String',handles.log,'Value',length(handles.log));
    end
    % Select target dir
    saveto = handles.fpath{ f_idx(fik) };
    if strcmp(saveto(end), filesep )~=1
        saveto = [ saveto filesep];
    end
    % Change file selection & update figs
    set(handles.ui_listbox1,'Value',f_lbx(fik));
    gui_plots_update(handles.figure1, handles, 0 );
    % ------------------------------------------------------------
    % Plot types one-by-one
    % ------------------------------------------------------------
    for mik=1:length( pos_active_figs )
        %export_name
        save_filename = [ handles.fname{ f_idx(fik) }, ...
            handles.plots.expname_suffix{pos_active_figs(mik)}, fil_ext];
        %overwrite check:
        if isempty( dir([ saveto, save_filename ]))~=1
            sure=questdlg(sprintf(...
                'Image file %s%s already exist. %s',...
                handles.fpath{ f_idx(fik) }, save_filename,...
                'Overwrite?' ),...
                'Overwrite?','Overwrite','Cancel','Cancel');
            if strcmp(sure,'Overwrite')~=1
                continue
            end
            tmp_log_done = 'overwritten';
        else
            tmp_log_done = 'done';
        end
        % Export
        if strcmp( export_typ, 'tif-file')==1
            print( handles.plots.fig_h( pos_active_figs(mik) ),...
                '-r300', '-dtiff',...
                [saveto, save_filename]);
        end
        if strcmp( export_typ, 'fig-file')==1
            savefig( handles.plots.fig_h( pos_active_figs(mik) ),...
                [saveto, save_filename],'compact')
        end
        % ------------------------------------------------------------
        %Display log:
        handles.log{end+1,1} =...
            sprintf('    * %s: %s ',...
            handles.plots.expname_suffix{ pos_active_figs(mik) },...
            tmp_log_done);
        if isvalid( handles.lbox )==1
            set(handles.lbox,'String',handles.log,'Value',length(handles.log));
        end
    end
end
%Update handles structure
guidata(handles.figure1, handles);
%Update GUI:
gui_update_ui_values(handles.figure1, handles)
gui_update_ui_enables(handles.figure1, handles)
gui_plots_update(handles.figure1, handles, 0 );
%Display log:
handles.log{end+1,1} =...
	sprintf('### Finished');
if isvalid( handles.lbox )==1
    set(handles.lbox,'String',handles.log,'Value',length(handles.log));
end
warning on
%Closing waitbar:
% close(wbar);





% -------------------------------------------------------------------------
% Exporting statistical results dialog: ui callbacks
% -------------------------------------------------------------------------
function gui_exportdlg_ts_ui_callbacks(src,eventdata, dlg_data, curr_ui)
% -------------------------------------------------------------------------
% Radiobutton switches
% -------------------------------------------------------------------------
ui_new_repl = zeros(3,1);
if curr_ui == dlg_data.ui_repl_999
    ui_new_repl( 1 ) = 1;
elseif curr_ui == dlg_data.ui_repl_nan
    ui_new_repl( 2 ) = 1;
elseif curr_ui == dlg_data.ui_repl_int
    ui_new_repl( 3 ) = 1;
end
if isempty( find( ui_new_repl, 1 ))==0
    set(dlg_data.ui_repl_999, 'Value', ui_new_repl(1) )
    set(dlg_data.ui_repl_nan, 'Value', ui_new_repl(2) )
    set(dlg_data.ui_repl_int, 'Value', ui_new_repl(3) )
end
ui_new_fil = zeros(5,1);
if curr_ui == dlg_data.ui_fil_asc
    ui_new_fil( 1 ) = 1;
elseif curr_ui == dlg_data.ui_fil_mat
    ui_new_fil( 2 ) = 1;
end
if isempty( find( ui_new_fil, 1 ))==0
    set(dlg_data.ui_fil_asc, 'Value', ui_new_fil(1) )
    set(dlg_data.ui_fil_mat, 'Value', ui_new_fil(2) )
end
if get(dlg_data.ui_fil_asc, 'Value') == 1
    set(dlg_data.ui_repl_999, 'Enable', 'on' )
    set(dlg_data.ui_repl_nan, 'Enable', 'on' )
    set(dlg_data.ui_repl_int, 'Enable', 'on' )
else
    set(dlg_data.ui_repl_999, 'Enable', 'off' )
    set(dlg_data.ui_repl_nan, 'Enable', 'off' )
    set(dlg_data.ui_repl_int, 'Enable', 'off' )
end

% -------------------------------------------------------------------------
% Exporting statistical results dialog: dialog callbacks
% -------------------------------------------------------------------------
function gui_exportdlg_ts_callbacks(src,eventdata,handles,dlg_data,b_cl,...
    f_idx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case Cancel -> close dlg
if b_cl==0
    delete( dlg_data.ui_fig )
    return
end
%Update Log:
handles.log{end+1,1} = sprintf('>>> Exporting velocity time series...');
if isvalid( handles.lbox )==1
    set(handles.lbox,'String',handles.log,'Value',length(handles.log));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get user selections
if get( dlg_data.ui_fil_asc, 'Value')==1
    export_typ = 'asc-ts';
end
if get( dlg_data.ui_fil_mat, 'Value')==1
    export_typ = 'mat-ts';
end
if get( dlg_data.ui_repl_999, 'Value')==1
    export_repl = 999;
end
if get( dlg_data.ui_repl_nan, 'Value')==1
    export_repl = nan;
end
if get( dlg_data.ui_repl_int, 'Value')==1
    export_repl = [];
end
%Close dialog
delete( dlg_data.ui_fig )
%Update Log:
handles.log{end+1,1} = sprintf('>>> Exporting statistical results...');
if isvalid( handles.lbox )==1
    set(handles.lbox,'String',handles.log,'Value',length(handles.log));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
export_unit_str = 'm, m/s';
scaling = handles.gui.units_id/1000;%cm->m: 1/100
% Exporting: collect data AND EXPORT
if strcmp( export_typ, 'mat-ts')==1
    %-----------------------------------------------------
    %Generate dataset to export:
    %Export units:
    temp_data.units = export_unit_str;
    %Descripitom:
    temp_data.description = {...
        'vel_nans: columns [u,v,w], bad samples replaced by NaN';...
        'vel_intp: columns [u,v,w], bad samples replaced by linear interpolation'};
    %Collect data (all files are enabled - checked in previous gui)
    for fik=1:length(f_idx)
        %Filename
        temp_data.ts( fik ).fullname = handles.fullname{ f_idx(fik),1};
        temp_data.ts( fik ).fname = handles.fname{ f_idx(fik), 1};
        temp_data.ts( fik ).fpath = handles.fpath{ f_idx(fik), 1};
        %time
        temp_data.ts( fik ).t = handles.ts( f_idx(fik) ).seri_x.t;
        %Interpolated time-series
        %raw@valid pos, interpolated@erroneous samples
        temp_data.ts( fik ).vel_intp =...
            scaling*(...
            handles.ts( f_idx(fik) ).seri_data{1} +...
            repmat(( handles.ts( f_idx(fik) ).valu_data(1:3) +...
            handles.ts( f_idx(fik) ).intl.mean_correctn ),...
            [ size( handles.ts( f_idx(fik) ).seri_data{1},1),1 ]));
        %Not interpolated time-series
        %raw@valid pos, nan@erroneous samples
        temp_raw = scaling*handles.ts( f_idx(fik) ).raw_veldata;
        temp_raw( handles.ts( f_idx(fik) ).vals==0, : ) = nan;
        temp_data.ts( fik ).vel_nans =...
            temp_raw( handles.ts( f_idx(fik) ).intl.used_ts_seg, 1:3 );
    end
    %-----------------------------------------------------
    %Ask a file name for saving:
    [Filename,Pathname,Ftypeidx]=uiputfile(...
        {'*.mat','Matlab file (*.mat)'},...
        'Please select file name for exporting','viper_timeseries');
    if isequal(Filename,0)
        errordlg('No file selected, export cancelled!',' Error! ');
        return
    end
    [path,Fname,ext] = fileparts(Filename);
    if strcmp(ext,'.mat')==0
        errordlg('Wrong file extension, export cancelled!',' Error! ');
        return
    end
    wbar = waitbar(0.01,'Exporting data to a MAT-file, please wait...',...
        'WindowStyle','modal');
    %-----------------------------------------------------
    %Writing data:
    save([Pathname,Filename], '-struct', 'temp_data', '-v7.3' );
end
if strcmp( export_typ, 'asc-ts')==1
    %-----------------------------------------------------
    %Check for overwrite:
    f_idx_notoverwriting = [];
    list_overwriting  = {};
    for fik=1:length(f_idx)
        if exist([ handles.fpath{f_idx(fik)},...
                handles.fname{f_idx(fik)}, '.asc'],'file')==2
            list_overwriting{1,end+1} = handles.fullname{ f_idx(fik)};
        else
            f_idx_notoverwriting(end+1,1) = f_idx(fik);
        end
    end
    if isempty( list_overwriting )==0
        sure = questdlg(['Overwrite asc-files of following files:',...
            list_overwriting],...
            'Overwrite?',...
            'Overwrite files','Skip files','Cancel','Skip files');
        if strcmp(sure,'Overwrite files')==1
            f2exp = f_idx;
        elseif strcmp(sure,'Skip files')==1
            f2exp = f_idx_notoverwriting;
        else
            handles.log{end+1,1} =...
                sprintf('  # Exporting cancelled');
            %Update Log:
            if isvalid( handles.lbox )==1
                set(handles.lbox,'String',handles.log,'Value',length(handles.log));
            end
            errordlg('Export cancelled!',' Error! ');
            return
        end
    else
        f2exp = f_idx_notoverwriting;
    end
    %-----------------------------------------------------
    %Export
    wbar = waitbar(0.01,'Exporting data to ascii files, please wait...',...
        'WindowStyle','modal');
    for fik=1:length(f2exp)
        %Generate data:
        if isempty( export_repl )==1%interolated
            %Interpolated time-series
            %raw@valid pos, interpolated@erroneous samples
            vel2export =...
                scaling*(...
                handles.ts( f2exp(fik) ).seri_data{1} +...
                repmat(( handles.ts( f2exp(fik) ).valu_data(1:3) +...
                handles.ts( f2exp(fik) ).intl.mean_correctn ),...
                [ size( handles.ts( f2exp(fik) ).seri_data{1},1),1 ]));
        else
            %Not interpolated time-series
            %raw@valid pos, nan@erroneous samples
            temp_raw = scaling*handles.ts( f_idx(fik) ).raw_veldata;
            temp_raw( handles.ts( f_idx(fik) ).vals==0, : ) = export_repl;
            vel2export =...
                temp_raw( handles.ts( f_idx(fik) ).intl.used_ts_seg, 1:3 );
        end
        data2export = [handles.ts( f2exp(fik) ).seri_x.t, vel2export] ;
        %Write data to file
        dlmwrite(...
            [handles.fpath{f2exp(fik)},handles.fname{f2exp(fik)},'.asc'],...
            data2export,...
             '\t')
        waitbar( fik/length(f2exp), wbar)
%         %Collect
%         %time
%         temp_data.ts( fik ).t = handles.ts( f_idx(fik) ).seri_x.t;
%         %Time-series - raw@valid pos, interpolated@erroneous samples
%         temp_data.ts( fik ).uvw_intp =...
%             scaling*(...
%             handles.ts( f_idx(fik) ).seri_data{1} +...
%             repmat(( handles.ts( f_idx(fik) ).valu_data(1:3) +...
%             handles.ts( f_idx(fik) ).intl.mean_correctn ),...
%             [ size( handles.ts( f_idx(fik) ).seri_data{1},1),1 ]));
%         %Time-series - raw@valid pos, nan@erroneous samples
%         %Find segment of raw in seri
%         used4seri = false(size( handles.ts(f_idx(fik)).vals ));
%         used4seri( find( handles.ts(f_idx(fik)).vals==1,1,'first') : ...
%             find( handles.ts(f_idx(fik)).vals==1,1,'last' )) = true;
%         %Read raw for used segment:
%         temp_data.ts( fik ).uvw_nans =...
%             handles.ts( f_idx(fik) ).raw_veldata( used4seri, 1:3 );
%         %Replace invalid/erroneus samples by nan:
%         temp_data.ts( fik ).uvw_nans(...
%             handles.ts( f_idx(fik) ).vals( used4seri )==0, 1:3 )= nan;        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finish
handles.log{end+1,1} =...
	sprintf('  # Finished exporting data');
%Update handles structure
guidata(handles.figure1, handles);
%Update Log:
if isvalid( handles.lbox )==1
    set(handles.lbox,'String',handles.log,'Value',length(handles.log));
end
%Closing waitbar:
close(wbar);















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 
function chapter_11_file_read_and_write_functions
% 
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ------------------------------------------------------------------------
% vno-dat-file: findout vectrino parameters 
% ------------------------------------------------------------------------
function [ temp_dat_props, input_v_unit, temp_log ] =...
    fcn_addfile_findout_dat_props(...
    temp_type, temp_data, temp_hdr, curr_filename, curr_t_col, curr_u_col)
%temp_type: 'vno' 'vec' 'usr'
input_v_unit = 1;
% ------------------------------------------------------------
%%% 1. type:
% ------------------------------------------------------------
temp_dat_props.type = temp_type;

% ------------------------------------------------------------
%%% 2. t-columns, v-columns AND vel_units!
% ------------------------------------------------------------
u_col = [];%initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A. Try automatically for temp_types with known pattern
if strcmp(temp_type, 'vno')==1
    %input_v_unit = 1;
    %.dat-columns are: 1&2= opt. 3=ensemble 4=status
    % try to find out status-column, where values seem to be constant
    % within a .dat-file => find col with constant value
    status_col = find( [...
        numel( unique( temp_data(:,1) )),...
        numel( unique( temp_data(:,2) )),...
        numel( unique( temp_data(:,3) )),...
        numel( unique( temp_data(:,4) )) ] ==1);
    if numel( status_col )==1%only one col with constant value
        u_col = status_col+1;
        %time-col = status_col-2
        if status_col>1
            t_col = status_col-2;
        else
            t_col = 0;
        end
        %text for msg-box: source of parameter
        temp_msg_ext = 'basis: status column';
    end
    %if no column with constant value -> leave u_col empty
end
if strcmp(temp_type, 'vec')==1
    %input_v_unit = 1;
    %.dat-columns are fixed: 1=burst counter, 2=ensemble, 3:5=velocities
    t_col = 0;
    u_col = 3;
    %text for msg-box: source of parameter
    temp_msg_ext = 'basis: dat-format';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%B. Ask user - then check!
if isempty( u_col )==1 %|| isempty( input_v_unit )==1;
    % Display data in a table
    figtab = figure('Position', [100 400 752 200]);
    tab_sel =...
        uitable('Parent', figtab, 'Position', [25 25 700 150]);
    set(tab_sel, 'Data', temp_data(1:5,:));
    % Ask user for column numbers
    col_input_str{1} = '';
    col_input_str{2} = '';
    repeat_input = 0;
    while isempty(col_input_str)==1 || isempty(col_input_str{1})==1 ||...
            isempty(col_input_str{2})==1 || repeat_input == 0
        col_input_str = inputdlg(...
            {'Column containing time (0 means: no such variable):';...
            'First column containing velocities:'},...
            'File format definition',1,...
            {num2str( curr_t_col ) ;...
            num2str( curr_u_col )});
        repeat_input = 1;
        col_input = str2double( regexprep(col_input_str, ',', '.') );
        %TODO: Cancel loading option! -> err_msg output!
        % Check user inputs:
        %If numbers are non-integer || too large || too small:
        if rem( col_input(1),1)~=0 || rem( col_input(2),1)~=0 ||...
                col_input(1)<0 || col_input(2)<0 ||...
                isnan( col_input(1) )==1 || isnan( col_input(2) )==1
            %invalid number
            dlg_h = errordlg(...
                'Invalid number','Input Error','modal');
            waitfor(dlg_h);
            repeat_input = 0;
        elseif col_input(1)>size( temp_data ,2) ||...
                col_input(2)>size( temp_data ,2)
            %such column not existing
            dlg_h = errordlg(...
                'Invalid column number','Input Error','modal');
            waitfor(dlg_h);
            repeat_input = 0;
        end
        %Dat type specific:
        if strcmp( temp_type, 'vno' )==1
            %Ensemble counter column = u_col-2
            if u_col<3
                dlg_h = errordlg(sprintf(...
                    'Number too small.\n%s%s',...
                    'For Vectrino data, the number of the ',...
                    'first velocity column has to be larger than 2!'),...
                    'Input Error','modal');
                waitfor(dlg_h);
                repeat_input = 0;
            end
            %Number of columns required
            if u_col > size( temp_data, 2 )-15
                dlg_h = errordlg(sprintf(...
                    'Number too large.\n%s%s',...
                    'For Vectrino data, the first velocity column ',...
                    'has to be followed by at least 15 columns!'),...
                    'Input Error','modal');
                waitfor(dlg_h);
                repeat_input = 0;
            end
        end
        if strcmp( temp_type, 'vno' )==1
            %Ensemble counter column = u_col-2
            if u_col<3
                dlg_h = errordlg(sprintf(...
                    'Number too small.\n%s%s',...
                    'For Vector data, the number of the ',...
                    'first velocity column has to be larger than 2!'),...
                    'Input Error','modal');
                waitfor(dlg_h);
                repeat_input = 0;
            end
            %Number of columns required
            if u_col > size( temp_data, 2 )-11
                dlg_h = errordlg(sprintf(...
                    'Number too large.\n%s%s',...
                    'For Vector data, the first velocity column ',...
                    'has to be followed by at least 11 columns!'),...
                    'Input Error','modal');
                waitfor(dlg_h);
                repeat_input = 0;
            end
        end
        if strcmp( temp_type, 'usr' )==1
            if col_input(1)<1
                dlg_h = errordlg(sprintf(...
                    'For User-defined data, time column is required!'),...
                    'Input Error','modal');
                waitfor(dlg_h);
                repeat_input = 0;
            end
            %Number of columns required
            if u_col > size( temp_data, 2 )-2
                dlg_h = errordlg(sprintf(...
                    'Number too large.\n%s%s',...
                    'For User-defined data, the first velocity column ',...
                    'has to be followed by at least 2 columns!'),...
                    'Input Error','modal');
                waitfor(dlg_h);
                repeat_input = 0;
            end
        end
    end
    t_col = col_input(1);
    u_col = col_input(2);
    %text for msg-box: source of parameter
    temp_msg_ext = 'basis: user input';
else
    figtab = [];
end
close(figtab);
% Store answers
temp_dat_props.column_t = t_col;
temp_dat_props.column_u = u_col;
%Create messages:
temp_log = sprintf('      * columns t:%d, u:%d (%s);', t_col, u_col, temp_msg_ext);
% ------------------------------------------------------------
%%% 3. freq
frq = [];%initialize
%Try automatically: load from hdr for known temp_types
if strcmp(temp_type, 'vno')==1 || strcmp(temp_type, 'vec')==1
    if temp_hdr.existing == 1
        frq = str2double( temp_hdr.sampling_rate{1} );
        %text for msg-box: source of parameter
        temp_msg_ext = 'from .hdr-file';
    end
end
%Ask user:
if isempty( frq )==1
    if t_col ~= 0 % if time_col exists:
        %Frequency = How many data during one second?
        frq = find( round(temp_data(:,t_col)*1000) ==...
            round((temp_data(1,t_col)+1)*1000) )-1;
        sure=questdlg({...
            'Please confirm the calculated frequency for the file: ';...
            curr_filename; ' ';...
            sprintf( 'The calulated frequency is: %d Hz', frq)},...
            'Correct frequency?','Yes','No','Yes');
        if strcmp(sure,'Yes')==1
            %text for msg-box: source of parameter
            temp_msg_ext = 'user input';
        end
    end
    if t_col==0 || strcmp(sure,'Yes')~=1
        frq_input_str{1} = '1.5';
        frq_input = str2double( regexprep(frq_input_str, ',', '.') );
        while isempty(frq_input_str) == 1 ||...
                isnan( frq_input )== 1 ||...
                rem( frq_input,1) ~= 0 ||...
                frq_input <= 0
            frq_input_str = inputdlg(...
                {'Please give the sampling rate of the measurement:'},...
                curr_filename,1,...
                {'200'});
            frq_input = str2double( regexprep(frq_input_str, ',', '.') );
        end
        frq = frq_input;
        %text for msg-box: source of parameter
        temp_msg_ext = 'user input';
    end
end
temp_dat_props.frq  = frq;
%Create messages:
temp_log = sprintf('%s fs: %d Hz (%s)', temp_log, frq, temp_msg_ext);
% ------------------------------------------------------------
%%% 4.a vel-range
temp_velrng = [];
%Try automatically: load from hdr for known temp_types
if strcmp(temp_type, 'vno')==1 || strcmp(temp_type, 'vec')==1
    if temp_hdr.existing == 1
        % Vel-range
        temp_velrng = temp_hdr.velocity_rnge;
        %text for msg-box: source of parameter
        temp_msg_ext = 'from .hdr-file';
    end
end
%Otherwise:
if isempty( temp_velrng )==1
    % Vel-range
    temp_velrng = {'0','-'};
    %text for msg-box: source of parameter
    temp_msg_ext = 'n.a.';
end
temp_dat_props.probe_nvr = temp_velrng;
%Create messages: none
% ------------------------------------------------------------
%%% 4.b probe_cfg
temp_probe_cfg = [];
%Try automatically for known temp_types
if strcmp(temp_type, 'vno')==1
    if isequal( temp_data(:,u_col+3) , zeros(size( temp_data, 1),1) )==1
        temp_probe_cfg = 'side-looking';
    else
        temp_probe_cfg = 'down-looking';
    end
    %text for msg-box: source of parameter
    temp_msg_ext = 'basis: components';
end
%Otherwise:
if isempty( temp_probe_cfg )==1
    temp_probe_cfg = 'n.a.';
end
temp_dat_props.probe_cfg = temp_probe_cfg;
%Create messages: none
        




% ------------------------------------------------------------------------
% hdr-file: read
% ------------------------------------------------------------------------
function [ hdr_nfo ] = fcn_addfile_hdr_read_vno_vec(...
    hdr_fullfilename, default_hdr_nfo )
%Use only fields defined in handles.ts_defaults.hdr_nfo
hdr_nfo = default_hdr_nfo;
hdr_nfo(1).existing = true;
%%%% Open file for read - read line by line:
fid = fopen( hdr_fullfilename , 'rt');%read text
last_needed_line = 0;
% Get the first line:
tline = fgetl(fid);%tline = fgetl(fid);
while ischar(tline) %feof(fid) == 0
    %line10
    if strncmp('Sampling rate', tline, 13) ==1
        temp_str = textscan( tline, '%*s %*s %s %s');
        hdr_nfo(1).sampling_rate{1} = temp_str{1}{1};
        hdr_nfo(1).sampling_rate{2} = temp_str{2}{1};
    end
    %line11
    if strncmp('Nominal velocity range', tline, 22) ==1
        temp_str = textscan( tline,'%*s %*s %*s %s %s');
        hdr_nfo(1).velocity_rnge{1} = temp_str{1}{1};
        hdr_nfo(1).velocity_rnge{2} = temp_str{2}{1};
    end
    %line12
    if strncmp('Probe type', tline, 10)  ==1
        temp_str = textscan( tline,'%*s %*s %s');
        hdr_nfo(1).probe_cfg = temp_str{1}{1};
    end
    %line13
    if strncmp('Transmit length', tline, 15)  ==1
        temp_str = textscan( tline, '%*s %*s %s %s');
        hdr_nfo(1).transmit_lgth{1} = temp_str{1}{1};
        hdr_nfo(1).transmit_lgth{2} = temp_str{2}{1};
    end
    %line14
    if strncmp('Sampling volume', tline, 15)  ==1
        temp_str = textscan( tline,'%*s %*s %s %s');
        hdr_nfo(1).sampling_volu{1} = temp_str{1}{1};
        hdr_nfo(1).sampling_volu{2} = temp_str{2}{1};
    end
    %line19
    if strncmp('Powerlevel', tline, 10)  ==1
        temp_str = textscan( tline,'%*s %s');
        hdr_nfo(1).power_level = temp_str{1}{1};
    end
    %line19
    if strncmp('Coordinate system', tline, 17)  ==1
        temp_str = textscan( tline,'%*s %*s %s');
        hdr_nfo(1).coord_system = temp_str{1}{1};
    end
    %line71: first one! (second serial no is for head!)
    if strncmp('Serial number', tline, 13) ==1 &&...
            isempty( hdr_nfo(1).serial_no )==1
            %&& isempty( strfind(tline, 'VNO' ) )==0
        temp_str = textscan( tline,'%*s %*s %s %s');
        hdr_nfo(1).serial_no=strcat(temp_str{1}{1},' ',temp_str{2}{1});
        %last_needed_line = 1;
    end
    %linex (only vector)
    if isempty(strfind(tline, 'Ensemble counter'))~=1
        temp_str = textscan( tline,'%s %*[^\n]');
        hdr_nfo(1).ensemble_column = temp_str{1}{1};
        last_needed_line = 1;
    end
    if last_needed_line ==0
        % Get the next line:
        tline = fgetl(fid);%tline = fgetl(fid);
    else
        % Break while loop
        tline = -1;%like if it was eof (in case fgetl gets eof -> tline=-1)
    end
end
status = fclose(fid);%status=... show 0 if succesful
% if status==0
%     fprintf(' ... ASCII-file reading successfull: %s \n', hdrfilename);
% end
clear fid tline status





% ------------------------------------------------------------------------
% parfile: Load parfile content: 
% ------------------------------------------------------------------------
function [ loaded_pars, loaded_pargroups, unknown_pars, log_msg_c ] =...
    fcn_parfile_load( par_fullfilename, needed_pargroups, ...
    gui_default_pars, gui_supported_errf_types)
% Loads pars of needed_pargroups from parfile:
% > into loaded_pars, which has the same records as in handles
% > into unknown_pars those records, which are unknown for current gui
%   * if errfilt contains only unsupported types -> no filt
% gui_default_pars = handles.ts_defaults.p
% (use gui_default_pars.parset_props.ver as gui_parset_ver)
% * parset_props always loaded,'parset_props' not added to loaded_pargroups
% IF Erroneous parfiles => empty(loaded_pars,loaded_pargroups) + err_msg
% IF none of needed par there => empty(loaded_pars,loaded_pargroups)
% IF some needed par there    => filled(loaded_pars,loaded_pargroups)
% IF parfile ver newer than code => loaded_pars, loaded_pargroups, err_msg
% ------------------------------------------------------------------------
%Init:
loaded_pars = [];%parameter set
loaded_pargroups = {};%list of paramtergroups within output
unknown_pars = [];%parameter set containing unknown parameters
log_msg_c = {};%messages - one message per cell
%Check whether file exists:
if exist( par_fullfilename ,'file')~=2
    log_msg_c{end+1} = 'Parfile doesn''t exist';
    return
end
% ------------------------------------------------------------------------
% Read all pars from ascii file to cellstr: (commented parts eliminated)
[ parfile_allpars, err_msg_1 ] =...
    fcn_parfile_read_all_ascii(par_fullfilename);
if isempty( err_msg_1 )~=1%error occured at read
    log_msg_c{end+1} = err_msg_1{1};
    return
end
% ------------------------------------------------------------------------
% Always load parset_props
if isfield( parfile_allpars ,'parset_props')~=1
    log_msg_c{end+1} = 'Parfile does not contain parset_props';
    return
else
    loaded_pars.parset_props = parfile_allpars.parset_props;
end
% Check version and
if parfile_allpars.parset_props.ver < 1.1
    log_msg_c{end+1} = 'Parfile too old and therefore not supported.';
    return
end
if parfile_allpars.parset_props.ver > gui_default_pars.parset_props.ver
    %FIRST LINE (IF ADDED) Message string used in other parts!!!
    log_msg_c{end+1} = 'Warning: parfile version is newer than this code';
end
% ------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. Convert old versions records to current parfile versions records:
if parfile_allpars.parset_props.ver == 1.1
    % Version 1.1:
    % handles.ts_defaults.p.parset_props.ver%num
    % handles.ts_defaults.p.parset_props.id%num
    % handles.ts_defaults.p.dat_props.type%'vno' | 'vec' | 'usr'
    % handles.ts_defaults.p.dat_props.t_col%num
    % handles.ts_defaults.p.dat_props.v_col%num
    % handles.ts_defaults.p.dat_props.fs%num
    % handles.ts_defaults.p.dat_props.probe_cfg%'side-looking' 'down-looking';
    % handles.ts_defaults.p.dat_props.probe_nvr%{};
    % handles.ts_defaults.p.coordsys.unit%=='m';
    % handles.ts_defaults.p.coordsys.file_coords%num(1,3)
    % handles.ts_defaults.p.coordsys.measpt_offset%num(1,3)
    % handles.ts_defaults.p.coordsys.measpt_coords%num(1,3)
    % handles.ts_defaults.p.coordsys.dir_probe_u%'+x' '-x' '+y' '-y'
    % handles.ts_defaults.p.errfilt.used_velcomp%num(1,3) 1|0
    % handles.ts_defaults.p.errfilt.used_types%cell, {} or e.g. {'pst'}
    % handles.ts_defaults.p.errfilt.typ(1).name%string;
    % handles.ts_defaults.p.errfilt.typ(1).usedcc%num(1,3) 1|0
    % handles.ts_defaults.p.errfilt.typ(1).lambda%num
    % handles.ts_defaults.p.errfilt.typ(1).use_hpf%0|1
    % handles.ts_defaults.p.errfilt.typ(1).hpf_freq%num
    %dat_props:
    parfile_allpars.dat_props.column_t =...
        parfile_allpars.dat_props.t_col;
    parfile_allpars.dat_props.column_u =...
        parfile_allpars.dat_props.v_col;
    parfile_allpars.dat_props.frq =...
        parfile_allpars.dat_props.fs;
    parfile_allpars.dat_props = rmfield( parfile_allpars.dat_props,...
        {'t_col','v_col','fs'});
    %errfilt:
    if isfield( parfile_allpars, 'errfilt' ) && ...
            isfield( parfile_allpars.errfilt, 'typ' )
        if isfield( parfile_allpars.errfilt.typ, 'usedcc' )
            for i_t=1:length(parfile_allpars.errfilt.typ)
                parfile_allpars.errfilt.typ(i_t).used_comp =...
                    parfile_allpars.errfilt.typ(i_t).usedcc;
            end
            parfile_allpars.errfilt.typ = rmfield(...
                parfile_allpars.errfilt.typ,...
                {'usedcc'});
        end
        if isfield( parfile_allpars.errfilt.typ, 'use_hpf' )
        for i_t=1:length(parfile_allpars.errfilt.typ)
            %case highpass not to use -> change hipass_frq->[]
            if isempty(parfile_allpars.errfilt.typ(i_t).use_hpf)==1 || ...
                    parfile_allpars.errfilt.typ(i_t).use_hpf == 0
                parfile_allpars.errfilt.typ(i_t).hipass_frq = [];
            else
                parfile_allpars.errfilt.typ(i_t).hipass_frq =...
                    parfile_allpars.errfilt.typ(i_t).hpf_freq;
            end
        end
        parfile_allpars.errfilt.typ = rmfield(...
            parfile_allpars.errfilt.typ,...
            {'use_hpf','hpf_freq'});
        end
    end
    %Message:
    log_msg_c{end+1} = 'converted parameter set from parfile ver. 1.1';
elseif parfile_allpars.parset_props.ver == 1.2
    % Version 1.2:
    % handles.ts_defaults.p.parset_props.ver%num
    % handles.ts_defaults.p.parset_props.id%num
    % handles.ts_defaults.p.dat_props.type%'vno' | 'vec' | 'usr'
    % handles.ts_defaults.p.dat_props.column_t%num
    % handles.ts_defaults.p.dat_props.column_u%num
    % handles.ts_defaults.p.dat_props.frq%num
    % handles.ts_defaults.p.dat_props.probe_cfg%'side-looking' 'down-looking';
    % handles.ts_defaults.p.dat_props.probe_nvr%{}
    % handles.ts_defaults.p.coordsys.unit%=='m';
    % handles.ts_defaults.p.coordsys.file_coords%num(1,3)
    % handles.ts_defaults.p.coordsys.measpt_offset%num(1,3)
    % handles.ts_defaults.p.coordsys.measpt_coords%num(1,3)
    % handles.ts_defaults.p.coordsys.dir_probe_u%'+x' '-x' '+y' '-y'
    % handles.ts_defaults.p.errfilt.used_velcomp%num(1,3) 1|0
    % handles.ts_defaults.p.errfilt.used_types%cell, {} or e.g. {'pst'}
    % handles.ts_defaults.p.errfilt.typ(1).name%string;
    % handles.ts_defaults.p.errfilt.typ(1).used_comp%num(1,3) 1|0
    % handles.ts_defaults.p.errfilt.typ(1).lambda%num
    % handles.ts_defaults.p.errfilt.typ(1).use_hipass%0|1
    % handles.ts_defaults.p.errfilt.typ(1).hipass_frq%num
    %errfilt:
    if isfield( parfile_allpars, 'errfilt' ) && ...
            isfield( parfile_allpars.errfilt, 'typ' )
        if isfield( parfile_allpars.errfilt.typ, 'use_hipass' )
        for i_t=1:length(parfile_allpars.errfilt.typ)
            %case highpass not to use -> change hipass_frq->[]
            if isempty(parfile_allpars.errfilt.typ(i_t).use_hipass)==1 || ...
                    parfile_allpars.errfilt.typ(i_t).use_hipass == 0
                parfile_allpars.errfilt.typ(i_t).hipass_frq = [];
            end
        end
        parfile_allpars.errfilt.typ = rmfield(...
            parfile_allpars.errfilt.typ,...
            {'use_hipass'});
        end
    end
    %Message:
    log_msg_c{end+1} = 'converted parameter set from parfile ver. 1.2';
end
% ------------------------------------------------------------------------
% Remove unsupported erffilt types - go through types used in par-file
if isfield( parfile_allpars, 'errfilt' )==1
for i_e = 1:length( parfile_allpars.errfilt.used_types )
    %if an err.filt_type unknown - move to unknown and remove it
    if isempty(find(strcmpi(...
            parfile_allpars.errfilt.used_types{i_e},...
            gui_supported_errf_types),1))==1
        %check if unknown_pars has the field -> if not, create
        if isfield( unknown_pars, 'errfilt' )~=1
            unknown_pars.errfilt.used_types = {};
            temp_fields = fieldnames( parfile_allpars.errfilt.typ(i_e) );
            for i_f=1:length( temp_fields )
                unknown_pars.errfilt.typ.(temp_fields{i_f}) = [];
            end
            unknown_pars.errfilt.typ(1)=[];
        end
        %move them to unkown_pars
        unknown_pars.errfilt.used_types{end+1} =...
            parfile_allpars.errfilt.used_types{i_e};
        unknown_pars.errfilt.typ(end+1) = parfile_allpars.errfilt.typ(i_e);
        %Delete from allpars:
        parfile_allpars.errfilt.used_types(i_e) = [];
        parfile_allpars.errfilt.typ(i_e) = [];
    end
    %if there was not other err.filt type .> remove fields of .typ
    if isempty(parfile_allpars.errfilt.used_types)==1
        parfile_allpars.errfilt = rmfield( parfile_allpars.errfilt,'typ' );
    end
end
end
% ------------------------------------------------------------------------
% Get pars of needed pargoups: Go through pargroups, if a pargroup needed:
% > check if present in parfile
% > do conversions if needed due to older file version
% > move unsupported parameters to unknown_pars
% Go through needed_pargroups and load if parfile contains it
for i_g = 1:length( needed_pargroups )
%If file doesnt/does contain such a pargroup
if isfield( parfile_allpars, needed_pargroups{i_g} )~=1
    log_msg_c{end+1} = sprintf('Parfile does not contain ''%s'' ',...
        needed_pargroups{i_g});
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check types of errfilt (whether supported)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Go through records of current pargroup
    % and check whether current gui version knows them
    curr_rnames = fieldnames( parfile_allpars.(needed_pargroups{i_g}));
    for i_r = 1:length( curr_rnames )
        if isfield( gui_default_pars.(needed_pargroups{i_g}),...
                curr_rnames{i_r} )==1
            loaded_pars.(needed_pargroups{i_g}).(curr_rnames{i_r}) =...
                parfile_allpars.(needed_pargroups{i_g}).(curr_rnames{i_r});
            %add to loaded pargroups - if not already there
            if isempty(find(strcmp(...
                    loaded_pargroups,needed_pargroups{i_g}),1))==1
                loaded_pargroups{end+1} = needed_pargroups{i_g};
            end
        else
            unknown_pars.(needed_pargroups{i_g}).(curr_rnames{i_r}) =...
                parfile_allpars.(needed_pargroups{i_g}).(curr_rnames{i_r});
        end
    end
end
end
if isempty( unknown_pars )~=1
    %Message:
    log_msg_c{end+1} = sprintf(...
        'unsupported parameters in par-file');
end





% ------------------------------------------------------------------------
% parfile: Read 
% ------------------------------------------------------------------------
function [ temp_pars, err_msg_1 ] = fcn_parfile_read_all_ascii(parfilename)
% Open file for read, (read lines to cellstr), load content to variable
% if err_msg_c not empty == no variable loaded
err_msg_1 = {};%error : one line!!!
temp_pars = struct;
%temp_strs = cell(0);%line contents as strings
fid = fopen( parfilename , 'rt');%read text
% Get the first line
tline = fgetl(fid);
% Check line (exclude major errors indtroduced by users) then read:
while ischar(tline)%feof(fid) == 0%
    % Discard commented parts:
    pos_sigcmt = strfind( tline, '//');
    if isempty(pos_sigcmt)~=1
        tline(pos_sigcmt:end)=[];
    end
    if isempty( tline ) ~=1
        % Syntax check (security relevant)
        % Allow only lines with a single assignment ('=')
        pos_sigeql = strfind( tline, '=');
        if isempty(pos_sigeql)==1%value assignment not present
            err_msg_1{1} = 'Parfile syntax error (''='') ';
            return
        elseif length(pos_sigeql)>1%second assignment invalid
            err_msg_1{1} = 'Parfile syntax error (''='') ';
            return
        end
        % Allow only lines with struct-variables (containing '.') before'='
        if isempty( strfind( tline(1:pos_sigeql), '.') )==1
            err_msg_1{1} = 'Parfile syntax error (''.'') ';
            return
        end
        % Do not allow ';' before '='
        if isempty( strfind( tline(1:pos_sigeql), ';') )~=1
            err_msg_1{1} = 'Parfile syntax error ('';'') ';
            return
        end
        % Do not allow 'eval(' before '='
        if isempty( strfind( tline(1:pos_sigeql), 'eval(') )~=1
            err_msg_1{1} = 'Parfile syntax error (''eval('') ';
            return
        end
        % Do not allow 'char(' before '='
        if isempty( strfind( tline(1:pos_sigeql), 'char(') )~=1
            err_msg_1{1} = 'Parfile syntax error (''char('') ';
            return
        end
        % Remove spaces before '=' to avoid commands separated by space
        tline(isspace(tline(1:pos_sigeql))) = [];
        % Get line contents as string
        % temp_strs{end+1} = tline;
        % Load content to variable:
        eval( ['temp_pars.'  tline ]);
    end
    %Get new line:
    tline = fgetl(fid);
end
fclose(fid);
% status = fclose(fid);%status=... show 0 if succesful
% if status==0
%     fprintf(' ... ASCII-file reading successfull: %s \n', hdrfilename);
% end






% ------------------------------------------------------------------------
% parfile: write 
% ------------------------------------------------------------------------
function [ err_msg_1 ] = fcn_parfile_write_ascii(...
    parfilename, pars2exp)%file_perm, 
% Write content of temp_pars to ascii file:
% 1. Convert content of temp_pars to matlab-evaluatable strings
% 2. Open file for write (with permission file_perm: append oder replace)
% 3. Write strings line-by-line to file
% if err_msg_1 not empty == no variable loaded
err_msg_1 = {};%error : one line!!!
% ------------------------------------------------------------
% ------------------------------------------------------------
%pargroups in temp_pars:
pgroupnames = fieldnames(pars2exp);
% convert pargoups to strings and print to file
fid = fopen( parfilename , 'w');%write - append or overwrite
fprintf( fid, '%s\r\n', '//');
%Pargroups one-by-one:
for i_g = 1:length( pgroupnames )
%start with a commented line:
fprintf( fid, '%s\r\n', '//');
% fprintf( '%s\r\n', '//');
%current pargroup:
curr_pgroup = pars2exp.(pgroupnames{i_g});
%get par_names & number of fields in curr_pargroup:
curr_parnames = fieldnames( curr_pgroup );
%get content of pars:
for i_p = 1:length( curr_parnames )
    %current par within pargroup:
    curr_par = curr_pgroup.(curr_parnames{i_p});
    if isstruct( curr_par )~=1
        %Convert variable content to string
        [curr_par_str,err_msg] = fcn_parfile_conv_var2str(curr_par);
        if isempty( err_msg )~=1%error occured at read
            err_msg_1{end+1} = err_msg{1};
            return
        end
        %Add string to text
        fprintf( fid, '%s.%s=%s;\r\n',...
            pgroupnames{i_g},curr_parnames{i_p},...
            curr_par_str);
    else%its a struct - get each field!!!
        %get fieldnames of struct:
        curr_par_fieldnames = fieldnames( curr_par );
        %get content of pars in struct:
        for p_i = 1:length( curr_par )
            for p_f = 1:length( curr_par_fieldnames )
                %Convert variable content to string
                [curr_par_str,err_msg] = fcn_parfile_conv_var2str(...
                    curr_par(p_i).(curr_par_fieldnames{p_f}));
                if isempty( err_msg )~=1%error occured at read
                    err_msg_1{end+1} = err_msg{1};
                    return
                end
                %Add string to text
                fprintf( fid, '%s.%s(%d).%s=%s;\r\n',...
                    pgroupnames{i_g},curr_parnames{i_p},...
                    p_i,curr_par_fieldnames{p_f},...
                    curr_par_str );
            end
        end
    end
end
end
% ------------------------------------------------------------
fclose(fid);
% % status = fclose(fid);%status=... show 0 if succesful
% % if status==0
% %     fprintf(' ... ASCII-file reading successfull: %s \n', hdrfilename);
% % end





% ------------------------------------------------------------------------
% parfile: write preprocess
% ------------------------------------------------------------------------
function [ curr_par_str, err_msg ] = fcn_parfile_conv_var2str( curr_par )
%Converts variable content to str
%Supported: char, double, cell (Not supported: struct!)
err_msg = {};
% string conversion data class dependent:
curr_cl = class( curr_par );
if strcmp('char',curr_cl)==1
    curr_par_str = sprintf('''%s''',...
        curr_par);
elseif strcmp('double',curr_cl)==1
    curr_par_str = sprintf('%s',...
        mat2str( curr_par ) );
elseif strcmp('cell',curr_cl)==1
    curr_par_str = '{';
    for i_r = 1:size( curr_par, 1)
        if i_r~=1
            curr_par_str(end+1)=';';
        end
        for i_c = 1: size( curr_par, 2)
            if i_c~=1
                curr_par_str(end+1)=',';
            end
            curr_par_str = sprintf('%s''%s''',...
                curr_par_str,curr_par{i_r,i_c});
            % curr_par_str =...
            %     strcat(curr_par_str,curr_par{i_r,i_c});
        end
    end
    curr_par_str(end+1) = '}';
    curr_par_str = sprintf('%s',...
        curr_par_str);
else
    curr_par_str = '';
    err_msg = {'Paramter set not supported!'};
end




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 
function chapter_16_other_functions
% 
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------
% Calculates blackmannharris window
% ------------------------------------------------------------
function [ bhwin ] = fcn_blackmanharris( N )
a0 = 0.35875;
a1 = 0.48829;
a2 = 0.14128;
a3 = 0.01168;
n = 0:N-1;
bhwin = a0 - ...
    a1*cos(2*pi*n/(N-1)) + ...
    a2*cos(4*pi*n/(N-1)) - ...
    a3*cos(6*pi*n/(N-1));



% ------------------------------------------------------------
% Calculates rounded intervall limits for min and max!
% ------------------------------------------------------------
function [ rounded_lims ] = fcn_sokoray_fit_lims2data(min_i, max_i)
% Author: sokoray (bsokoray@gmail.com)
rounded_lims = nan(size( min_i, 1),2 );
for d_i = 1:size( min_i, 1)
    %Exponent-range of data:
    int_omg_pow = floor(log10( max_i(d_i) - min_i(d_i) ));
    if int_omg_pow<0 %interval_omg < 1
        %Order of magnitude of rounding 
        edges_omg = 10^(floor(log10( max_i(d_i) - min_i(d_i) ))+1);
    elseif int_omg_pow<=1%interval_omg = 1 or 10
        %Order of magnitude of rounding
        edges_omg = 10^(floor(log10( max_i(d_i) - min_i(d_i) )));
    else%interval_omg >= 100
        %Order of magnitude of rounding
        edges_omg = 10^(floor(log10( max_i(d_i) - min_i(d_i) ))-1);
    end
    %Rounded lower & upper limits:
    rounded_lims(d_i,:) = [...
        edges_omg*floor( min_i(d_i) / edges_omg ),...
        edges_omg*ceil(  max_i(d_i) / edges_omg )];
end























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 
function chapter_12_importing_data
% 
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% -------------------------------------------------------------------------
% Import time series data
% -------------------------------------------------------------------------
function [ curr_rawdset, curr_log, addlist_err_loading, addlist_parsloaded,...
    addlist_wrn_atload, addlist_wrn_parsupport,addlist_wrn_enucoords] =...
    gui_addfile_check_load_files( Pathname, curr_filename, handles )
%Collect file by file fields of rawts as:
% curr_rawts.fpath
% curr_rawts.filename
% curr_rawts.par_pars
% curr_rawts.par_groups
% curr_rawts.hdr_nfo
% curr_rawts.type
% curr_rawts.data
% curr_rawts.dat_props
%Initialize:
temp_log_case_loaded = {};%collected logs for files that are ok to be loaded
curr_log = {};%log. "temp_log_loading" is included at the end
addlist_parsloaded     = false;
addlist_wrn_atload     = false;
addlist_wrn_parsupport = false;
addlist_wrn_enucoords  = false;
addlist_err_loading    = false;
% ------------------------------------------------------------
% Filename info
% ------------------------------------------------------------
% if iscell( curr_filename )==1
%     Filename{1} = 1
%     i=1;
% end
curr_rawdset.fpath    = Pathname;
curr_rawdset.filename = curr_filename;
[path,fname,ext] = fileparts( curr_filename );
% ------------------------------------------------------------
% Check whether .dat file (only .dat files allowed)
% It has to be a dat-file, otherwise viper-exported ascii-files
% next to .hdr-file of the original .dat produce error
% ------------------------------------------------------------
if strcmp(ext,'.dat')==0
    curr_log{end+1,1} =...
        sprintf('!!! Error: %s NOT loaded! (%s)',...
        curr_filename, 'Not supported file extension');
    addlist_err_loading = 1;
    curr_rawdset = struct([]);
    return
end
% ------------------------------------------------------------
% Read par file contents (.par) - if exist
% ------------------------------------------------------------
if exist( [Pathname fname '.par'],'file')==2
    pargroups = fieldnames( handles.ts_defaults.p );%== all groups
    [ temp_pars, temp_pargroups, unknown_p, log_msg_c ]=...
        fcn_parfile_load( [Pathname fname '.par'], pargroups, ...
        handles.ts_defaults.p,...
        handles.viper_props.errf_types_supported(:,1));
    % ---------------------------------------------------------
    %If error loading par
    %else check whether warnings at par loading
    if isempty( temp_pars )==1
        temp_log_case_loaded{end+1,1} =...
            sprintf('    ! Error loading par-file');
        addlist_wrn_atload = true;
    else
        temp_log_case_loaded{ end+1, 1 } =...
            sprintf('    # par-file loaded');
        %unsupported pars:
        if isempty( unknown_p )~=1
            temp_log_case_loaded{end+1,1} =...
                sprintf('      ! ignoring unsupported parameters!');
            addlist_wrn_parsupport = true;
        end
    end
    %further warnings and errors
    if isempty( log_msg_c )~=1
        for i_r=1:length(log_msg_c)
            temp_log_case_loaded{end+1,1} =...
                sprintf('      * %s', log_msg_c{i_r});
        end
    end
else
    temp_pars = [];
    temp_pargroups = [];
end
curr_rawdset.par_pars   = temp_pars;
curr_rawdset.par_groups = temp_pargroups;
% ------------------------------------------------------------
% Read hdr file contents (.hdr) - if exist
% ------------------------------------------------------------
if exist([Pathname fname '.hdr'],'file')==2
    [ temp_hdr_read ]= fcn_addfile_hdr_read_vno_vec(...
        [ Pathname, fname, '.hdr'], handles.ts_defaults.hdr_nfo );
    %check if vec or vno
    if isempty( temp_hdr_read.serial_no )==0 && (...
            strcmp( temp_hdr_read.serial_no(1:3), 'VNO' )==1 ||...
            strcmp( temp_hdr_read.serial_no(1:3), 'VEC' )==1 ) && ...
            isempty( temp_hdr_read.sampling_rate )==0
        %Store hdr content
        temp_hdr = temp_hdr_read;
        temp_log_case_loaded{ end+1, 1 } =...
            sprintf('    # hdr-file loaded');
    else%not vec or vno
        temp_hdr.existing = false;
    end
    %check supported coord system:
    if strcmp( temp_hdr_read.coord_system, 'XYZ' )~=1 &&...
            strcmp( temp_hdr_read.coord_system, 'ENU' )~=1
        curr_log{end+1,1} =...
            sprintf('!!! Error: %s NOT loaded! (%s)',...
            curr_filename, 'Not supported probe coordinate system');
        addlist_err_loading = 1;
        curr_rawdset = struct([]);
        return
    elseif strcmp( temp_hdr_read.coord_system, 'ENU' )==1
        addlist_wrn_enucoords = true;
        temp_log_case_loaded{ end+1, 1 } =...
            sprintf('     * Vectors in ENU coordinates!');
    end
else
    temp_hdr.existing = false;
end
curr_rawdset.hdr_nfo = temp_hdr;
% ------------------------------------------------------------
% Determine OR ask data type: 'vno' | 'vec' | 'usr'
% ------------------------------------------------------------
% If parfile exists and loaded - get from there
% Else: if hdrfile exists and loaded - get from there
% Else: ask user
if isempty( find( strcmp( temp_pargroups, 'dat_props'),1 ))~=1
    temp_type = temp_pars.dat_props.type;
elseif temp_hdr.existing==1
    if strcmp( temp_hdr.serial_no(1:3), 'VNO' )==1
        temp_type = 'vno';
    end
    if strcmp( temp_hdr.serial_no(1:3), 'VEC' )==1
        temp_type = 'vec';
    end
else
    input = questdlg(sprintf(...
        'Please specify data type for the file:\n%s%s',...
        Pathname, curr_filename),...
        'Data type?','Vectrino','Vector','User defined',...
        'User defined');%'Vector'
    if strcmp(input,'Vectrino')==1
        temp_type = 'vno';
    elseif strcmp(input,'Vector')==1
        temp_type = 'vec';
    elseif strcmp(input,'User defined')==1
        temp_type = 'usr';
    else % strcmp(input,'Skip file')==1 || empty
        curr_log{end+1,1} =...
            sprintf('!!! Error: %s NOT loaded! (%s)',...
            curr_filename, 'Undefined data type');
        addlist_err_loading = 1;
        curr_rawdset = struct([]);
        return
    end
end
curr_rawdset.type = temp_type;
% ------------------------------------------------------------
% Read file content -> check whether ascii
% * so far supporting only rectangular numeric data without header
% ------------------------------------------------------------
temp_file_ascii = 0;
try%load == ascii file must contain a rectangular numeric table
    temp_data = dlmread( [Pathname, curr_filename ] );%,'',[0 0 3 4]);
    temp_file_ascii = 1;
catch%not ascii or not rectangular table: error!
    temp_log_ext = 'Non-ascii formats are not supported';
end
% if ascii: Check whether data & probe type consistent (num columns)
if temp_file_ascii == 1
    if strcmp( temp_type, 'vno' )==1 && size( temp_data, 2 )<18
        input = questdlg(sprintf(...
            'Not enough columns in file for Vectrino data\n%s%s',...
            Pathname, curr_filename),...
            'Wrong data type?',...
            'Load as ''User defined'' data','Skip file',...
            'Skip file');%'Vector'
        if strcmp(input,'Load as ''User defined'' data')==1
            temp_type = 'usr';
        else % strcmp(input,'Skip file')==1 || empty
            curr_log{end+1,1} =...
                sprintf('!!! Error: %s NOT loaded! (%s)',...
                curr_filename, 'Insufficient number of columns');
            addlist_err_loading = 1;
            curr_rawdset = struct([]);
            return
        end
    end
    if strcmp( temp_type, 'vec' )==1 && size( temp_data, 2 )<14
        input = questdlg(sprintf(...
            'Not enough columns in file for Vector data\n%s%s',...
            Pathname, curr_filename),...
            'Wrong data type?',...
            'Load as ''User defined'' data','Skip file',...
            'Skip file');%'Vector'
        if strcmp(input,'Load as ''User defined'' data')==1
            temp_type = 'usr';
        else % strcmp(input,'Skip file')==1 || empty
            curr_log{end+1,1} =...
                sprintf('!!! Error: %s NOT loaded! (%s)',...
                curr_filename, 'Insufficient number of columns');
            addlist_err_loading = 1;
            curr_rawdset = struct([]);
            return
        end
    end
    if strcmp( temp_type, 'usr' )==1 && size( temp_data, 2 )<4
        dlg_h = errordlg(sprintf(...
            'Not enough columns in file for time-series data\n%s%s',...
            Pathname, curr_filename),...
            'Wrong data file?','Skip file');%'Vector'
        waitfor(dlg_h)
        curr_log{end+1,1} =...
            sprintf('!!! Error: %s NOT loaded! (%s)',...
            curr_filename, 'Insufficient number of columns');
        addlist_err_loading = 1;
        curr_rawdset = struct([]);
        return
    end
else
    curr_log{end+1,1} =...
        sprintf('!!! Error: %s NOT loaded! (%s)',...
        curr_filename, temp_log_ext);
    addlist_err_loading = 1;
    curr_rawdset = struct([]);
    return
end
curr_rawdset.data = temp_data;
temp_log_case_loaded{ end+1, 1 } =...
    sprintf('    # dat-file loaded');
% ------------------------------------------------------------
%  Get data properties parameters - needed to load data
%  %dat_props: .type .column_t .column_u .frq .probe_cfg .probe_nvr
% ------------------------------------------------------------
% if parfile loaded -> get from parset loaded from parfile
% else (parfile not loaded) -> findout from .dat, .hdr or user
%      -> If user input -> check consistency
if isempty( find( strcmp( temp_pargroups, 'dat_props'),1 ))==0
    addlist_parsloaded = true;%mark that paramter has been loaded:
    %Create messages:
    temp_log_ext = 'from parfile';%text for msg-box: source of par
    temp_log_case_loaded{end+1,1} = sprintf(...
        '    * dat columns t:%d, u:%d (%s); fs: %d Hz (%s)',...
        temp_pars.dat_props.column_t,...
        temp_pars.dat_props.column_u,...
        temp_log_ext, temp_pars.dat_props.frq, temp_log_ext);
    %Store:
    temp_dat_props = temp_pars.dat_props;
else
    % get based on .dat-file (temp_data) or hdr_info or user input:
    [temp_dat_props, input_v_unit, temp_log] =...
        fcn_addfile_findout_dat_props(...
        temp_type, temp_data, temp_hdr, curr_filename,...
        handles.viper_props.userdef.dat_props.column_t,...
        handles.viper_props.userdef.dat_props.column_u );
    temp_log_case_loaded{end+1,1} = temp_log;
end
curr_rawdset.dat_props = temp_dat_props;
% Update log at success:
curr_log{end+1,1} =...
    sprintf('  > Loading: %s',curr_filename);
%Add temporary log
if isempty( temp_log_case_loaded )~=1
    curr_log = vertcat( curr_log, temp_log_case_loaded );
end









% -------------------------------------------------------------------------
% Import time series data
% -------------------------------------------------------------------------
function [ curr_ts, curr_fpath, curr_fname, curr_fullname, curr_dispname,...
    curr_log, addlist_parsloaded, addlist_err_errfiltering ] =...
    gui_addfile_import_ts_data(...
    curr_raw_dset, addlist_parsloaded, handles )
% Import data - nothing to ask
% handles.ts obtains here the final fields based on curr_ts fields
%            further fields must not be added later in viper!
curr_fpath = {};
curr_fname = {};
curr_fullname = {};
curr_dispname = {};
curr_log = {};
addlist_err_errfiltering = false;
%Log:
curr_log{end+1,1} =...
    sprintf('  > Processing %s',...
    curr_raw_dset.filename );

% ------------------------------------------------------------
% I. REQUIRED FROM RAW
% ------------------------------------------------------------
% 1. Initializing TS variables with defaults
% Initialize intl data:
curr_ts.intl        = handles.ts_defaults.intl;
% Further ts-variables - defaults
curr_ts.seri_data   = handles.ts_defaults.seri_data;
curr_ts.seri_enable = handles.ts_defaults.seri_enable;
curr_ts.valu_data   = handles.ts_defaults.valu_data;
curr_ts.valu_enable = handles.ts_defaults.valu_enable;

% ------------------------------------------------------------
% 2. IMPORTING hdr_nfo (built at loading files)
% Store loaded hdr-data
curr_ts.hdr_nfo = curr_raw_dset.hdr_nfo;

% ------------------------------------------------------------
% 3. PARAMETERSET BUILT HERE based on available loaded pars
% Store dat_props to handles:
curr_ts.p.dat_props = curr_raw_dset.dat_props;
% Store coordsys parameters - needed if data has to be rotated
% (corrdinates & direction) from parfile
if isempty( find( strcmp( curr_raw_dset.par_groups, 'coordsys'),1 ))~=1
    %store raw to handles:
    curr_ts.p.coordsys = curr_raw_dset.par_pars.coordsys;
    %adjust units for internal use:
    curr_ts.intl.coordsys = curr_ts.p.coordsys;
    curr_ts.intl.coordsys.unit = handles.ts_defaults.intl.coordsys.unit;
    curr_ts.intl.coordsys.file_coords = ...
        curr_ts.p.coordsys.file_coords*handles.gui.coordunitscale;
    curr_ts.intl.coordsys.measpt_offset = ...
        curr_ts.p.coordsys.measpt_offset*handles.gui.coordunitscale;
    curr_ts.intl.coordsys.measpt_coords = ...
        curr_ts.p.coordsys.measpt_coords*handles.gui.coordunitscale;
    if isempty(find(isnan(curr_ts.intl.coordsys.measpt_coords),1))==1
        curr_ts.enable_coords = {'On'};
    else
        curr_ts.enable_coords = {'Off'};
    end
    addlist_parsloaded = true;%mark that paramter has been loaded:
else%no coordsys from parfile -> set defaults
    curr_ts.p.coordsys         = handles.ts_defaults.p.coordsys;
    curr_ts.intl.coordsys      = handles.ts_defaults.p.coordsys;
    curr_ts.intl.coordsys.unit = handles.ts_defaults.intl.coordsys.unit;
    curr_ts.enable_coords = {'Off'};
end
                
% ------------------------------------------------------------
% 4. CREATE TIME variable: Generate/load time data
if strcmp( curr_raw_dset.type, 'vno' )==1 ||...
        strcmp( curr_raw_dset.type, 'vec' )==1
    % VNO & VEC: time column does not have the needed accuracy!
    % AND if missing time-steps -> positions to correct
    % This can be done using ensemble (sample) counter of vectrino!
    if strcmp( curr_raw_dset.type, 'vno' )==1
        e_col = curr_ts.p.dat_props.column_u-2;
    end
    if strcmp( curr_raw_dset.type, 'vec' )==1
        if isempty( curr_raw_dset.hdr_nfo.ensemble_column )~=1
            e_col = str2double( curr_raw_dset.hdr_nfo.ensemble_column );
        else
            e_col = 2;
        end
    end
    % ------------------------------------------------------------
    ens_cnt = curr_raw_dset.data(:,e_col);
    %Accurate time values - independent from missing samples
    curr_ts.raw_t =...
        (0: 1/curr_ts.p.dat_props.frq : (ens_cnt(end)-ens_cnt(1))/...
        curr_ts.p.dat_props.frq)';
    %Data length: data_length = length(ens_cnt);
    %raw_tslength: if missing samples => raw_tslength > length(ens_cnt)
    raw_tslength = length( curr_ts.raw_t );
    %Positions of data samples within corrected time-series
    pos_spl = ens_cnt-(ens_cnt(1)-1);
    curr_ts.raw_vals = false(raw_tslength,1);
    curr_ts.raw_vals( pos_spl ) =1;
else
    %Time colums is required in dat of (non-vec|non-vno)-> use it!
    curr_ts.raw_t = curr_raw_dset.data(:,curr_ts.p.dat_props.column_t);
    raw_tslength = length( curr_ts.raw_t );
    pos_spl = 1:raw_tslength;
    curr_ts.raw_vals = true( raw_tslength, 1);
end
%Valid raw-positions:
curr_ts.vals = curr_ts.raw_vals;

% ------------------------------------------------------------
% 5. IMPORT MEASURED DATA:
%   > valid samples: (a) h=nan(raw_tslength,1) (b) h(pos_samples)=data
%   > depending on probe type: leave columns nan
% ------------------------------------------------------------
%Allocate variable for [u,v,w,w1,w2](case non-vno: [w1,w2] == nan)
curr_ts.raw_veldata = nan( raw_tslength, 5 );
%Allocate snr-cor variable: [1:4,5:8]: for non-vno: [4,8]==nan
if strcmp( curr_raw_dset.type, 'vno' )==1 || ...
        strcmp( curr_raw_dset.type, 'vec' )==1
    curr_ts.raw_snrcor = nan( raw_tslength, 8 );
else
    curr_ts.raw_snrcor = [];
end
%import variables:
u_col = curr_ts.p.dat_props.column_u;
vsc = 1000/handles.gui.units_id;%velocity scaling
%u- and v-velocities:
curr_ts.raw_veldata(pos_spl,1) = vsc*curr_raw_dset.data(:,u_col);
curr_ts.raw_veldata(pos_spl,2) = vsc*curr_raw_dset.data(:,u_col+1);
%w-velocities:
curr_ts.raw_veldata(pos_spl,3) = vsc*curr_raw_dset.data(:,u_col+2);
%VNO: w-velocities % comp. dir at load (parallel or normal to trasmiter)
if strcmp( curr_raw_dset.type, 'vno' )==1
    if strcmp(curr_ts.p.dat_props.probe_cfg,...
            'down-looking')==1
        curr_ts.raw_veldata(pos_spl,3) =...
            vsc*(curr_raw_dset.data(:,u_col+2) + curr_raw_dset.data(:,u_col+3 ) )/2;
        %no change
        curr_ts.intl.cmp_probe_np = {'n' 'n' 'p'};
    else%side-looking -> w2=w1
        curr_ts.raw_veldata(pos_spl,3) =...
            vsc*curr_raw_dset.data(:,u_col+2);
        curr_ts.intl.cmp_probe_np = {'n' 'p' 'n'};
    end
    %raw w1 & w2:
    curr_ts.raw_veldata(pos_spl,4) = vsc*curr_raw_dset.data(:,u_col+2);
    curr_ts.raw_veldata(pos_spl,5) = vsc*curr_raw_dset.data(:,u_col+3);
end
if strcmp( curr_raw_dset.type, 'vec' )==1
    curr_ts.intl.cmp_probe_np = {'n' 'n' 'p'};
end
%SNR & COR
if strcmp( curr_raw_dset.type, 'vno' )==1
    %SNR:
    curr_ts.raw_snrcor(pos_spl,1) = curr_raw_dset.data(:,u_col+8);
    curr_ts.raw_snrcor(pos_spl,2) = curr_raw_dset.data(:,u_col+9);
    curr_ts.raw_snrcor(pos_spl,3) = curr_raw_dset.data(:,u_col+10);
    curr_ts.raw_snrcor(pos_spl,4) = curr_raw_dset.data(:,u_col+11);
    %Correlations:
    curr_ts.raw_snrcor(pos_spl,5) = curr_raw_dset.data(:,u_col+12);
    curr_ts.raw_snrcor(pos_spl,6) = curr_raw_dset.data(:,u_col+13);
    curr_ts.raw_snrcor(pos_spl,7) = curr_raw_dset.data(:,u_col+14);
    curr_ts.raw_snrcor(pos_spl,8) = curr_raw_dset.data(:,u_col+15);
end
if strcmp( curr_raw_dset.type, 'vec' )==1
    %SNR:
    curr_ts.raw_snrcor(pos_spl,1) = curr_raw_dset.data(:,u_col+6);
    curr_ts.raw_snrcor(pos_spl,2) = curr_raw_dset.data(:,u_col+7);
    curr_ts.raw_snrcor(pos_spl,3) = curr_raw_dset.data(:,u_col+8);
    %Correlations:
    curr_ts.raw_snrcor(pos_spl,5) = curr_raw_dset.data(:,u_col+9);
    curr_ts.raw_snrcor(pos_spl,6) = curr_raw_dset.data(:,u_col+10);
    curr_ts.raw_snrcor(pos_spl,7) = curr_raw_dset.data(:,u_col+11);
end

% ------------------------------------------------------------
% 6. ROTATE UVW components: ONLY in RAW!(if parfile requires)
% This is before calculating error_filtering and statistics
% ------------------------------------------------------------
if strcmp( curr_ts.intl.coordsys.dir_probe_u, '+x')~=1
    dir_old = '+x';
    dir_new = curr_ts.intl.coordsys.dir_probe_u;
    only_raw = 1;
    %rotate u-v vector in affected variables
    %@addfile:
    %  - rotate only raw_veldata_uvw & intl
    %  - parameters are contained in the file rotated=correctly
    %  - seri & valu and errfilt res are generated LATER correctly
    [ curr_ts.raw_veldata,...
        unused1,...
        unused2,...
        curr_ts.intl,...
        unused3,...
        unused4,...
        alpha ] =...
        fcn_calc_rotate_uv(...
        only_raw, dir_old, dir_new,...
        curr_ts.raw_veldata,...
        [],...
        [],...
        curr_ts.intl,...
        [],...
        [],...
        [],...
        [] );
    %Create log:
    curr_log{end+1,1} =...
        sprintf('    # velocity vector rotated [ u_probe -> %s ]',...
        curr_ts.intl.coordsys.dir_probe_u );
end

% ------------------------------------------------------------
% II. OPTIONAL FROM RAW - if exist: errfilt
% ------------------------------------------------------------
% ------------------------------------------------------------
% 7. Assign err.filt pars from parfile (if exists) or use defaults
% ------------------------------------------------------------
if isempty( find( strcmp( curr_raw_dset.par_groups, 'errfilt'),1))~=1
    %Check whether errfilt pars compatible with data:
    % in principle this is unneded, as .par-files are verified before save
    % but if user creates .par file by copying an other one?!
    [ err_msg ] = fcn_errfiltpars_check_validity(...
        curr_raw_dset.par_pars.errfilt,...
        curr_ts.p.dat_props, isfield(handles,'dvl') );
    if isempty( err_msg )==0
        %Error -> do not use any filter methods!
        curr_ts.p.errfilt = [];
        curr_ts.p.errfilt.used_velcomp =...
            handles.ts_defaults.p.errfilt.used_velcomp;
        curr_ts.p.errfilt.used_types = {};
        %Create log:
        curr_log{end+1,1} =...
            sprintf('    ! err.filtering parameters from par-file not compatible with data');
        for i_e=1:length(err_msg)
            curr_log{end+1,1} =...
                sprintf('      * %s', err_msg{i_e});
        end
        curr_log{end+1,1} =...
            sprintf('      * assigned parameters without err.filtering');
    else
        %assign:
        curr_ts.p.errfilt = curr_raw_dset.par_pars.errfilt;
        addlist_parsloaded = true;%mark that paramter has been loaded:
        %Create log:
        curr_log{end+1,1} =...
            sprintf('    # err.filtering parameters assigned from par-file');
    end    
else%no errfilt_paramemeters in parfile -> set defaults
    %Use err.filtering parameter defaults
    curr_ts.p.errfilt.used_velcomp =...
        handles.ts_defaults.p.errfilt.used_velcomp;
    curr_ts.p.errfilt.used_types = {};
    %Use those that are ok for dat_type -> Get enabled ones for dat_type:
    ef_vals = strcmpi(...
        handles.viper_props.errf_types_enable_matrix(:,1),...
        curr_raw_dset.type );
    errf_types_enabled = ...
        handles.viper_props.errf_types_supported( ...
        handles.viper_props.errf_types_enable_matrix{ef_vals ,2}, 1);
    %go through enabled ones and if selected in defaults -> copy it
    for i_ef = 1:length( errf_types_enabled )
        ef_num = length( curr_ts.p.errfilt.used_types ) + 1;
        ef_mrk = strcmpi(...
            handles.ts_defaults.p.errfilt.used_types,...
            errf_types_enabled{ i_ef } );
        if sum( ef_mrk )>0
            curr_ts.p.errfilt.used_types{ ef_num } = ...
                errf_types_enabled{ i_ef };
            curr_ts.p.errfilt.typ( ef_num ) = ...
                handles.ts_defaults.p.errfilt.typ( ef_mrk );
        end
    end
    %Check errfilt pars for validity == applicability to selected file:
    [ err_msg ] = fcn_errfiltpars_check_validity(...
        curr_ts.p.errfilt, curr_ts.p.dat_props, isfield(handles,'dvl') );
    if isempty( err_msg )==0
        %Error -> do not use any filter methods!
        curr_ts.p.errfilt = [];
        curr_ts.p.errfilt.used_velcomp =...
            handles.ts_defaults.p.errfilt.used_velcomp;
        curr_ts.p.errfilt.used_types = {};
        %Create log:
        curr_log{end+1,1} =...
            sprintf('    ! err.filtering parameter defaults not compatible with data');
        for i_e=1:length(err_msg)
            curr_log{end+1,1} =...
                sprintf('      * %s', err_msg{i_e});
        end
        curr_log{end+1,1} =...
            sprintf('      * assigned parameters without err.filtering');
    else
        %Create log:
        curr_log{end+1,1} =...
            sprintf('    # err.filtering parameter defaults assigned');
    end
end
% -----------------------------------------------------------------
% DVL: After Rotating! Before Errfilt! TODO: check errfilt support
% -----------------------------------------------------------------
if isfield( handles, 'dvl')==1
    [ ~, curr_ts ] = viper_dvl(...
        'dvl_addfile', handles, { curr_raw_dset.type, curr_ts } );
end

% ------------------------------------------------------------
% III. DERIVED ELEMENTS: %.raw_corsnr_xt
% ------------------------------------------------------------

% ------------------------------------------------------------
% 8. Calculate errfilt lambda
% ------------------------------------------------------------
% not needed as done after crop!

% ------------------------------------------------------------
% 9. Generate quality parameter variables (min and avg COR & SNR)
% ------------------------------------------------------------
%Preallocate variables:
if strcmp( curr_raw_dset.type, 'vno' )==1 || ...
        strcmp( curr_raw_dset.type, 'vec' )==1
    curr_ts.raw_corsnr_xt = nan( raw_tslength, 4 );
else
    curr_ts.raw_corsnr_xt = [];
end
if strcmp( curr_raw_dset.type, 'vno' )==1
    %Minimal correlations among the beams
    curr_ts.raw_corsnr_xt(:,1) = min( horzcat(...
        curr_ts.raw_snrcor(:,5),...
        curr_ts.raw_snrcor(:,6),...
        curr_ts.raw_snrcor(:,7),...
        curr_ts.raw_snrcor(:,8) ), [], 2);
    %Average correlations among the beams
    curr_ts.raw_corsnr_xt(:,2) = (...
        curr_ts.raw_snrcor(:,5) +...
        curr_ts.raw_snrcor(:,6) +...
        curr_ts.raw_snrcor(:,7) +...
        curr_ts.raw_snrcor(:,8) )/4;
    %Minimal SNR among the beams
    curr_ts.raw_corsnr_xt(:,3) = min( horzcat(...
        curr_ts.raw_snrcor(:,1),...
        curr_ts.raw_snrcor(:,2),...
        curr_ts.raw_snrcor(:,3),...
        curr_ts.raw_snrcor(:,4) ), [], 2);
    %Average SNR among the beams
    curr_ts.raw_corsnr_xt(:,4) = (...
        curr_ts.raw_snrcor(:,1) +...
        curr_ts.raw_snrcor(:,2) +...
        curr_ts.raw_snrcor(:,3) +...
        curr_ts.raw_snrcor(:,4) )/4;
end
if strcmp( curr_raw_dset.type, 'vec' )==1
    %Minimal correlations among the beams
    curr_ts.raw_corsnr_xt(:,1) = min( horzcat(...
        curr_ts.raw_snrcor(:,5),...
        curr_ts.raw_snrcor(:,6),...
        curr_ts.raw_snrcor(:,7) ), [], 2);
    %Average correlations among the beams
    curr_ts.raw_corsnr_xt(:,2) = (...
        curr_ts.raw_snrcor(:,5) +...
        curr_ts.raw_snrcor(:,6) +...
        curr_ts.raw_snrcor(:,7) )/3;
    %Minimal SNR among the beams
    curr_ts.raw_corsnr_xt(:,3) = min( horzcat(...
        curr_ts.raw_snrcor(:,1),...
        curr_ts.raw_snrcor(:,2),...
        curr_ts.raw_snrcor(:,3) ), [], 2);
    %Average SNR among the beams
    curr_ts.raw_corsnr_xt(:,4) = (...
        curr_ts.raw_snrcor(:,1) +...
        curr_ts.raw_snrcor(:,2) +...
        curr_ts.raw_snrcor(:,3) )/3;
end

% ------------------------------------------------------------
% IV. PROCESSED DATA
% ------------------------------------------------------------

% ------------------------------------------------------------
% 10. ERROR FILTERING of raw data (to get valid data)
% ------------------------------------------------------------
[ efilt_res, err_msg ] = fcn_errfilt__main(...
    0, curr_ts.p.errfilt,...
    curr_ts.p.dat_props,...
    curr_ts.raw_veldata,...
    curr_ts.raw_vals,...
    curr_ts.raw_t,...
    curr_ts.raw_corsnr_xt,...
    isfield( handles, 'dvl'));
if isempty( err_msg )~=1
    % ---------------------------------------------------------
    %Error -> not store:
    curr_log{end+1,1} =...
        sprintf('    ! err.filtering failed: %s',...
        'invalid parameters');
    for i_e=1:length(err_msg)
        curr_log{end+1,1} =...
            sprintf('    * %s', err_msg{i_e});
    end
    %Create message:
    addlist_err_errfiltering = true;
    %raw_vals as .vals stays unchanged:
    %curr_ts.vals = curr_ts.raw_vals;
else
    % ---------------------------------------------------------
    %Store new efilt_res:
    curr_ts.efilt_res = efilt_res;
    % At this point: seri_enable: everything off
    % ---------------------------------------------------------
    % Store NEW vals due to errorfiltering:
    curr_ts.vals = curr_ts.efilt_res.valmarker(:,1);
    %Create message:
    curr_log{end+1,1} =...
        sprintf('    # err.filtering done');
end

% ------------------------------------------------------------
% 11. Calculate COR & SNR distribution for valid data
% > if no err.filtering done: uses rawvals
% ------------------------------------------------------------
curr_ts.corsnr_dist = [];
if strcmp( curr_raw_dset.type, 'vno' )==1 || ...
         strcmp( curr_raw_dset.type, 'vec' )==1
    %Calculate quality-plots:
    [ temp_corsnr, temp_boxp ] =...
        fcn_calc_corsnrdist_vals_adv_sokoray(...
        curr_ts.vals,...
        curr_ts.raw_corsnr_xt );
    %Correlation distrib.+cumulative & SNR distrib.+cumulative:
    curr_ts.corsnr_dist = temp_corsnr.corsnr_dist;
    %Boxplot data:
    curr_ts.valu_data( 19:23 ) = temp_boxp(1:5);
    %Enable
    curr_ts.valu_enable( 19:23 )=1;%COR-distrib boxplot
end

% ------------------------------------------------------------
% 12. Calculate velocity histograms
% > if no err.filtering done: uses rawvals
% ------------------------------------------------------------
[ curr_ts.hist.p_per_bin,...
    curr_ts.hist.vel_of_bin ] =...
    fcn_calc_histograms_of_vals(...
    curr_ts.raw_veldata,...
    curr_ts.vals, ...
    curr_ts.p.errfilt.used_velcomp );

% -----------------------------------------------------------------
% 13. Calc. err.filtered fluctuation & basic statistics using .vals:
% -----------------------------------------------------------------
[ curr_ts.seri_data,...
    curr_ts.seri_enable,...
    curr_ts.intl.mean_correctn,...
    curr_ts.intl.used_ts_seg,...
    new_f_seri_x,...
    new_f_valu_data,...
    new_f_valu_enable,...
    temp_f_valuupdate ] =...
    fcn_calc_statistical_func_n_valu(...
    curr_ts.raw_veldata,...
    curr_ts.vals,...
    curr_ts.p,...
    curr_ts.raw_t,...
    curr_ts.raw_corsnr_xt,...
    handles.ts_defaults,...
    size( curr_ts.valu_data ));
%Store additional seri_x
curr_ts.seri_x.t    = new_f_seri_x.t;
curr_ts.seri_x.lag  = new_f_seri_x.lag;
curr_ts.seri_x.tlag = new_f_seri_x.tlag;
curr_ts.seri_x.f    = new_f_seri_x.f;
%Store newly calculated valu (non-nan values)
curr_ts.valu_data( temp_f_valuupdate==1 ) =...
    new_f_valu_data( temp_f_valuupdate==1 );
curr_ts.valu_enable( temp_f_valuupdate==1 ) =...
    new_f_valu_enable( temp_f_valuupdate==1 );
%Update log - if bad data interpolation occured:
if isempty(find(isnan( curr_ts.intl.mean_correctn )==0,1))==0
    curr_log{end+1,1} =...
        sprintf('      * u-v-w corrected by [%2.2f, %2.2f, %2.2f] %s%s',...
        curr_ts.intl.mean_correctn,...
        'in order to get zero mean fluctuation after interpolation',...
        'for internal calculations ');
end
curr_log{end+1,1} =...
    sprintf('    # statistical features calculated');

% ------------------------------------------------------------
% V. Store filename data
% ------------------------------------------------------------
[~,fname,~] = fileparts( curr_raw_dset.filename );
curr_fpath{end+1,1}    = curr_raw_dset.fpath;
curr_fname{end+1,1}    = fname;
curr_fullname{end+1,1} = [ curr_raw_dset.fpath, curr_raw_dset.filename ];
curr_dispname{end+1,1} = [ curr_raw_dset.filename '   [' curr_raw_dset.fpath ']' ];


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 
function chapter_13_error_filtering
% 
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -------------------------------------------------------------------------
% ADV Err.filt main function: check pars and determine valid samples 
% by applying selected err.filters
% -------------------------------------------------------------------------
function [ efilt_res, err_msg ] = fcn_errfilt__main(...
    just_testing, errfilt_pars, curr_f_datprops, curr_f_raw_veldata,...
    curr_f_raw_vals, curr_f_raw_t, curr_f_raw_corsnr_xt, use_dvl)
err_msg = {};%add msg: err_msg{end+1,1} = ' text ';
%%%% 1. - Not checking parameter validity!
%%%% 2. Then creates default efilt_res output == each sample valid
%%%% 3. Then applies each filter in stored paramaterset
%%%%    *if no filter selected -> no change in default efilt_res
%%%%    *else: each filter types creates a valmarker set
%%%%           *after each filter: valmarkers are combined
%%%%           *erroneous positions are also collected for each type
% > valmarker = valid marker (same size as raw_valpos, 1 at valid values
% > errmarker == erroneous position marker (1 at erroneous values)
%%%% Err.filt_paramters to use are in errfilt_pars
% - if stored err.filt_pars needed, then input:
%   errfilt_pars = handles.ts( f_idx ).p.errfilt;
% - if displayed err.filt_pars needed, then input:
%   [ errfilt_pars, err_str ] = gui_collect_errfilt_pars( f_idx, handles);
% ------------------------------------------------------------
% 1. Check err.filt pars for validity:
% ! should happen before calling fcn_errfilt__main
% ------------------------------------------------------------
% 2. Determine valid samples by applying error_filtering methods
% Mark valid and bad elements based on the individual filters
% Caution: consider that (.raw_vals == 0) should become (.vals = 0 )
% -------------------------------------------------------------------------
raw_tslength = size( curr_f_raw_vals, 1 );
% -------------------------------------------------------------------------
% - Default error filtering results == Valid in raw:
raw_valpos   = find( curr_f_raw_vals ==1 );
efilt_res.eftyp         = {'res'};
efilt_res.valmarker     = false( raw_tslength, 1);
efilt_res.valmarker( raw_valpos, 1) = 1;
efilt_res.errmarker     = false( raw_tslength, 1);
efilt_res.badpos_cc_vct = [];
efilt_res.x_ell_vct     = [];
efilt_res.y_ell_vct     = [];
efilt_res.badpos_cc_pst = [];
efilt_res.x_ell_pst     = [];
efilt_res.y_ell_pst     = [];
% -------------------------------------------------------------------------
% - crop time series
ef_mrk=strcmpi( errfilt_pars.used_types,'crop');
if sum( ef_mrk )>0
    %--- Get filtering parameters -----------------------------------------
    crop_lims = errfilt_pars.typ( ef_mrk ).lambda;
    %--- Determine valid elements -----------------------------------------
    valmarker_cr = true( raw_tslength ,1);
    valmarker_cr(...
        curr_f_raw_t < crop_lims(1) |...
        curr_f_raw_t > crop_lims(2) )=0;
    %--- Store results of methods -----------------------------------------
    ef_col = size( efilt_res.valmarker, 2) + 1;
    efilt_res.eftyp{ ef_col }   = 'crop';
    efilt_res.valmarker( :, ef_col ) = false( raw_tslength, 1);
    efilt_res.valmarker( :, ef_col ) = valmarker_cr;
    efilt_res.errmarker( :, ef_col ) = false( raw_tslength, 1);
    efilt_res.errmarker( valmarker_cr==0, ef_col ) = 1;
    %--- Combine elmement validity of different filters -------------------
    efilt_res.valmarker( :, 1 ) = ...
        efilt_res.valmarker( :, 1 ) .* efilt_res.valmarker( :, ef_col );
end
%%%% Raw_valpos elements within crop limits (entire data, if no crop)
curr_valpos = find( efilt_res.valmarker( :, 1 )==1 );
% -------------------------------------------------------------------------
% - PST based spike filter
ef_mrk=strcmpi( errfilt_pars.used_types,'pst');
if sum( ef_mrk )>0 && numel( curr_valpos ) > 4
    % Ellipsoid - number of line elements
    ell_steps = 25;
    %Calculation of lambda:
    pst_lambda = sqrt( 2*log(  raw_tslength  ) );%2*log(datalength)
    %--- Determine valid elements -----------------------------------------
    [ x_ell_pst, y_ell_pst, badpos_cc_pst, valmarker_pst ] =...
        fcn_errfilt_pst_apply(...
        curr_f_raw_veldata( :, 1:3 ),...
        curr_f_raw_vals,...
        curr_valpos, ...
        curr_f_datprops.frq,...
        pst_lambda, ell_steps,...
        errfilt_pars.typ( ef_mrk ).used_comp,...
        errfilt_pars.typ( ef_mrk ).hipass_frq);
    %--- Store results of methods -----------------------------------------
    ef_col = size( efilt_res.valmarker, 2) + 1;
    efilt_res.eftyp{ ef_col }   = 'PST';
    efilt_res.valmarker( :, ef_col ) = false( raw_tslength, 1);
    efilt_res.valmarker( curr_valpos, ef_col ) = logical( valmarker_pst );
    efilt_res.errmarker( :, ef_col ) = false( raw_tslength, 1);
    efilt_res.errmarker( curr_valpos(valmarker_pst==0), ef_col ) = 1;
    efilt_res.badpos_cc_pst{1}    = curr_valpos(badpos_cc_pst{1});
    efilt_res.badpos_cc_pst{2}    = curr_valpos(badpos_cc_pst{2});
    efilt_res.badpos_cc_pst{3}    = curr_valpos(badpos_cc_pst{3});
    efilt_res.x_ell_pst      = x_ell_pst;
    efilt_res.y_ell_pst      = y_ell_pst;
    %--- Combine elmement validity of different filters -------------------
    efilt_res.valmarker( :, 1 ) = ...
        efilt_res.valmarker( :, 1 ) .* efilt_res.valmarker( :, ef_col );
end
% -------------------------------------------------------------------------
% - Correlation filters: Filter out samples where any cor is lower than
ef_mrk = strcmpi( errfilt_pars.used_types,'cormin');
if sum( ef_mrk )>0 && numel( curr_valpos ) > 0
    %--- Get filtering parameters -----------------------------------------
    corthresh = errfilt_pars.typ( ef_mrk ).lambda;
    %--- Determine valid elements -----------------------------------------
    valmarker_mc = true( length( curr_valpos ) ,1);
    valmarker_mc( curr_f_raw_corsnr_xt(curr_valpos ,1) < corthresh )=0;
    %--- Store results of methods -----------------------------------------
    ef_col = size( efilt_res.valmarker, 2) + 1;
    efilt_res.eftyp{ ef_col }   = 'cormin';
    efilt_res.valmarker( :, ef_col ) = false( raw_tslength, 1);
    efilt_res.valmarker( curr_valpos, ef_col ) = valmarker_mc;
    efilt_res.errmarker( :, ef_col ) = false( raw_tslength, 1);
    efilt_res.errmarker( curr_valpos(valmarker_mc==0), ef_col ) = 1;
    %--- Combine elmement validity of different filters -------------------
    efilt_res.valmarker( :, 1 ) = ...
        efilt_res.valmarker( :, 1 ) .* efilt_res.valmarker( :, ef_col );
end
% - Filter out samples where the average cor is lower than
ef_mrk = strcmpi( errfilt_pars.used_types,'coravg');
if sum( ef_mrk )>0 && numel( curr_valpos ) > 0
    %--- Get filtering parameters -----------------------------------------
    corthresh = errfilt_pars.typ( ef_mrk ).lambda;
    %--- Determine valid elements -----------------------------------------
    valmarker_ac = true( length( curr_valpos ) ,1);
    valmarker_ac( curr_f_raw_corsnr_xt(curr_valpos ,2) < corthresh )=0;
    %--- Store results of methods -----------------------------------------
    ef_col = size( efilt_res.valmarker, 2) + 1;
    efilt_res.eftyp{ ef_col }   = 'coravg';
    efilt_res.valmarker( :, ef_col ) = false( raw_tslength, 1);
    efilt_res.valmarker( curr_valpos, ef_col ) = valmarker_ac;
    efilt_res.errmarker( :, ef_col ) = false( raw_tslength, 1);
    efilt_res.errmarker( curr_valpos(valmarker_ac==0), ef_col ) = 1;
    %--- Combine elmement validity of different filters -------------------
    efilt_res.valmarker( :, 1 ) = ...
        efilt_res.valmarker( :, 1 ) .* efilt_res.valmarker( :, ef_col );
end
% -------------------------------------------------------------------------
% - Further methods: DVL
if use_dvl==1
    ell_steps = 25;
    [ ~, efilt_res ] = viper_dvl('dvl_errfilt__main_adv',...
        [], { efilt_res, errfilt_pars, ell_steps, ...
        curr_f_raw_veldata, curr_valpos, curr_f_datprops });
end
% -------------------------------------------------------------------------
% Final combined errmarker:
% -------------------------------------------------------------------------
efilt_res.errmarker( efilt_res.valmarker( :, 1 )==0, 1 ) = 1;






% ------------------------------------------------------------
% Check err.filt parameters for validity - before err.filtering_ts
% ------------------------------------------------------------
function [ err_msg ] = fcn_errfiltpars_check_validity(...
    curr_p_errfilt, curr_dat_props, dvl_ext)
%curr_dat_pos:  [] == not checking hipass_filter against data_frq
%               * otherwise: checking
% Case modifying this: a copy of this fcn is also in viper_dvl!
err_msg = {};
% - Used velcomp - its controlled by the gui, but check here as well
if sum(curr_p_errfilt.used_velcomp(:))==0
    if isempty(err_msg)==1
        err_msg{1,1} = 'Error in error-filter parameters: ';
    end
    err_msg{end+1,1} = '- No velocity component selected';
end
% - Crop raw time-series data
ef_mrk=strcmpi( curr_p_errfilt.used_types,'crop');
if sum( ef_mrk )>0
    %check whether NaNs (in theory, this could come from a parfile)
    if sum( isnan(  curr_p_errfilt.typ( ef_mrk ).lambda  ))>0
        if isempty(err_msg)==1
            err_msg{1,1} = 'Error in error-filter parameters: ';
        end
        err_msg{end+1,1} = '- Invalid crop threshold value';
    end
    %check whether lambda2>lambda1
    if curr_p_errfilt.typ( ef_mrk ).lambda(1) > ...
            curr_p_errfilt.typ( ef_mrk ).lambda(2)
        if isempty(err_msg)==1
            err_msg{1,1} = 'Error in error-filter parameters: ';
        end
        err_msg{end+1,1} = '- Invalid crop threshold values';
    end
    %DO NOT check whether lambda2 < ts_duration
end
% - PST based spike filter
ef_mrk=strcmpi( curr_p_errfilt.used_types,'pst');
if sum( ef_mrk )>0
    %check selected components (not controlled by gui)
    if sum(curr_p_errfilt.typ( ef_mrk ).used_comp(:))==0
        if isempty(err_msg)==1
            err_msg{1,1} = 'Error in PST filter parameters: ';
        end
        err_msg{end+1,1} = '- At least one component has to be selected';
    end
    if isempty( curr_p_errfilt.typ( ef_mrk ).hipass_frq )==0
        %nan or negative
        if isnan(  curr_p_errfilt.typ( ef_mrk ).hipass_frq  )==1 ||...
                curr_p_errfilt.typ( ef_mrk ).hipass_frq <=0
            if isempty(err_msg)==1
                err_msg{1,1} = 'Error in error-filter parameters: ';
            end
            err_msg{end+1,1} =...
                '- Invalid high-pass filter frequency value';
        end
        %file specific check (case dat_props available)
        if isempty( curr_dat_props )==0 && ...
                curr_p_errfilt.typ( ef_mrk ).hipass_frq>curr_dat_props.frq/2
            if isempty(err_msg)==1
                err_msg{1,1} = 'Error in error-filter parameters: ';
            end
            err_msg{end+1,1} =...
                '- High-pass filter frequency not applicable to data';
        end
    end
end
% - Filter out samples where any cor is lower than
ef_mrk=strcmpi( curr_p_errfilt.used_types,'cormin');
if sum( ef_mrk )>0
    if isnan(  curr_p_errfilt.typ( ef_mrk ).lambda  )==1
        if isempty(err_msg)==1
            err_msg{1,1} = 'Error in error-filter parameters: ';
        end
        err_msg{end+1,1} = '- Invalid COR threshold value';
    end
end
ef_mrk=strcmpi( curr_p_errfilt.used_types,'coravg');
if sum( ef_mrk )>0
    if isnan(  curr_p_errfilt.typ( ef_mrk ).lambda  )==1
        if isempty(err_msg)==1
            err_msg{1,1} = 'Error in error-filter parameters: ';
        end
        err_msg{end+1,1} = '- Invalid COR threshold value';
    end
end
if dvl_ext==1%isfield( handles, 'dvl')==1
    [ ~, ~ ] = viper_dvl('dvl_errfiltpars_check_validity',...
        [], { curr_p_errfilt });
end



function [ x_ell, y_ell, badpos_cc, valmarker_pst ] =...
    fcn_errfilt_pst_apply( raw_veldata, raw_vals, curr_valpos, frq,...
    lambda, steps, used_comp, hipass_frq)
% Note: do not output raw_veldata, as it is changed only temporarily
% -------------------------------------------------------------------------
% Allocate variables:
%valmarker - marks the samples that are invalid by 0 ! size == raw_valpos
valmarker_pst  = ones(size( curr_valpos, 1), 1);
valmarker9 = ones(size( curr_valpos, 1), 3, 3);%( samples, comp, pst_plot_id)
badpos_cc = cell(1,3);
x_ell = cell(3,3);
y_ell = cell(3,3);
num_bad3 = nan(1,3);
%Variables required if filtering on
if isempty( hipass_frq )==0%==use highpass
    %For frequency filtering, invalid samples in raw data are interpolated
    pos2int = find( raw_vals==0 );%raw_invalids
    pos_val = find( raw_vals==1 );%raw_valid
    %Interpolate invalids in raw
    raw_veldata(pos2int,:) =...
        interp1( pos_val, raw_veldata( pos_val, :), pos2int );
    %Filtering window
    winsiz = round( frq/hipass_frq );
    win = fcn_blackmanharris( winsiz );
end
vel_data = raw_veldata( curr_valpos, :);
% -------------------------------------------------------------------------
for cmp=1:3
if used_comp( cmp )==1
    % -------------------------------------------------------------------------
    % 0. High-pass filter: excluding low freq. oscillations from filtering
    % HERE using conv instead of filter, because filter introduces a delay!
    if isempty( hipass_frq )==0%==use highpass
        %Removing vel_lowpass (=obtained from convolution)
        raw_veldata(:,cmp) = raw_veldata(:,cmp) -...
            conv( raw_veldata(:,cmp), win'/sum(win),'same');
        %Replace data to filtered:
        vel_data(:,cmp) = raw_veldata( curr_valpos, cmp);
    end
    % -------------------------------------------------------------------------
    % 1a. Remove median from fluctutations:
    vel_fluct = vel_data(:,cmp) - median( vel_data(:,cmp) );
    % -------------------------------------------------------------------------
    % 1b. Prefilter: Remove obvious spikes (Parsheh et al.)
    spike_lim = 1.483*median(abs( vel_fluct ))*lambda;%vel_fluct-median(vel_fluct)
    posbad = find( abs(vel_fluct)>spike_lim );
    posval = find( abs(vel_fluct)<=spike_lim );
    valmarker9( posbad, cmp, 1 ) = 0;
    valmarker9( posbad, cmp, 2 ) = 0;
    valmarker9( posbad, cmp, 3 ) = 0;
    % -------------------------------------------------------------------------
    % 2. Determine valid positions by ellipse fitting - Goring-Nikora
    %subtract median again? not needed as median not sensitiv
    %vel_fluct = vel_data(:,cmp) - median( vel_fluct(posval) );
    curr_num_bad = -1;
    curr_iter = 0;
    while curr_num_bad ~=0 && numel( posval ) > 4
        curr_iter = curr_iter+1;
        %Calc. ellipse parameters from valids and find samples outside:
        [ curr_ell_fi,curr_ell_a2,curr_ell_b2,curr_badpos_cc,curr_fluct ] =...
            fcn_errfilt_pst_ellipse_fit( vel_fluct(posval), lambda );
        % -----------------------------------------------------------------
        % 2b. Rehabilitate false-erroneous ones: (Parsheh et al.)
        % not used, because PDF dependent C-value needed
        % -----------------------------------------------------------------
        % 2c. If second | second last invalid - set first | last -> invalid
        %[ numel(posval) min(curr_badpos_cc{2}) max(curr_badpos_cc{2})]
        for cc=1:3
            if isempty(curr_badpos_cc{cc})~=1
                if curr_badpos_cc{cc}(1)==2
                    curr_badpos_cc{cc}= [1; curr_badpos_cc{1}];
                end
                if curr_badpos_cc{cc}(end)==length(posval)-1
                    curr_badpos_cc{cc}= [curr_badpos_cc{1};length(posval)];
                end
            end
        end
        %badpos in current iteration
        curr_badpos = union(curr_badpos_cc{1},curr_badpos_cc{2});
        curr_badpos = union(curr_badpos,curr_badpos_cc{3});
        curr_num_bad = numel(curr_badpos);
        %Mark bad positions, which were outside ellipse, in this iteration:
        valmarker9( posval( curr_badpos_cc{1} ), cmp, 1 ) = 0;
        valmarker9( posval( curr_badpos_cc{2} ), cmp, 2 ) = 0;
        valmarker9( posval( curr_badpos_cc{3} ), cmp, 3 ) = 0;
        %Combine validity:
        valmarker = prod( valmarker9( :, cmp, :), 3);
        %Redefine valids of the three plots after this iteration:
%         posval4curr = posval;%store previous
        posval = find( valmarker==1 );
        num_bad = length( find( valmarker==0 ) );
    end
%     %Case iteration have not even started:
%     if curr_iter == 0
%         valmarker = prod( valmarker9( :, cmp, :), 3);
%         num_bad = length( find( valmarker==0 ) );
%     end
    num_bad3(cmp) = num_bad;
    % ---------------------------------------------------------------------
    % 3. Determine finally valids
    valmarker_pst = valmarker_pst.*valmarker;
    badpos_cc{cmp} = find( valmarker==0 );
    % ---------------------------------------------------------------------
    % 4. Generate ellipse to plot using final ell_fi & ell_ab2
    [ x_ell(:,cmp), y_ell(:,cmp) ] = fcn_errfilt_pst_ellipse2plot(...
        curr_ell_fi, curr_ell_a2, curr_ell_b2, steps);
else
    badpos_cc{cmp} = [];
    %valmarker - nothing to do
    x_ell{1,cmp} = 0;%zeros(2*steps+1,1);
    x_ell{2,cmp} = 0;%zeros(2*steps+1,1);
    x_ell{3,cmp} = 0;%zeros(2*steps+1,1);
    y_ell{1,cmp} = 0;%zeros(2*steps+1,1);
    y_ell{2,cmp} = 0;%zeros(2*steps+1,1);
    y_ell{3,cmp} = 0;%zeros(2*steps+1,1);
end
end
% fprintf('--- Err.filt pst num_bad [u,v,w]: %d %d %d \n', num_bad3)


function [ x_ell, y_ell ] = fcn_errfilt_pst_ellipse2plot(...
    ell_fi, ell_a2, ell_b2, steps)
%Allocate variables:
x_ell = cell(3,1);
y_ell = cell(3,1);
%%%% Calculate ellipse coordinates for plot:
%Calculate ellipsoids for u-du -> rotation not needed as fi==0
elli_x{1}  =...
    ( -sqrt(ell_a2(1)):(sqrt(ell_a2(1))/steps):sqrt(ell_a2(1)) )';
elli_x{1}(1)=0.999*elli_x{1}(1);
elli_x{1}(end)=0.999*elli_x{1}(end);
elli_y{1} = sqrt(...
    repmat( ell_b2(1) ,[2*steps+1,1]) -...
    repmat( ell_b2(1)/ell_a2(1) ,[2*steps+1,1]) .* (elli_x{1}.^2) );
elli_x{1}=[ elli_x{1}; flipud( elli_x{1} ) ];
elli_y{1}=[ elli_y{1}; flipud(-elli_y{1} ) ];
%Rotate ellipse-coord.sys back to data-coord.sys:
x_ell{1}=elli_x{1};
y_ell{1}=elli_y{1};
%Calculate ellipsoids for % du-ddu -> rotation not needed as fi==0
elli_x{2}  =...
    [ -sqrt( ell_a2(2) ) : (sqrt( ell_a2(2) )/steps) : sqrt( ell_a2(2) )]';
elli_x{2}(1)=0.999*elli_x{2}(1);
elli_x{2}(2*steps+1)=0.999*elli_x{2}(2*steps+1);
elli_y{2} = sqrt(...
    repmat( ell_b2(2) ,[2*steps+1,1]) -...
    repmat( ell_b2(2)/ell_a2(2,1) ,[2*steps+1,1]) .*( elli_x{2} .^2) );
elli_x{2}=[ elli_x{2}; flipud( elli_x{2} ) ];
elli_y{2}=[ elli_y{2}; flipud(-elli_y{2} ) ];
%Rotate ellipse-coord.sys back to data-coord.sys:
x_ell{2}=elli_x{2};
y_ell{2}=elli_y{2};
%Calculate ellipsoids for % u-ddu
elli_x{3}  =...
    [ -sqrt( ell_a2(3) ) : (sqrt( ell_a2(3) )/steps) : sqrt( ell_a2(3) )]';
elli_x{3}(1)=0.999*elli_x{3}(1);
elli_x{3}(2*steps+1)=0.999*elli_x{3}(2*steps+1);
elli_y{3} = sqrt(...
    repmat( ell_b2(3) ,[2*steps+1,1]) -...
    repmat( ell_b2(3)/ell_a2(3) ,[2*steps+1,1]) .*( elli_x{3} .^2) );
elli_x{3}=[ elli_x{3};  flipud( elli_x{3} ) ];
elli_y{3}=[ elli_y{3} ; flipud(-elli_y{3} )];
%Rotate ellipse-coord.sys back to data-coord.sys:
x_ell{3}=elli_x{3}*cos(-ell_fi(3))+elli_y{3}*sin(-ell_fi(3));
y_ell{3}=elli_y{3}*cos(-ell_fi(3))-elli_x{3}*sin(-ell_fi(3));



function [ ell_fi, ell_a2, ell_b2, badpos, fluct ] =...
    fcn_errfilt_pst_ellipse_fit( fluct, lambda )
%Allocate variable
fluct(:,2:3) = zeros( size( fluct, 1), 2 );
ell_fi = zeros(3,1);
% -------------------------------------------------------------------------
% Calculate du' and ddu' series and their variances
fluct(2:end-1,2) = ( fluct(3:end,  1)-fluct(1:end-2,1) )/2;%Goring Nikora
fluct(3:end-2,3) = ( fluct(4:end-1,2)-fluct(2:end-3,2) )/2;%Goring Nikora
% fluct(2:end-1,3) = fluct(3:end,1)-2*fluct(2:end-1,2)+fluct(1:end-2,1);
%Mod: Islam et al
% fluct_bw = zeros( size(fluct,1),1 );
% fluct_fw = zeros( size(fluct,1),1 );
% fluct_bw(2:end-1) = fluct(2:end-1,  1) - fluct(1:end-2,1);
% fluct_fw(2:end-1) = fluct(2:end-1,  1) - fluct(3:end,1);
% fluct_abs = min( abs( fluct_bw ), abs( fluct_fw ) );
% idx_bw = find( fluct_abs == abs( fluct_bw ));
% idx_fw = find( fluct_abs == abs( fluct_fw ));
% fluct(idx_bw,2) = fluct_bw(idx_bw);
% fluct(idx_fw,2) = fluct_fw(idx_fw);
% fluct_bw = zeros( size(fluct,1),1 );
% fluct_fw = zeros( size(fluct,1),1 );
% fluct_bw(3:end-2) = fluct(3:end-2,  1) - fluct(2:end-3,1);
% fluct_fw(3:end-2) = fluct(3:end-2,  1) - fluct(4:end-1,1);
% fluct_abs = min( abs( fluct_bw ), abs( fluct_fw ) );
% idx_bw = find( fluct_abs == abs( fluct_bw ));
% idx_fw = find( fluct_abs == abs( fluct_fw ));
% fluct(idx_bw,3) = fluct_bw(idx_bw);
% fluct(idx_fw,3) = fluct_fw(idx_fw);
% -------------------------------------------------------------------------
% Calculate median of absolute deviations (MAD) instead of variances
stds(1) = 1.483*median(abs( fluct(:,1)-median(fluct(:,1)) ));%mod: Wahl
stds(2) = 1.483*median(abs( fluct(:,2)-median(fluct(:,2)) ));%mod: Wahl
stds(3) = 1.483*median(abs( fluct(:,3)-median(fluct(:,3)) ));%mod: Wahl
% stds(1) = std( fluct(:,1), 1);%original Goring Nikora
% stds(2) = std( fluct(:,2), 1);%original Goring Nikora
% stds(3) = std( fluct(:,3), 1);%original Goring Nikora
% -------------------------------------------------------------------------
% Calculate rotation angle of principal axes:
ell_fi(1) = 0;%0 due to symmetry
ell_fi(2) = 0;%0 due to symmetry
ell_fi(3) = atan(...
    sum( fluct(:,1).*fluct(:,3) )/sum( fluct(:,1).^2 ) );
% -------------------------------------------------------------------------
% Calculate major & minor axis lengths using std_u, std_du, std_ddu
ell_a2(1,1) = ( lambda*stds(1) )^2;%for {1}-{2} % u-du
ell_b2(1,1) = ( lambda*stds(2) )^2;%for {1}-{2} % u-du
ell_a2(2,1) = ( lambda*stds(2) )^2;%for {2}-{3} % du-ddu
ell_b2(2,1) = ( lambda*stds(3) )^2;%for {2}-{3} % du-ddu
%for {3}-{1} % u-ddu
ell_a2(3,1) = ...
    ( (lambda*stds(3))^2 * sin(ell_fi(3))^2 - ...
    (lambda*stds(1))^2 * cos(ell_fi(3))^2 )/...
    ( sin(ell_fi(3))^4 - cos(ell_fi(3))^4 );%abs(
ell_b2(3,1) = ...
    ( (lambda*stds(3))^2 * cos(ell_fi(3))^2 - ...
    (lambda*stds(1))^2 * sin(ell_fi(3))^2 )/...
    ( cos(ell_fi(3))^4 - sin(ell_fi(3))^4 );%abs(
%If a<b, swap ellipsoid axes:
% for plt=1:3
% if ell_a2(plt)<ell_b2(plt)
%     tmp_ab2 = ell_a2(plt);
%     ell_a2(plt) = ell_b2(plt);
%     ell_b2(plt) = tmp_ab2;
% end
% end
% -------------------------------------------------------------------------
% Determine positions outside ellipse:
badpos{1}=find( ...
    ( fluct(:,1).^2/ell_a2(1) + fluct(:,2).^2/ell_b2(1) )>1);
badpos{2}=find( ...
    ( fluct(:,2).^2/ell_a2(2) + fluct(:,3).^2/ell_b2(2) )>1);
% Here consider angle -> rotate data to ellipsoid coordinate-system:
if ell_a2(3,1)>0 && ell_b2(3,1)>0
    badpos{3}=find(...
        ((fluct(:,1)*cos(ell_fi(3))+fluct(:,3)*sin(ell_fi(3))).^2/...
           ell_a2(3)+...
         (fluct(:,3)*cos(ell_fi(3))-fluct(:,1)*sin(ell_fi(3))).^2/...
           ell_b2(3)) > 1 );
else
    badpos{3}=[];
end





















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 
function chapter_14_calculating_statistical_features
% 
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -------------------------------------------------------------------------
% Generate quality plots of valid samples. Case no filtering: use all
% -------------------------------------------------------------------------
function [ temp_corsnr,temp_boxp ] = fcn_calc_corsnrdist_vals_adv_sokoray(...
    curr_f_vals, curr_f_rawcorsnr_xt )
%%% Num of samples:
raw_tslength = size(  curr_f_vals, 1 );
num_valid = sum( curr_f_vals  );
%%% Valid positions:
valpos = find(  curr_f_vals==1  );
% -------------------------------------------------------------------------
% CORRELATION DISTRIBUTION: c; SNR DISTRIBUTION: s
% -------------------------------------------------------------------------
% For the distribution: data is classified to bins (steps in quality)
numbins = 50;
%Calculate bin edges and midpoints
bin_edge_c = ( 0:( 100/ numbins ):100 )';% bin edges: numbins+1
bin_mid_c = ( 0+(100/numbins)/2 : (100/numbins): 100)';%bin mids: numbins
bin_edge_s  = ( 0:( 25/ numbins ):25 )';% bin edges: numbins+1
bin_mid_s = ( 0+(25/numbins)/2 : (25/numbins): 25 )';%bin mids: numbins
% -------------------------------------------------------------------------
% Number of elements within bins: Mark elements that are within the 
% bin-edges - each columns is a valmarker for one bin
marker_bins_mc = zeros( raw_tslength, numbins+1 );
marker_bins_ac = zeros( raw_tslength, numbins+1 );
marker_bins_ms = zeros( raw_tslength, numbins+1 );
marker_bins_as = zeros( raw_tslength, numbins+1 );
for fls=1:numbins+1
    %case minimal-cor of the 4 beams:
    pos_val =...
        curr_f_rawcorsnr_xt(valpos ,1) <=bin_edge_c(fls);
    marker_bins_mc( valpos( pos_val) ,fls ) = 1;
    %case average-cor of the 4 beams:
    pos_val =...
        curr_f_rawcorsnr_xt(valpos ,2) <=bin_edge_c(fls);
    marker_bins_ac( valpos( pos_val) ,fls ) = 1;
    %case minimal-snr of the 4 beams:
    pos_val =...
        curr_f_rawcorsnr_xt(valpos ,3) <=bin_edge_s(fls);
    marker_bins_ms( valpos( pos_val) ,fls ) = 1;
    %case average-snr of the 4 beams:
    pos_val =...
        curr_f_rawcorsnr_xt(valpos ,4) <=bin_edge_s(fls);
    marker_bins_as( valpos( pos_val) ,fls ) = 1;
end
% -------------------------------------------------------------------------
%%%Percentage of samples in bins
%(Percent calculated based on number of existing samples in raw)
% COR Cumulative-DISTRIBUTION:
pct_in_bins_c_mc = 100*sum( marker_bins_mc ,1 )./ num_valid;
pct_in_bins_c_ac = 100*sum( marker_bins_ac ,1 )./ num_valid;
% COR DISTRIBUTION:
pct_in_bins_mc = [ 0,diff( pct_in_bins_c_mc )];
pct_in_bins_ac = [ 0,diff( pct_in_bins_c_ac )];
% SNR Cumulative-DISTRIBUTION:
pct_in_bins_c_ms = 100*sum( marker_bins_ms ,1 )./ num_valid;
pct_in_bins_c_as = 100*sum( marker_bins_as ,1 )./ num_valid;
% SNR DISTRIBUTION:
pct_in_bins_ms = [ 0,diff( pct_in_bins_c_ms )];
pct_in_bins_as = [ 0,diff( pct_in_bins_c_as )];
%%%% SAVE
temp_corsnr.corsnr_dist = cell(4,4);
%Correlation distribution:
temp_corsnr.corsnr_dist(1,1:4) =...
    {bin_edge_c, bin_mid_c, pct_in_bins_mc, pct_in_bins_ac};
%Correlation cumulative:
temp_corsnr.corsnr_dist(2,1:4) =...
    {bin_edge_c, bin_mid_c, 100-pct_in_bins_c_mc, 100-pct_in_bins_c_ac};
%SNR distribution:
temp_corsnr.corsnr_dist(3,1:4) =...
    {bin_edge_s, bin_mid_s, pct_in_bins_ms, pct_in_bins_as};
%SNR cumulative:
temp_corsnr.corsnr_dist(4,1:4) =...
    {bin_edge_s, bin_mid_s, 100-pct_in_bins_c_ms, 100-pct_in_bins_c_as};
% -------------------------------------------------------------------------
% Calculate Cor-Distribution boxplot (cd_boxp) values:
% -------------------------------------------------------------------------
cd_boxp = cell(1,5);
%min:
cd_boxp{1} = min( curr_f_rawcorsnr_xt(valpos ,1) );
%median:
cd_boxp{3} = median( curr_f_rawcorsnr_xt(valpos ,1) );
%max:
cd_boxp{5} = max( curr_f_rawcorsnr_xt(valpos ,1) );
%%% 25% quartil:
ind25_1 = find( pct_in_bins_c_mc <= 25 ,1,'last');
ind25_2 = find( pct_in_bins_c_mc >= 25 ,1,'first');
if pct_in_bins_c_mc(ind25_1)==pct_in_bins_c_mc(ind25_2)
    cd_boxp{2} = bin_mid_c(ind25_1);
else
    cd_boxp{2} = ...
        bin_edge_c(ind25_1) + (25-pct_in_bins_c_mc(ind25_1))/...
        (pct_in_bins_c_mc(ind25_2)-pct_in_bins_c_mc(ind25_1)) *...
        (bin_edge_c(ind25_2)-bin_edge_c(ind25_1));
end
%75% quartil:
ind75_1 = find( pct_in_bins_c_mc <= 75 ,1,'last');
ind75_2 = find( pct_in_bins_c_mc >= 75 ,1,'first');
if pct_in_bins_c_mc(ind75_1)==pct_in_bins_c_mc(ind75_2)
    cd_boxp{4} = bin_mid_c(ind75_1);
else
    cd_boxp{4} = ...
        bin_edge_c(ind75_1) + (75-pct_in_bins_c_mc(ind75_1))/...
        (pct_in_bins_c_mc(ind75_2)-pct_in_bins_c_mc(ind75_1)) * ...
        (bin_edge_c(ind75_2)-bin_edge_c(ind75_1));
end
%Check! If empty -> replace by 0
if isempty(cd_boxp{1})==1 || isempty(cd_boxp{2})==1 ||...
        isempty(cd_boxp{3})==1 || isempty(cd_boxp{4})==1 ||...
        isempty(cd_boxp{5})==1
    for c=1:5
        cd_boxp{c}=0;
    end
end
temp_boxp = cell2mat( cd_boxp );



% -------------------------------------------------------------------------
% Generate quality plots of valid samples. Case no filtering: use all
% -------------------------------------------------------------------------
function [ p_per_bin, vel_of_bin ] = fcn_calc_histograms_of_vals(...
    curr_f_raw_veldata, curr_f_vals, curr_f_used_velcomp )
%%% Valid positions - are determined by err.filtering scripts
valpos = curr_f_vals==1;
numval = sum(valpos);
num_bins = 25;
%Initialize variables: nans
p_per_bin = nan( num_bins,3 );
vel_of_bin= nan( num_bins,3 );
% ---------------------------------------------------------------------
% Calculate histogram
for cmp_idx=1:3
if curr_f_used_velcomp( cmp_idx )==1
    [ p_per_bin(:,cmp_idx) , vel_of_bin(:,cmp_idx) ] = hist(...
        curr_f_raw_veldata( valpos ,cmp_idx), num_bins );
    p_per_bin(:,cmp_idx) = p_per_bin(:,cmp_idx)/numval*100;
end
end



% -------------------------------------------------------------------------
% Generates error filtered fluctuations and valid data based statistics:
% -------------------------------------------------------------------------
function [ new_f_seri_data, new_f_seri_enable, new_f_meancrtn, new_f_usedseg,...
    temp_f_seri_x, temp_f_valu_data, temp_f_valu_enable, temp_f_valuupdate ] =...
    fcn_calc_statistical_func_n_valu( ...
    curr_f_raw_veldata, curr_f_vals, curr_f_p, curr_f_raw_t,...
    curr_f_raw_corsnr_xt,  default_ts, siz_valu_data )
% INFO: Raw time-series length:
% raw_tslength = size( curr_f_vals, 1 );
% INFO: Valid and erroneous positions - are determined by err.filtering scripts
% valpos = find( curr_f_vals==1  );
% errpos = find( curr_f_vals==0  );
% Used segment: (without initial and trailing erroneous samples)
new_f_usedseg = false(size( curr_f_vals ));
new_f_usedseg(...
    find( curr_f_vals==1,1,'first') : find( curr_f_vals==1,1,'last') ) =...
    true;
seri_tslength = sum( new_f_usedseg );
% INFO: Valid and erroneous positions within segment of seri
% valpos_in_seri = find( curr_f_vals( used4seri )==1  );
% errpos_in_seri = find( curr_f_vals( used4seri )==1  );
%Time axes for seri (time, Autocorrelation time lag):
temp_f_seri_x.t = curr_f_raw_t( new_f_usedseg );
temp_f_seri_x.lag  = (0:seri_tslength-1)';
temp_f_seri_x.tlag = temp_f_seri_x.lag/curr_f_p.dat_props.frq;
temp_f_seri_x.f = [0 0];
%Initialize valu-variables as fullsize, in order to store data to right pos
temp_f_valu_data   = nan( siz_valu_data );%store valu data to right pos
temp_f_valu_enable = false( siz_valu_data );
temp_f_valuupdate = false( siz_valu_data );%update marker
%Initialize seri-variables: reset all to defaults and recalc again
new_f_seri_data   = default_ts.seri_data; %cell(8,3);
new_f_seri_enable = default_ts.seri_enable;
% ---------------------------------------------------------------------
% Case valid samples - No sense to calculate anything if  less than 2
% ---------------------------------------------------------------------
if sum( curr_f_vals )>1
    %Initialize variables
    new_f_seri_data{1} = nan( seri_tslength, 3 );
    new_f_seri_data{3} = nan( seri_tslength, 3 );
    new_f_seri_data{4} = nan( seri_tslength, 3 );
    %PSD length = (nfft/2 + 1).
    %NFFT used here is length(ts), if length(ts) odd -> ignore last sample
    new_f_seri_data{8} = nan( floor(seri_tslength/2)+1, 3 );
    %mean correction
    new_f_meancrtn = nan(1,3);
    % ---------------------------------------------------------------------
    % Calculate mean, variance and fluctuations-series: (valpos needed)
    % ---------------------------------------------------------------------
    for cmp_idx=1:3
    if curr_f_p.errfilt.used_velcomp( cmp_idx )==1
        %%% Features using raw data:
        % Calculate mean values:
        temp_f_valu_data( cmp_idx ) =...
            mean( curr_f_raw_veldata( curr_f_vals ,cmp_idx) );
        % Calculate variance values:
        temp_f_valu_data( cmp_idx +3 ) =...
            var( curr_f_raw_veldata( curr_f_vals ,cmp_idx ),1 );
        % Calculate std values:
        temp_f_valu_data( cmp_idx +6 ) =...
            sqrt( temp_f_valu_data( cmp_idx +3 ) );
        % Calculate 4th central moment values:
        temp_f_valu_data( cmp_idx +27 ) =...
            sum( ( curr_f_raw_veldata( curr_f_vals ,cmp_idx ) -...
                   temp_f_valu_data( cmp_idx ) ).^4 ) /...
            sum( curr_f_vals );
        % Calculate turbulence intensity numerator (I== u_std/|U|)
        temp_f_valu_data( cmp_idx +42 ) = temp_f_valu_data( cmp_idx +6 );
        % Generate Fluctuation time-series to use (nans end up in nans)
        new_f_seri_data{1}(:, cmp_idx ) =...
            curr_f_raw_veldata( new_f_usedseg, cmp_idx ) -...
            temp_f_valu_data( cmp_idx );
        %%% Features involving interpolated values of seri!
        if isempty( find( curr_f_vals( new_f_usedseg )==0, 1  ) )~=1
            % Replace bad data by linear interpolation
            % No extrapolation - those samples -> nan
            new_f_seri_data{1}( curr_f_vals( new_f_usedseg )==0 , cmp_idx ) =...
                interp1(...
                temp_f_seri_x.t( curr_f_vals( new_f_usedseg )==1 ),...
                new_f_seri_data{1}( curr_f_vals( new_f_usedseg )==1, cmp_idx ),...
                temp_f_seri_x.t( curr_f_vals( new_f_usedseg )==0 ),...
                'linear');%outside valid -> nan
            % Correction of mean - needed due to interpolated values:
            new_f_meancrtn( cmp_idx ) = ...
                mean( new_f_seri_data{1}(:, cmp_idx ));
            %Set mean of interpolated_fluct to 0:
            new_f_seri_data{1}(:, cmp_idx )=...
                new_f_seri_data{1}(:, cmp_idx ) - ...
                new_f_meancrtn( cmp_idx );%
        else
            new_f_meancrtn( cmp_idx ) = 0;
        end
        % Generating  auto-correlation function (Original Autor: P. Heneka)
        % using fft(x,2*length) you get what xcorr gives as follows:
        % For autocorr function -> normalize by std^2!!
        temp_seri = ...
            1/seri_tslength * 1/temp_f_valu_data( cmp_idx +3 ) *...
            ifft(fft(new_f_seri_data{1}(:, cmp_idx ), 2*seri_tslength).*...
            conj(fft(new_f_seri_data{1}(:, cmp_idx ), 2*seri_tslength)));
        new_f_seri_data{3}(:, cmp_idx ) = temp_seri(1:seri_tslength);
        % Generate cumulative mean:
        new_f_seri_data{4}(:, cmp_idx ) =...
            cumsum(new_f_seri_data{1}(:, cmp_idx ))./(1:seri_tslength)'+...
            temp_f_valu_data( cmp_idx );
        %%% Calculate power spectrum density estimation using Welch
        %WARNING: pwelch outputs 0Hz and corresponding psd value -> remove
        %NFFT==num of points of DFT. - USE NFFT that is even!
        %VIPER: ts_length %%% Alternative: 2^nextpow2(ts_length)
        %if your NFFT~=ts_length -> adjust initial seri size before loop
        if rem( seri_tslength , 2 )==1
            NFFT = seri_tslength-1;
        else
            NFFT = seri_tslength;
        end
        %Number of segments:
        %more segments reduce variance of spectrum but also frq. resolution
        nseg = 4;%4/5
        %Using blackmanharris window
        %Overlap ratio of segments:
        roverlap = 0.5;%use 0.5 as taper of blackmanharris==0.5
        %Segment length:
        if rem( seri_tslength , 2 )==1
            seg_len = floor( (seri_tslength-1) / (nseg + roverlap - roverlap*nseg) );
        else
            seg_len = floor( (seri_tslength) / (nseg + roverlap - roverlap*nseg) );
        end
        %Overlap input - CAUTION:
        % - Matlab TSProcessing Toolbox needs window length!
        % - external function: ratio
        % But both accept [] and have 0.5 as default!
        overlap = [];
%         overlap = floor(roverlap*seg_len);%Matlab
%         overlap = roverlap;%External
        [ new_f_seri_data{8}(:, cmp_idx ), temp_f_seri_x.f ] =...
            pwelch( detrend( new_f_seri_data{1}(:, cmp_idx ) ),...
            fcn_blackmanharris( seg_len ), overlap,...
            NFFT, curr_f_p.dat_props.frq );%output = col. vectors
        % Calculate TEI based on Spectra (Sokoray Thesis)
%         curr_dft_raw = fft(new_f_seri_data{1}(:,cmp_idx),ts_length);
%         curr_psd = 2*abs(curr_dft_raw(2:ts_length/2+1)).^2/...
%            (curr_f_p.dat_props.frq*ts_length);
%         temp_f_valu_data( cmp_idx +24 ) = pi/2*...
%             curr_psd(1)/...
%             temp_f_valu_data( cmp_idx +3 );
        %using  seri{8} instead of curr_psd(1):
%         temp_f_valu_data( cmp_idx +24 ) = pi/2*...
%             new_f_seri_data{8}(1, cmp_idx )/...
%             temp_f_valu_data( cmp_idx +3 );
        % Calculate TEI based on Autocorr_fcn (old_tool and also Heneka)
        zero_pos = find(  new_f_seri_data{3}(:, cmp_idx )<=0  ,1);
        temp_f_valu_data( cmp_idx +24 ) = ...
            trapz( temp_f_seri_x.tlag(1:zero_pos),...
            new_f_seri_data{3}(1:zero_pos, cmp_idx )  );
        % Calculate s(u_mean) and e(u_mean) (Sokoray Thesis)*2 for 95%conf
        temp_f_valu_data( cmp_idx +30 ) = ...
            2*sqrt( 2 * temp_f_valu_data( cmp_idx +24 )/...
                (seri_tslength/curr_f_p.dat_props.frq)*...
                temp_f_valu_data( cmp_idx +3 ) );
        temp_f_valu_data( cmp_idx +33 ) = ...
            temp_f_valu_data( cmp_idx +30 )/...
            abs(temp_f_valu_data( cmp_idx ))*100;
        % Calculate s(u_var) and e(u_var) (Sokoray Thesis)*2 for 95%conf
        temp_f_valu_data( cmp_idx +36 ) = ...
            2*sqrt( 2 * temp_f_valu_data( cmp_idx +24 )/...
                (seri_tslength/curr_f_p.dat_props.frq)*...
                (temp_f_valu_data( cmp_idx +27 ) - ...
                 temp_f_valu_data( cmp_idx +3 )^2) );
        temp_f_valu_data( cmp_idx +39 ) = ...
            temp_f_valu_data( cmp_idx +36 )/...
            temp_f_valu_data( cmp_idx +3 )*100;
    else
        temp_f_valu_data( cmp_idx ) = nan;
        temp_f_valu_data( cmp_idx +3 ) = nan;
        temp_f_valu_data( cmp_idx +6 ) = nan;
        temp_f_valu_data( cmp_idx +24 ) = nan;
        temp_f_valu_data( cmp_idx +27 ) = nan;
        temp_f_valu_data( cmp_idx +30 ) = nan;
        temp_f_valu_data( cmp_idx +33 ) = nan;
        temp_f_valu_data( cmp_idx +36 ) = nan;
        temp_f_valu_data( cmp_idx +39 ) = nan;
        temp_f_valu_data( cmp_idx +42 ) = nan;
        new_f_meancrtn( cmp_idx ) = 0;
    end
    end
    cmp_idx = find(curr_f_p.errfilt.used_velcomp==1)';
    % ---------------------------------------------------------------------
    % Correcting pwelch results: remove 0 Hz and corresponding PSD-value
    % ---------------------------------------------------------------------
    new_f_seri_data{8}(1, : ) = [];
    temp_f_seri_x.f(1) = [];
%     %Remove rows unneded for psd:
%     psd_length = length( curr_x_var.f );
%     if size( new_f_seri_data{8}, 1) > psd_length
%         new_f_seri_data{8}( psd_length+1:end, : )= [];
%     end
    % ---------------------------------------------------------------------
    % Others including all components:
    % ---------------------------------------------------------------------
    %vel_magn
    temp_f_valu_data( 10 )=...
        sqrt(   sum( temp_f_valu_data( 0+cmp_idx ).^2 ));
    %tke
    temp_f_valu_data( 11 )=...
        sum( temp_f_valu_data( 3+cmp_idx ) )/2;
    %vel_rms
    temp_f_valu_data( 12 )=...
        sum( temp_f_valu_data( 6+cmp_idx ) )/3;
    %covariance
    cov_mat = cov( curr_f_raw_veldata( curr_f_vals==1 ,:), 1);
    % u*v:
    if curr_f_p.errfilt.used_velcomp(1)==1 &&...
            curr_f_p.errfilt.used_velcomp(2)==1
        temp_f_valu_data( 13 )= cov_mat(1,2);
    else
        temp_f_valu_data( 13 )= NaN;
    end
    % u*w:
    if curr_f_p.errfilt.used_velcomp(1)==1 &&...
            curr_f_p.errfilt.used_velcomp(3)==1
        temp_f_valu_data( 14 )= cov_mat(1,3);
    else
        temp_f_valu_data( 14 )= NaN;
    end
    % v*w:
    if curr_f_p.errfilt.used_velcomp(2)==1 &&...
            curr_f_p.errfilt.used_velcomp(3)==1
        temp_f_valu_data( 15 )= cov_mat(2,3);
    else
        temp_f_valu_data( 15 )= NaN;
    end
    %valid ratio%
    temp_f_valu_data( 16 ) =...
        sum( curr_f_vals )/seri_tslength*100;
    %avg-cor, avg-snr:
    if strcmp( curr_f_p.dat_props.type, 'usr' )==1
        temp_f_valu_data( 17 ) = NaN;
        temp_f_valu_data( 18 ) = NaN;
    else
        temp_f_valu_data( 17 ) = ...
            mean( curr_f_raw_corsnr_xt( curr_f_vals==1 ,2) );
        temp_f_valu_data( 18 ) = ...
            mean( curr_f_raw_corsnr_xt( curr_f_vals==1 ,4) );
    end
    %valid duration%
    temp_f_valu_data( 24 ) = seri_tslength/curr_f_p.dat_props.frq;
    % Turbulence intensity components' denominators (I== u_std/|U|)
    temp_f_valu_data( 43 ) = temp_f_valu_data( 43 )/temp_f_valu_data( 10 );
    temp_f_valu_data( 44 ) = temp_f_valu_data( 44 )/temp_f_valu_data( 10 );
    temp_f_valu_data( 45 ) = temp_f_valu_data( 45 )/temp_f_valu_data( 10 );
    %turbulence intensity (whole vector)= vel_rms/vel_magn:
    temp_f_valu_data( 46 ) = temp_f_valu_data( 12 )/temp_f_valu_data( 10 );
    % ---------------------------------------------------------------------
    % Enabling
    % ---------------------------------------------------------------------
    %Enabling seri:
    new_f_seri_enable(1)=1;%fluctuation
    new_f_seri_enable(3)=1;%autocorr
    new_f_seri_enable(4)=1;%cumm_mean
    new_f_seri_enable(8)=1;%spectrum
    %Enabling statistics:
    temp_f_valu_enable( 1:3 )=1;%mean
    temp_f_valu_enable( 4:6 )=1;%var
    temp_f_valu_enable( 7:9 )=1;%std
    temp_f_valu_enable( 10:12 )=1;%3D ertekek
    temp_f_valu_enable( 13:15 )=1;%cov
    temp_f_valu_enable( 16:18 )=1;%quality
    temp_f_valu_enable( 24 )=1;%quality
    temp_f_valu_enable( 25:27 )=1;%TEI
    temp_f_valu_enable( 28:30 )=1;%Kurtosis
    temp_f_valu_enable( 31:42 )=1;%duration err
    temp_f_valu_enable( 43:46 )=1;%turbulence intensity
    %Set as updated:
    temp_f_valuupdate( 1:3 )=1;%mean
    temp_f_valuupdate( 4:6 )=1;%var
    temp_f_valuupdate( 7:9 )=1;%std
    temp_f_valuupdate( 10:12 )=1;%3D ertekek
    temp_f_valuupdate( 13:15 )=1;%cov
    temp_f_valuupdate( 16:18 )=1;%quality
    temp_f_valuupdate( 24 )=1;%quality
    temp_f_valuupdate( 25:27 )=1;%TEI
    temp_f_valuupdate( 28:30 )=1;%Kurtosis
    temp_f_valuupdate( 31:42 )=1;%duration err
    temp_f_valuupdate( 43:46 )=1;%turbulence intensity
    
% ---------------------------------------------------------------------
% Case no valid samples (less than 2):
% ---------------------------------------------------------------------
else
    %Enabling seri:
    new_f_seri_enable = default_ts.seri_enable;
    %defaults valu - except for quality values!!!
    temp_f_valu_data( 1:18)  = default_ts.valu_data( 1:18);
    temp_f_valu_data( 24:46) = default_ts.valu_data(24:46);
    %valu enable off
    temp_f_valu_enable( 1:18 )  = default_ts.valu_enable( 1:18 );
    temp_f_valu_enable( 24:46 ) = default_ts.valu_enable( 24:46 );
    % ---------------------------------------------------------------------
    % Quality values:
    % ---------------------------------------------------------------------
    % cor-dist boxplot values calculated by calc_qal_plot
    temp_f_valu_data( 16 ) = sum( curr_f_vals )/seri_tslength*100;%valid%
    temp_f_valu_data( 17 ) = ...
        mean( curr_f_raw_corsnr_xt( curr_f_vals==1 ,2) );
    temp_f_valu_data( 18 ) = ...
        mean( curr_f_raw_corsnr_xt( curr_f_vals==1 ,4) );
    %set as updated
    temp_f_valuupdate( 1:18 )  = 1;
    temp_f_valuupdate( 24:46 ) = 1;
    %internal data:
    new_f_meancrtn = [0 0 0];
end





















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 
function chapter_15_coordinates_calculation
% 
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ curr_f_raw_veldata, curr_f_seri_data, curr_f_valu_data,...
    curr_f_intl, curr_f_p, curr_f_efiltres, alpha ] = fcn_calc_rotate_uv(...
    only_raw, dir_old, dir_new,...
    curr_f_raw_veldata, curr_f_seri_data, curr_f_valu_data,...
    curr_f_intl, curr_f_p, curr_f_efiltres, curr_f_seri_enable, field_props )
%only_raw [ 0 | 1 ]: do rotation only for .raw_veldata and .cmp_probe_uvw
%      * needed at adding files
% -------------------------------------------------------------------------
%I. Create pos vector for rotation
%1.  old direction: INPUT
if strcmp( dir_old, '+x' )==1
    deg_old = 0;
elseif strcmp( dir_old, '-x' )==1
    deg_old = 180;
elseif strcmp( dir_old, '+y' )==1
    deg_old = 90;
elseif strcmp( dir_old, '-y' )==1
    deg_old = 270;
end
%2. get new direction
if strcmp( dir_new, '+x' )==1
    deg_new = 0;
elseif strcmp( dir_new, '-x' )==1
    deg_new = 180;
elseif strcmp( dir_new, '+y' )==1
    deg_new = 90;
elseif strcmp( dir_new, '-y' )==1
    deg_new = 270;
end
%3. calc angle of needed rotation AND calc rotation matrix
alpha = deg_new - deg_old;
M_rot = [ cosd(alpha) -sind(alpha) 0;...
          sind(alpha)  cosd(alpha) 0;...
          0 0 1];
%4. calc rotated component order:
comp_trans = transpose( M_rot*[ 1 2 3]');
%5. calc rotated column positions and signs for u-v-w components:
pos_new = abs(comp_trans);
sgn_new = sign(comp_trans);
% -------------------------------------------------------------------------
%II. rotate raw data
% A. Rearrange u-v-w columns of affected variables (==rotate) of current ts
% B. adjust sign BUT not for every variable!!!
%1. .raw_veldata: only first 3 comlumns uvw - rearrange & adjust signs
curr_f_raw_veldata(:,1:3) = ...
    curr_f_raw_veldata(:,pos_new).*...
    repmat( sgn_new, [size( curr_f_raw_veldata(:,1:3),1) ,1] );
%2. current components in probe u-v-w and componentnp_attrib
curr_f_intl.cmp_probe_uvw =...
    curr_f_intl.cmp_probe_uvw(pos_new);
curr_f_intl.cmp_probe_np = ...
    curr_f_intl.cmp_probe_np(pos_new);
% -------------------------------------------------------------------------
%III. rotate data that is calculated from raw data (not at addfile)
% A. Rearrange u-v-w columns of affected variables (==rotate) of current ts
% B. adjust sign BUT not for every variable!!!
if only_raw == 0
    %0. modify cordsys data: .p and .intl
    curr_f_p.coordsys.dir_probe_u   = dir_new;
    curr_f_intl.coordsys.dir_probe_u = dir_new;
    %1. mean_correctn: rearrange comlumns:
    curr_f_intl.mean_correctn = ...
        curr_f_intl.mean_correctn(1,pos_new).*sgn_new;
    %2. .seri_data: DONE only if .seri is enabled (not at addfile)
    for c_i=1:length( curr_f_seri_data )%seri_data == cell
        if curr_f_seri_enable( c_i )==1
            %rearrange comlumns:
            curr_f_seri_data{c_i} = curr_f_seri_data{c_i}(:,pos_new);
            %adjust signs:
            if field_props.seri_adjrotsign( c_i )==1
                curr_f_seri_data{c_i} =...
                    curr_f_seri_data{c_i}.* repmat( sgn_new,...
                    [ size( curr_f_seri_data{c_i},1) ,1 ] );
            end
        end
    end
    %3. .valu_data: only positions where .valu_uvw is nonzero...
    pos_valu_uvw = find( field_props.valu_uvw~=0 );
    old_valu_uvw = reshape( curr_f_valu_data( pos_valu_uvw ),...
        [3,length(pos_valu_uvw)/3]);
    %FOR MOST: rearrange comlumns:
    curr_f_valu_data( pos_valu_uvw ) = old_valu_uvw(pos_new,:);
    %FOR MOST: adjust signs: only for mean values:
    pos_valu_uvw_sgn = find( field_props.valu_uvw~=0 &...
        field_props.valu_adjrotsign==1);
    old_valu_uvw_sgn = reshape( curr_f_valu_data( pos_valu_uvw_sgn ),...
        [3,length(pos_valu_uvw_sgn)/3]);
    curr_f_valu_data( pos_valu_uvw_sgn ) = old_valu_uvw_sgn .*...
        repmat( sgn_new', [1,size( old_valu_uvw_sgn, 2)] );
    %FOR cov_uv cov_uw cov_vw: rearrange:
    if rem( alpha, 180 )~=0%90,-90,270,-270
        curr_f_valu_data( [13,14,15] ) = curr_f_valu_data( [13,15,14] );
    %else: do nothing
    end
    %FOR cov_uv cov_uw cov_vw: signs
    curr_f_valu_data( [13,14,15] ) = curr_f_valu_data( [13,14,15] ).*...
        [sgn_new(1)*sgn_new(2) sgn_new(1) sgn_new(2)];
    %4. Errfilt pars and res!
    % errfilt pars: velcomp & used_comp (if existing)
    % errfilt res: badpos and x_ell y_ell
    curr_f_p.errfilt.used_velcomp =...
        curr_f_p.errfilt.used_velcomp(pos_new);
    if sum(strcmpi(curr_f_p.errfilt.used_types,'pst'))>0
        ef_mrk=strcmpi(curr_f_p.errfilt.used_types,'pst');
        curr_f_p.errfilt.typ( ef_mrk ).used_comp = ...
            curr_f_p.errfilt.typ( ef_mrk ).used_comp(pos_new);
        curr_f_efiltres.badpos_cc_pst =...
            curr_f_efiltres.badpos_cc_pst(pos_new);
        curr_f_efiltres.x_ell_pst = ...
            curr_f_efiltres.x_ell_pst(:,pos_new);
        curr_f_efiltres.y_ell_pst = ...
            curr_f_efiltres.y_ell_pst(:,pos_new);
    end
else%only_raw==1
    curr_f_efiltres =[];
end




% ------------------------------------------------------------
% Retrieve coordinates from string
% ------------------------------------------------------------
function [ cds ] = fcn_cds_read_from_str_sokoray(pos_cds, cds_fmt, str)
% str_cds is coordinate format string
% pos_cds are the positions of the coordinate markers (e.g. %x%) in str_cds
% replaces automatically ,-s to .-s
%Check:
if isempty(pos_cds) || isempty(cds_fmt) || isempty(str) ||...
        length(str)<pos_cds(3)
    cds = [NaN NaN NaN];
    return
end
try
    [pos_sorted,idx_cds] = sort(pos_cds);
    % ------------------------------------------------------------
    % Position of x and delimiters after x, y and z
    % pos of x coord, and delimiters after x and y coords:
    cd1_start = pos_cds( idx_cds(1) );
    dlm1_str   = cds_fmt( pos_cds( idx_cds(1) )+3:pos_cds( idx_cds(2) )-1);
    dlm2_str   = cds_fmt( pos_cds( idx_cds(2) )+3:pos_cds( idx_cds(3) )-1);
    % character after z - if any
    if length(cds_fmt)>pos_cds( idx_cds(3) )+2
        fin_str = cds_fmt( pos_cds( idx_cds(3) )+3 );
    else
        fin_str =[];
    end
    % ------------------------------------------------------------
    % Strings of coordinates:
    str_length = length( str );
    %Str of coordinate1 using Pos of delimiter12
    pos_dlm1  = cd1_start + strfind( str(cd1_start:str_length),dlm1_str)-1;
    str_cds{ idx_cds(1) } = str( cd1_start:pos_dlm1(1)-1 );
    %Str of coordinate2 using Pos of delimiter23
    cd2_start = pos_dlm1(1) + length(dlm1_str);
    pos_dlm2  = cd2_start + strfind( str(cd2_start:str_length),dlm2_str)-1;
    str_cds{ idx_cds(2) } = str( cd2_start:pos_dlm2(1)-1 );
    %Str of coordinate3 using fin_str
    cd3_start = pos_dlm2(1) + length(dlm2_str);
    if isempty(fin_str)~=1%there is a character after third coordinate!
        relpos_dlmfin =  strfind(str( cd3_start:str_length), fin_str );
        if isempty( relpos_dlmfin )~=1
            pos_dlmfin = cd3_start + relpos_dlmfin(1) -1;
            str_cds{ idx_cds(3) } = str( cd3_start:pos_dlmfin(1)-1 );
        else
            str_cds{ idx_cds(3) } = str( cd3_start:str_length );
        end
    else %isempty(dlm_fin)==1 || isempty( relpos_dlmfin )==1
        str_cds{ idx_cds(3) } = str( cd3_start:str_length );
    end
    str_cds = regexprep( str_cds , ',', '.');
    cds = str2double( str_cds );
catch
    cds = [NaN NaN NaN];
    return
end












