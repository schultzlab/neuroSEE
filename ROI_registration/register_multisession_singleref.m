function [A_union, A2, assignments, matchings, templates_shifted, cmatched_ROIs, cnonmatched_1, cnonmatched_2, R ] = ...
            register_multisession_singleref( A, options, templates, options_mc, figname_pref, figclose )
% REGISTER_MULTISESSION - register ROIs from multiple recording sessions
%
%   [A_UNION, ASSIGNMENTS, MATCHINGS] = REGISTER_MULTISESSION(A,...
%          OPTIONS, TEMPLATES, OPTIONS_MC)
%
% Register ROIs across multiple sessions using an intersection over union metric
% and the Hungarian algorithm for optimal matching. Registration occurs by 
% aligning session 2 to session 1, keeping the union of the matched and 
% non-matched components. Then session 3 is aligned to 1 and so on.
%
% INPUTS:
% A:                      cell array with spatial components from each
%                           session arranged chronologically
% options                 parameter structure with inputs:
%       d1:               number of rows in FOV
%       d2:               number of columns in FOV
%       d3:               number of planes in FOV (default: 1)
%       dist_maxthr:      threshold for turning spatial components into binary masks (default: 0.1)
%       dist_exp:         power n for distance between masked components: dist = 1 - (and(m1,m2)/or(m1,m2))^n (default: 1)
%       dist_thr:         threshold for setting a distance to infinity. (default: 0.5)
%       dist_overlap_thr: overlap threshold for detecting if one ROI is a subset of another (default: 0.8)
%       plot_reg:         create a contour plot of registered ROIs
% templates:              cell array with a reference image from each session
% options_mc:             normcorre options structure (for template alignment) 
%                           (structure array the length of AA)

% OUTPUTS:
% A_union:                union of ROIs aligned to session 1 (for matched
%                               pairs the ROIs from session # 1 are kept)
% assignments:            matrix of size # of total distinct components x # sessions
%                               element [i,j] = k if component k from session j is mapped to component
%                               i in the A_union matrix. If there is no match the value is NaN
% matchings:              cell array. Each entry matchings{i}(j) = k means that component j from session
%                               i is represented by component k in A_union

% Ann's addition:
if nargin < 7, figclose = true; end

defoptions = CNMFSetParms;
if ~exist('options','var'); options = defoptions; end

if ~exist('templates','var') || ~exist('options_mc','var')
    warning('Some required inputs for aligning ROIs before registering are missing. Skipping alignment');
    align_flag = false;
else
    align_flag = true;
end

if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); options.d1 = d1; end 
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); options.d2 = d2; end 
if ~isfield(options,'d3') || isempty(options.d3); options.d3 = 1; end
if ~isfield(options,'dist_maxthr') || isempty(options.maxthr); options.dist_maxthr = 0.15; end
if ~isfield(options,'dist_exp') || isempty(options.dist_exp); options.dist_exp = 1; end
if ~isfield(options,'dist_thr') || isempty(options.dist_thr); options.dist_thr = 0.5; end
if ~isfield(options,'dist_overlap_thr') || isempty(options.dist_overlap_thr); options.dist_overlap_thr = 0.8; end
options.plot_reg = true;

if ~align_flag
    templates = {[]};
end
n_sessions = length(A);
if length(templates) == 1
    for i = 2:n_sessions
        templates{i} = templates{1};
    end
end

% siz = [options.d1,options.d2,options.d3];
if length(options_mc) == 1
    O = struct('O', repmat(options_mc, length(A)-1, 1));
    options_mc = O.O;
end
for n = 1:length(A)-1
    options_mc(n).r.correct_bidir = false;
    options_mc(n).nr.correct_bidir = false;
end
A_union{1} = A{1};
matchings{1} = 1:size(A{1},2);

% This version aligns all to one global template (1). In the case of
% matched pairs between 1 and 2, it keeps 1.
for s = 2:n_sessions
   if ~isempty(figname_pref)
        fname_fig = [figname_pref '_s' num2str(s) '_regto_1'];
   else
       fname_fig = [];
   end
   disp(['session: ' num2str(s)])
   [ ~, ~, ~, As_shifted, ~, ~, template_shifted ] = ...
      register_ROIs( A_union{s-1}, A{s}, options, templates{1}, templates{s}, options_mc(s-1).r, [], figclose );
   disp('rigid done')
   [ matched_ROIs, nonmatched_1, nonmatched_2, A2{s-1}, R{s-1}, A_union{s}, templates_shifted{s-1} ] = ...
      register_ROIs( A_union{s-1}, As_shifted, options, templates{1}, template_shifted, options_mc(s-1).nr, fname_fig, figclose );
   new_match = zeros(1,size(A{s},2));
   new_match(matched_ROIs(:,2)) = matched_ROIs(:,1);
   new_match(nonmatched_2) = size(A_union{s},2)-numel(nonmatched_2)+1 : size(A_union{s},2);
   matchings{s} = new_match;
   
   cmatched_ROIs{s-1} = matched_ROIs;
   cnonmatched_1{s-1} = nonmatched_1;
   cnonmatched_2{s-1} = nonmatched_2;
end

assignments = NaN(size(A_union{n_sessions},2), n_sessions);
for s = 1:n_sessions
    assignments(matchings{s}, s) = 1:length(matchings{s});
end