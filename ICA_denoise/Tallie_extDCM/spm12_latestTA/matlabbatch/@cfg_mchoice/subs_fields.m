function fnames = subs_fields(item)

% function fnames = subs_fields(item)
% This function works as a "class-based switch" to return the value of
% the private mysubs_fields function for the appropriate class. 
% This function is identical for all classes derived from cfg_item, but
% it has to be in the class directory to access the proper private
% function mysubs_fields.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subs_fields.m 4867 2012-08-30 13:04:51Z volkmar $

rev = '$Rev: 4867 $'; %#ok

fnames = mysubs_fields;