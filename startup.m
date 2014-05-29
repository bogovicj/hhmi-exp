% clear classes; clear java; clc;
global FIJIPATH;        % locaction of fiji install
global IMGLIB2PATH;     % location of imglib2

FIJIPATH = '/groups/jain/home/bogovicj/packages/Fiji.app';
IMGLIB2PATH = '/groups/jain/home/bogovicj/packages/Fiji.app';

% these lines require customization:
global JBJAINPATH;      % location of john's code repository
global MGDMPATH;        % location of mgdm code repository

global EXPPATH;         % where to save experiment output
global VERBOSE;         % make functions chattier

EXPPATH    = '/groups/jain/home/bogovicj/saalfeldExp';

JBJAINPATH = '/groups/jain/home/bogovicj/workspace/janelia_jainlab';
MGDMPATH = '/groups/jain/home/bogovicj/workspace/mgdm';
 
VERBOSE    = false;

if(~isdeployed)
    addpath(genpath(JBJAINPATH));
    
    addpath('/groups/jain/home/bogovicj/libraries/matlab/export_fig');
    addpath('/groups/jain/home/bogovicj/libraries/matrix2latexMatlab');
    addpath('/groups/jain/home/bogovicj/libraries/ext');
    addpath(genpath('/groups/jain/home/bogovicj/lib_workspace/hhmi-exp/matlab'));
end

% addpath('/groups/jain/home/bogovicj/pubDrafts/paper_flash_features/code/handselected_exps');

% setdbprefs('FetchBatchSize','1000')
% setdbprefs('FetchInBatches','no')


% rng('shuffle')
% 
