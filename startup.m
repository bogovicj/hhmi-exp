% clear classes; clear java; clc;
global FIJIPATH;        % locaction of fiji install
global IMGLIB2PATH;     % location of imglib2

if(~isdeployed)
    javaaddpath('/groups/jain/home/bogovicj/lib_workspace/hhmi-exp/target/classes');
    javaaddpath('/groups/jain/home/bogovicj/dev/imglibFork/pom-imglib2/algorithms/core/target/classes');
    javaaddpath('/groups/saalfeld/home/bogovicj/dev/main/patchSearchStandalone/bin');
end

FIJIPATH = '/groups/jain/home/bogovicj/packages/Fiji.app';
IMGLIB2PATH = '/groups/jain/home/bogovicj/packages/Fiji.app';

% these lines require customization:
global JBJAINPATH;      % location of john's code repository
global MGDMPATH;        % location of mgdm code repository
global HHMIPATH;
global DICTPATH;
global DFEVAL_DIR;

global EXPPATH;         % where to save experiment output
global VERBOSE;         % make functions chattier

DFEVAL_DIR    = '/nobackup/saalfeld/john/tmp';
EXPPATH    = '/groups/saalfeld/home/bogovicj/reseach/exp';

JBJAINPATH = '/groups/jain/home/bogovicj/workspace/janelia_jainlab';
MGDMPATH = '/groups/jain/home/bogovicj/workspace/mgdm';

HHMIPATH = '/groups/saalfeld/home/bogovicj/dev/main/hhmi-exp';
DICTPATH = '/groups/saalfeld/home/bogovicj/dev/main/dictionary-feature-classification';

VERBOSE    = false;

if(~isdeployed)
    addpath(genpath(JBJAINPATH));
    
    addpath('/groups/jain/home/bogovicj/libraries/matlab/export_fig');
    addpath('/groups/jain/home/bogovicj/libraries/matrix2latexMatlab');
    addpath('/groups/jain/home/bogovicj/libraries/ext');
    addpath(genpath(fullfile(HHMIPATH,'matlab')));
    addpath(genpath(fullfile(DICTPATH,'matlab')));
    
end

% setdbprefs('FetchBatchSize','1000')
% setdbprefs('FetchInBatches','no')

% rng('shuffle')

%% build static java classpath from maven dependencies
% cmd = 'dependency:build-classpath | grep ''Dependencies classpath'' -A 1 | tr '':'' ''\n'''
% system( cmd );