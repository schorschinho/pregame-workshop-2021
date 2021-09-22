function [MRSCont] = OspreyPosthocModel(MRSCont)
%% [MRSCont] = OspreyPosthocModel(MRSCont)
%   This function takes modeling results from OspreyFit or OspreyQuantify,
%   and performs additional post-hoc modeling steps to derive quantities
%   like apparent diffusion coefficients.
%
%   USAGE:
%       MRSCont = OspreyPosthocModel(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2021-06-19)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)

outputFolder = MRSCont.outputFolder;
diary(fullfile(outputFolder, 'LogFile.txt'));

% Checking for version, toolbox, and previously run modules
%osp_CheckRunPreviousModule(MRSCont, 'OspreyPosthocModel');
%[~,MRSCont.ver.CheckOsp ] = osp_Toolbox_Check('OspreyPosthocModel',MRSCont.flags.isGUI);

collectedResults = [];
if isfield(MRSCont.opts, 'posthoc')
    % Collect the independent variable
    if isfield(MRSCont.opts.posthoc, 'x')
        if isfield(MRSCont.opts.posthoc.x, 'values')
            if ~isempty(MRSCont.opts.posthoc.x.values)
                if ischar(MRSCont.opts.posthoc.x.values)
                    eval(['x = ' MRSCont.opts.posthoc.x.values ';']);
                else
                    x = MRSCont.opts.posthoc.x.values;
                end
            end
        end
    end
    
    % Collect the dependent variable
    if isfield(MRSCont.opts.posthoc, 'y')
        if isfield(MRSCont.opts.posthoc.y, 'params')
            if ~isempty(MRSCont.opts.posthoc.y.params)
                switch MRSCont.opts.posthoc.y.params
                    case 'amplitudes'
                        targetField = 'ampl';
                end
                
                % Collect basis function names and their results
                nSpectra    = length(MRSCont.fit.results.off.fitParams{1});
                bfNames     = MRSCont.fit.results.off.fitParams{1}{1}.name;
                collectedResults.names = bfNames;
                for ll = 1:length(bfNames)
                    for nn = 1:nSpectra
                        collectedResults.(MRSCont.opts.posthoc.y.params)(ll,nn) = MRSCont.fit.results.off.fitParams{1}{nn}.(targetField)(ll);
                    end
                end
            end
        end
    end
end

% Do the post-hoc modeling
for ll = 1:length(bfNames)
    x_in = x;
    y_in = collectedResults.(MRSCont.opts.posthoc.y.params)(ll,:);
    myfittype = fittype(MRSCont.opts.posthoc.model, 'dependent', {'y_in'}, 'independent', {'x'});
    myfit{ll} = fit(x_in,y_in',myfittype);
end

% Save the output structure to the output folder
% Determine output folder
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

end
