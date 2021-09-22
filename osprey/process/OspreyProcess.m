function [MRSCont] = OspreyProcess(MRSCont)
%% [MRSCont] = OspreyProcess(MRSCont)
%   This function pre-processes MRS data from all major vendors.
%   Data is read from the provided input filenames. It is shaped,
%   preprocessed, aligned, etc. according to the type of sequence
%   (un-edited data, MEGA-edited (ON/OFF), HERMES/HERCULES (A/B/C/D),
%   etc.).
%
%   USAGE:
%       MRSCont = OspreyProcess(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-19)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-02-19: First version of the code.

outputFolder = MRSCont.outputFolder;
diary(fullfile(outputFolder, 'LogFile.txt'));

% Checking for version, toolbox, and previously run modules
osp_CheckRunPreviousModule(MRSCont, 'OspreyProcess');
[~,MRSCont.ver.CheckOsp ] = osp_Toolbox_Check('OspreyProcess',MRSCont.flags.isGUI);


% Post-process raw data depending on sequence type
if ~MRSCont.flags.skipProcess
    if ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
        if MRSCont.flags.isUnEdited
            [MRSCont] = osp_processUnEdited(MRSCont);
        elseif MRSCont.flags.isMEGA
            [MRSCont] = osp_processMEGA(MRSCont);
        elseif MRSCont.flags.isHERMES
            [MRSCont] = osp_processHERMES(MRSCont);
        elseif MRSCont.flags.isHERCULES
            % For now, process HERCULES like HERMES data
            [MRSCont] = osp_processHERCULES(MRSCont);
        elseif MRSCont.flags.isDWMRS
            [MRSCont] = osp_processDW(MRSCont);
        else
            msg = 'No flag set for sequence type!';
            fprintf(msg);
            error(msg);
        end
    else
        [MRSCont] = osp_processMultiVoxel(MRSCont);
    end
    
else
    % We want to be able to skip processing under certain circumstances
    % (for example, if the user provides pre-existing LCModel .raw files
    % that are ready to fit out of the box).
    [MRSCont] = osp_processSkip(MRSCont);
end




% Gather some more information from the processed data;
SubSpecNames = fieldnames(MRSCont.processed);
NoSubSpec = length(fieldnames(MRSCont.processed));
if iscell(MRSCont.processed.A{1})
    NoExtras = length(MRSCont.processed.A{1});
else
    for ss = 1 : NoSubSpec
        for kk = 1 : MRSCont.nDatasets
            temp_sz(1,kk)= MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1);
            temp_sz_sw{1,kk} = ['np_sw_' num2str(MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1)) '_' num2str(MRSCont.processed.(SubSpecNames{ss}){1,kk}.spectralwidth)];
        end
        [MRSCont.info.(SubSpecNames{ss}).unique_ndatapoint_spectralwidth,MRSCont.info.(SubSpecNames{ss}).unique_ndatapoint_spectralwidth_ind,~]  = unique(temp_sz_sw,'Stable');
        [MRSCont.info.(SubSpecNames{ss}).max_ndatapoint,MRSCont.info.(SubSpecNames{ss}).max_ndatapoint_ind] = max(temp_sz);
    end
end


%% If DualVoxel or MRSI we want to extract y-axis scaling
% Creates y-axis range to align the process plots between datasets

if MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI
    MRSCont.plot.processed.match = 1; % Scaling between datasets is turned off by default
else
    MRSCont.plot.processed.match = 0; % Scaling between datasets is turned off by default
end

if ~MRSCont.flags.isDWMRS
    MRSCont = osp_scale_yaxis(MRSCont,'OspreyProcess');
end

%% Clean up and save
% Set exit flags and reorder fields
MRSCont.flags.didProcess    = 1;
diary off
[MRSCont]                   = osp_orderProcessFields(MRSCont);

% Store data quality measures in csv file
if MRSCont.flags.isUnEdited
    names = {'NAA_SNR','NAA_FWHM','residual_water_ampl','freqShift'};
    subspec = {'A'};
    if MRSCont.flags.hasRef
        names = {'NAA_SNR','NAA_FWHM','water_FWHM','residual_water_ampl','freqShift'};
    end
elseif MRSCont.flags.isMEGA
    names = {'NAA_SNR','NAA_FWHM','residual_water_ampl','freqShift'};
    subspec = {'A'};
    if MRSCont.flags.hasRef
        names = {'NAA_SNR','NAA_FWHM','water_FWHM','residual_water_ampl','freqShift'};
    end
elseif MRSCont.flags.isHERMES
    names = {'NAA_SNR','NAA_FWHM','residual_water_ampl','freqShift'};
    subspec = {'sum'};
    if MRSCont.flags.hasRef
        names = {'NAA_SNR','NAA_FWHM','water_FWHM','residual_water_ampl','freqShift'};
    end
elseif MRSCont.flags.isHERCULES
    % For now, process HERCULES like HERMES data
    names = {'NAA_SNR','NAA_FWHM','residual_water_ampl','freqShift'};
    subspec = {'sum'};
    if MRSCont.flags.hasRef
        names = {'NAA_SNR','NAA_FWHM','water_FWHM','residual_water_ampl','freqShift'};
    end
elseif MRSCont.flags.isDWMRS
    names = {'NAA_SNR','NAA_FWHM','residual_water_ampl'};
    subspec = {'A'};
    if MRSCont.flags.hasRef
        names = {'NAA_SNR','NAA_FWHM','water_FWHM','residual_water_ampl'};
    end
else
    msg = 'No flag set for sequence type!';
    fprintf(fileID,msg);
    error(msg);
end

if ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI && ~MRSCont.flags.isDWMRS
    if ~MRSCont.flags.hasRef
        QM = horzcat(MRSCont.QM.SNR.(subspec{1})',MRSCont.QM.FWHM.(subspec{1})',MRSCont.QM.res_water_amp.(subspec{1})',MRSCont.QM.freqShift.(subspec{1})');
    else
        QM = horzcat(MRSCont.QM.SNR.(subspec{1})',MRSCont.QM.FWHM.(subspec{1})',MRSCont.QM.FWHM.ref',MRSCont.QM.res_water_amp.(subspec{1})',MRSCont.QM.freqShift.(subspec{1})');
    end
    MRSCont.QM.tables = array2table(QM,'VariableNames',names);
    writetable(MRSCont.QM.tables,[outputFolder '/QM_processed_spectra.csv']);
end

% Optional:  Create all pdf figures
if MRSCont.opts.savePDF && ~MRSCont.flags.isDWMRS
    osp_plotAllPDF(MRSCont, 'OspreyProcess')
end

% Optional: write edited files to LCModel .RAW files
if MRSCont.opts.saveLCM && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    [MRSCont] = osp_saveLCM(MRSCont);
end

% Optional: write edited files to jMRUI .txt files
if MRSCont.opts.savejMRUI && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI && ~MRSCont.flags.isDWMRS
    [MRSCont] = osp_saveJMRUI(MRSCont);
end

% Optional: write edited files to vendor specific format files readable to
% LCModel and jMRUI
% SPAR/SDAT if Philips
% RDA if Siemens
if MRSCont.opts.saveVendor && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI && ~MRSCont.flags.isDWMRS
    [MRSCont] = osp_saveVendor(MRSCont);
end

% Optional: write edited files to NIfTI-MRS format
if MRSCont.opts.saveNII && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI && ~MRSCont.flags.isDWMRS
    [MRSCont] = osp_saveNII(MRSCont);
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
