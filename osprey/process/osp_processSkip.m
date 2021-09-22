function [MRSCont] = osp_processSkip(MRSCont)
%% [MRSCont] = osp_processSkip(MRSCont)
%   This function skips all regular post-processing steps, and just
%   performs QA control on the raw, averaged datasets that were passed by
%   OspreyLoad.
%
%   USAGE:
%       [MRSCont] = osp_processSkip(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2021-08-16)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)

warning('off','all');

%% Loop over all datasets
refProcessTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end

for kk = 1:MRSCont.nDatasets
    [~] = printLog('OspreyProcess', kk, MRSCont.nDatasets ,progressText ,MRSCont.flags.isGUI, MRSCont.flags.isMRSI); 
    
    if ~(MRSCont.flags.didProcess == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'processed') && (kk > length(MRSCont.processed.A)))  || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)

        % Sham frequency/phase correction (just to have some trivial QA)
        raw = MRSCont.raw{kk};
        raw.flags.averaged  = 1;
        raw.dims.averages   = 0;
        raw.specReg.fs              = 0; % save align parameters
        raw.specReg.phs             = 0; % save align parameters
        raw.specReg.weights         = 1; % save align parameters
        raw.watersupp.amp = 0;
        driftPre = op_measureDrift(raw);
        driftPost = driftPre;
        MRSCont.processed.A{kk}     = raw;
        refShift = 0;
        
        
        %%% 2. GET REFERENCE DATA %%%
        % If there are reference scans, load them here.
        if MRSCont.flags.hasRef
            MRSCont.processed.ref{kk} = MRSCont.raw_ref{kk};        % Get the kk-th dataset
        end
        
        %%% 3. GET MM DATA %%%
        if MRSCont.flags.hasMM %re_mm
            MRSCont.processed.mm{kk}  = MRSCont.raw_mm{kk};     	% Get the kk-th dataset  %re_mm
        end

        %%% 4. GET SHORT-TE WATER DATA %%%
        if MRSCont.flags.hasWater
            MRSCont.processed.w{kk}   = MRSCont.raw_w{kk};          % Get the kk-th dataset

        end

        %%% 5. QUALITY CONTROL PARAMETERS %%%
        SubSpec = {'A'};       
        SNRRange = {[1.8,2.2]};
        if MRSCont.flags.hasMM
            SubSpec{end+1} = 'mm';
            SNRRange{end+1} = [0.7,1.1];
        end
        if MRSCont.flags.hasRef
            SubSpec{end+1} = 'ref';
            SNRRange{end+1} = [4.2,5.2];
        end
        if MRSCont.flags.hasWater
            SubSpec{end+1} = 'w';
            SNRRange{end+1} = [4.2,5.2];
        end
        % Calculate some spectral quality metrics here;
        for ss = 1 : length(SubSpec)          
            MRSCont.QM.SNR.(SubSpec{ss})(kk)    = op_getSNR(MRSCont.processed.(SubSpec{ss}){kk});       
            MRSCont.QM.FWHM.(SubSpec{ss})(kk)   = op_getLW(MRSCont.processed.(SubSpec{ss}){kk},SNRRange{ss}(1),SNRRange{ss}(2)); % in Hz       
            if ~(strcmp(SubSpec{ss},'ref') || strcmp(SubSpec{ss},'w') || strcmp(SubSpec{ss},'mm'))
                 MRSCont.QM.drift.pre.(SubSpec{ss}){kk}  = driftPre;
                MRSCont.QM.drift.post.(SubSpec{ss}){kk} = driftPost;
                MRSCont.QM.freqShift.(SubSpec{ss})(kk)  = refShift;       
                MRSCont.QM.res_water_amp.(SubSpec{ss})(kk) = sum(MRSCont.processed.(SubSpec{ss}){kk}.watersupp.amp);  
                if strcmp(SubSpec{ss},'diff1') || strcmp(SubSpec{ss},'sum')
                    MRSCont.QM.drift.pre.(SubSpec{ss}){kk}  = reshape([MRSCont.QM.drift.pre.A'; MRSCont.QM.drift.pre.B'], [], 1)';
                    MRSCont.QM.drift.post.(SubSpec{ss}){kk} = reshape([MRSCont.QM.drift.post.A'; MRSCont.QM.drift.post.B'], [], 1)';
                end
                MRSCont.QM.drift.pre.AvgDeltaCr.(SubSpec{ss})(kk) = mean(MRSCont.QM.drift.pre.(SubSpec{ss}){kk} - 3.02);
                MRSCont.QM.drift.post.AvgDeltaCr.(SubSpec{ss})(kk) = mean(MRSCont.QM.drift.pre.(SubSpec{ss}){kk} - 3.02);
            end
        end              
    end
end
time = toc(refProcessTime);
[~] = printLog('done',time,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 

%%% 10. SET FLAGS %%%
MRSCont.flags.avgsAligned       = 1;
MRSCont.flags.averaged          = 1;
MRSCont.flags.ECCed             = 1;
MRSCont.flags.waterRemoved      = 1;
MRSCont.runtime.Proc = time;
% Close any remaining open figures
close all;

end
