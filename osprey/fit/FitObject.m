classdef FitObject < handle
    
    
    properties
        % Everything we need to perform a fit and store the results.
        Data = struct('fids', [], 'DwellTime', [] , 'SpectralWidth', [], 'txfrq', [], 't', [], 'ppm', []);
        BasisSets = struct('fids', [], 'names', [], 'includeInFit', []);
        BaselineBasis = struct('specs', []);
        Options = struct;
        Model = struct;
        Results = struct;
    end
    
    
    
    methods
        
        function obj = FitObject(data, basis, options)
            % class constructor
            if(nargin > 0) % don't have to initialize the fields
                
                
                %%% DATA %%%
                % Copy the information necessary to create appropriate time-
                % and frequency-domain signals:
                obj.Data.fids            = data.fids;
                obj.Data.DwellTime       = data.dwelltime;
                obj.Data.SpectralWidth   = data.spectralwidth;
                obj.Data.txfrq           = data.txfrq;
                obj.Data.t               = data.t;
                
                % Calculate the ppm axis
                npts            = size(data.fids, 1);
                spectralwidth   = data.spectralwidth;
                f   = (-spectralwidth/2) + (spectralwidth/(2*npts)) : spectralwidth/npts : (spectralwidth/2) - (spectralwidth/(2*npts));
                obj.Data.ppm = f / (data.txfrq * 1e-6) + 4.68;
            
                
                %%% BASIS SET %%%
                % Check that basis set and data have the same resolution
                if abs(basis.dwelltime - obj.Data.DwellTime) > eps
                    warning('Dwell time does not agree between basis set (%5.2e) and data (%5.2e).', obj.Data.DwellTime, basis.dwelltime);
                    fprintf('Resampling the basis set for you. \n');
                    basis = fit_resampleBasis(data, basis);
                end
                if round(basis.spectralwidth) ~= round(obj.Data.SpectralWidth)
                    warning('Spectral width does not agree between basis set (%5.2e) and data (%5.2e).', obj.Data.DwellTime, basis.spectralwidth);
                    fprintf('Resampling the basis set for you. \n');
                    basis = fit_resampleBasis(data, basis);
                end
                
                % Copy the necessary information
                obj.BasisSets.fids  = basis.fids;
                obj.BasisSets.names = basis.name;
                obj.BasisSets.includeInFit = ones(size(basis.name));
                
                
                %%% OPTIONS %%%
                % Initialize an empty container
                if ~isfield(options, 'optimDomain')
                    options.optimDomain = 'FD'; % FD, TD, or FDTD
                else
                    if ~(strcmpi(options.optimDomain,'FD') ||...
                            strcmpi(options.optimDomain,'TD') ||...
                            strcmpi(options.optimDomain,'FDTD'))
                        error('Invalid optimization domain specification (options.optimDomain): (%s).', options.optimDomain)
                    end
                end
                
                if ismember(options.optimDomain, {'FD', 'FDTD'})
                    if ~isfield(options, 'optimFreqFitRange')
                        options.optimFreqFitRange = [0.5 4.0];
                    end
                end
                
                if ismember(options.optimDomain, {'TD', 'FDTD'})
                    if ~isfield(options, 'optimTimeFitRange')
                        options.optimTimeFitRange = [0 1];
                    end
                end
                
                if ~isfield(options, 'optimSignalPart')
                    options.optimSignalPart = 'R'; % R, I, or RI
                end
                
                
                if ~isfield(options, 'parametrizations')
                    % Initialize phi0 as constant with value 0
                    options.parametrizations.ph0.fun     = 'free';
                    options.parametrizations.ph0.gradfun = 'free';
                    options.parametrizations.ph0.lb      = -pi;
                    options.parametrizations.ph0.ub      = pi;
                    options.parametrizations.ph0.init    = 0;
                    
                    % Initialize phi1 as constant with value 0
                    options.parametrizations.ph1.fun     = 'free';
                    options.parametrizations.ph1.gradfun = 'free';
                    options.parametrizations.ph1.lb      = -pi/30;
                    options.parametrizations.ph1.ub      = pi/30;
                    options.parametrizations.ph1.init    = 0;
                    
                    % Initialize Gaussian LB as constant with value [0.04 *
                    % hz/ppm]
                    options.parametrizations.gaussLB.fun     = 'free';
                    options.parametrizations.gaussLB.gradfun = 'free';
                    options.parametrizations.gaussLB.lb      = 0;
                    options.parametrizations.gaussLB.ub      = sqrt(5000);
                    options.parametrizations.gaussLB.init    = 0.04 * obj.Data.txfrq*1e-6;
                    
                    % Initialize Lorentzian LB as constant with value 2 Hz
                    options.parametrizations.lorentzLB.fun     = 'free';
                    options.parametrizations.lorentzLB.gradfun = 'free';
                    options.parametrizations.lorentzLB.lb      = 0;
                    options.parametrizations.lorentzLB.ub      = 100;
                    options.parametrizations.lorentzLB.init    = 2;
                    
                    % Initialize frequency shifts as constant with value 0 Hz
                    options.parametrizations.freqShift.fun     = 'free';
                    options.parametrizations.freqShift.gradfun = 'free';
                    options.parametrizations.freqShift.lb      = -20;
                    options.parametrizations.freqShift.ub      = 20;
                    options.parametrizations.freqShift.init    = 0;
                    
                    % Initialize metabolite amplitudes as free with value 0
                    options.parametrizations.metAmpl.fun     = 'free';
                    options.parametrizations.metAmpl.gradfun = 'free';
                    options.parametrizations.metAmpl.lb      = 0;
                    options.parametrizations.metAmpl.ub      = Inf;
                    options.parametrizations.metAmpl.init    = 0;
                    
                    % Initialize baseline amplitudes as free with value 0
                    options.parametrizations.baseAmpl.fun     = 'free';
                    options.parametrizations.baseAmpl.gradfun = 'free';
                    options.parametrizations.baseAmpl.lb      = -Inf;
                    options.parametrizations.baseAmpl.ub      = Inf;
                    options.parametrizations.baseAmpl.init    = 0;
                    
                    % Initialize x (the external dependency vector) as natural numbers
                    options.parametrizations.x.values = [1:1:size(obj.Data.fids,2)];
                    options.parametrizations.x.name   = 'independentVariable';
                    
                end
                
                % Save the property struct
                obj.Options = options;
                
                
                %%% CREATE BASELINE SPLINE BASIS %%%
                % Combine real and imaginary part to form a complex spline array.
                % Use the new, corrected function from here on
                fitRange    = obj.Options.optimFreqFitRange;
                dkntmn      = obj.Options.dkntmn;
                [splineArray] = fit_makeSplineBasis_New(data, fitRange, dkntmn);

                % Save into the property
                obj.BaselineBasis = splineArray;
            
            end
            
        end


        function includeBasisFunctionInFit(obj, input)
            % Sets the 'includeInFit' flag in the 'BasisSets' property to
            % 1, meaning that the corresponding basis function will be
            % included in the fit.
            
            if ischar(input)
                input = {input};
            end
            
            for ii = length(input)
                % Check which metabolites are available in the basis set
                % and match the input
                [metsToInclude, ~, ~] = intersect(obj.BasisSets.names, input, 'stable');
                for rr = 1:length(metsToInclude)
                    idxToInclude = find(strcmp(metsToInclude{rr}, obj.BasisSets.names));
                    obj.BasisSets.includeInFit(idxToInclude) = 1;
                end
                
            end
        end
        
        function excludeBasisFunctionFromFit(obj, input)
            % Sets the 'includeInFit' flag in the 'BasisSets' property to
            % 0, meaning that the corresponding basis function will be
            % excluded from in the fit.
            
            if ischar(input)
                input = {input};
            end
            
            if strcmpi(input, 'all')
                input = obj.BasisSets.names;
            end
            
            for ii = length(input)
                % Check which metabolites are available in the basis set
                % and match the input
                [metsToInclude, ~, ~] = intersect(obj.BasisSets.names, input, 'stable');
                for rr = 1:length(metsToInclude)
                    idxToInclude = find(strcmp(metsToInclude{rr}, obj.BasisSets.names));
                    obj.BasisSets.includeInFit(idxToInclude) = 0;
                end
                
            end
        end
        
        function tableOut = showBasisFunctions(obj)
            % Returns a table showing the basis functions included in the
            % basis set, and whether they will be used in the fit.
            tableOut = table(obj.BasisSets.names', obj.BasisSets.includeInFit', ...
                'VariableNames', {'Basis function', 'Include in fit?'});
        end
        
        
        
        function [ph0, gaussLB, lorentzLB, freqShift, metAmpl, baseAmpl] = initFit(obj)
            
            % This method runs a quick and dirty fit to get decent starting values.
            
            % Collect the basis functions
            basisSet        = obj.BasisSets;
            baselineBasis   = obj.BaselineBasis;
            data            = obj.Data.fids;
            ppm             = obj.Data.ppm;
            t               = obj.Data.t;
            fitRange        = obj.Options.optimFreqFitRange;
            [indMin, indMax] = FitObject.ppmToIndex(ppm, fitRange);
            
            % Only use basis functions that are included
            basisSet.fids   = basisSet.fids(:, logical(basisSet.includeInFit));
            
            % Create x0
            % ph0, gaussLB, lorentzLB, freqShift
            gaussInit = 0.04*obj.Data.txfrq * 1e-6;
            lorentzInit = 2;
            x0 = [0,   log(gaussInit), log(lorentzInit), 0];
            lb = [-pi, 0,              0,                -Inf];
            ub = [pi,  Inf,            Inf,              Inf];
            
            
            tstart=tic;
            % Prepare the function wrapper
            fcn  = @(x) FitObject.initialFitLossFunction(x, data, basisSet, baselineBasis, ppm, t, fitRange);
            xk = fmincon(fcn, x0, [], [], [], [], lb, ub);
            
            time = toc(tstart);
            
            % Translate output
            fids = basisSet.fids;
            nBasisFcts = size(fids, 2);
            nBaselineComps = size(baselineBasis, 2);
            
            ph0 = xk(1);
            gaussLB = exp(xk(2));
            lorentzLB = exp(xk(3));
            freqShift = xk(4);
            
            specs = FitObject.transformBasis(fids, gaussLB, lorentzLB, freqShift, t);
            
            % append basis set with the spline baseline basis functions
            fullBasis = cat(2, specs, baselineBasis);
            % phase
            fullBasisPhased = exp(-1j*ph0) .* fullBasis;
            
            % crop to fit range
            data = fftshift(fft(data, [], 1), 1);
            dataCrop = data(indMin:indMax,:);
            fullBasisCrop = fullBasisPhased(indMin:indMax,:);
            
            % concatenate real and imaginary
            dataCrop = cat(1, real(dataCrop), imag(dataCrop));
            fullBasisCrop = cat(1, real(fullBasisCrop), imag(fullBasisCrop));

            % solve the linear equation system
            beta = real(pinv(fullBasisCrop) * dataCrop);
            
            % project metabolite coefficients to >=0
            beta(beta(1:nBasisFcts) < 0) = 0;
            
            % make spectrum
            prediction = fullBasis * beta;
            
            % separate metAmpl and baseAmpl
            metAmpl = beta(1:nBasisFcts);
            baseAmpl = beta(nBasisFcts+1:nBasisFcts+nBaselineComps);
            
            % Save modeling results
            dataReversePhased = exp(1j*ph0).*data;
            parsOut.ph0 = ph0;
            parsOut.ph1 = 0;
            parsOut.gaussLB = gaussLB;
            parsOut.lorentzLB = lorentzLB;
            parsOut.freqShift = freqShift;
            parsOut.metAmpl = metAmpl;
            parsOut.baseAmpl = baseAmpl;
            obj.Model.init.fit.data = dataReversePhased;
            obj.Model.init.fit.fit = prediction;
            obj.Model.init.fit.baseline = baselineBasis*baseAmpl;
            obj.Model.init.fit.residual = dataReversePhased-prediction;
            obj.Model.init.fit.metabs   = fullBasis(:,1:nBasisFcts).*metAmpl';
            obj.Model.init.time = time;
            obj.Model.init.parsOut = parsOut;
            
        end
            
        function obj = createModel(obj)
            
            % This is the heart and soul of this class. This function takes
            % the parametrizations defined in the options, and translates
            % them into the appropriate loss and gradient functions.
            
            % Collect the basis functions
            basisSet        = obj.BasisSets;
            baselineBasis   = obj.BaselineBasis;
            data            = obj.Data.fids;
            ppm             = obj.Data.ppm;
            t               = obj.Data.t;
            fitRange        = obj.Options.optimFreqFitRange;
            
            % Only use basis functions that are included
            basisSet.fids   = basisSet.fids(:, logical(basisSet.includeInFit));
            
            % Create x0, lb, ub vectors by iteratively calling the
            % parameter-class-specific initialization
            pars = {'ph0', 'ph1', 'gaussLB', 'metAmpl', 'freqShift', 'lorentzLB', 'baseAmpl'};
            parsInit = [];
            parslb = [];
            parsub = [];
            for pp = 1:length(pars)
                [~, parslb, parsub] = initializeParameters(obj, parsInit, parslb, parsub, pars{pp});
                [parsInit] = getStartingValuesFromInitialFit(obj, parsInit, pars{pp});
            end
            
            x0 = FitObject.pars2x(parsInit);
            lb = FitObject.pars2x(parslb);
            ub = FitObject.pars2x(parsub);
            
            % Prepare the function wrapper
            fcn  = @(x) FitObject.lossFunction(x, data, basisSet, baselineBasis, ppm, t, fitRange);
            grad = @(x) FitObject.forwardGradient(x, data, basisSet, baselineBasis, ppm, t, fitRange);
            fun  = @(x) FitObject.fminunc_wrapper(x, fcn, grad);
            % Request very high accuracy:
            opts            = struct('factr', 1e4, 'pgtol', 1e-5, 'm', 10);
            opts.printEvery = 0;
            % Run the algorithm:
            % Feed initial guess from the input parameters
            opts.x0 = x0;
            tstart = tic;
            [xk, ~, info] = lbfgsb(fun, lb, ub, opts );
            time = toc(tstart);
            

            

            % Save modeling results
            nBasisFcts      = size(basisSet.fids, 2);
            parsOut         = FitObject.x2pars(xk, nBasisFcts);
            
            % Save data (with the reverse phasing applied for visualization
            fdData = fftshift(fft(data, [], 1), 1);
            [fit, baseline, metabs] = FitObject.forwardModel(xk, basisSet, baselineBasis, ppm, t);
            ph0 = parsOut.ph0;
            ph1 = parsOut.ph1;
            T_phReverse = exp(1j .* (ph0 + ph1.*ppm)).';
            T_phReverse_mets = repmat(T_phReverse, [1 nBasisFcts]);
            
            dataReversePhased       = T_phReverse .* fdData;
            fitReversePhased        = T_phReverse .* fit;
            baselineReversePhased   = T_phReverse .* baseline.';
            metabsReversePhased     = T_phReverse_mets .* metabs;
            
            
            obj.Model.fit.data      = dataReversePhased;
            obj.Model.fit.fit       = fitReversePhased;
            obj.Model.fit.baseline  = baselineReversePhased;
            obj.Model.fit.residual  = dataReversePhased - fitReversePhased;
            obj.Model.fit.metabs    = metabsReversePhased;
            obj.Model.time          = time;
            obj.Model.parsOut       = parsOut;
            obj.Model.info          = info;
            
            % Calculate CRLB
            [jac, ~] = FitObject.forwardGradient(xk, data, basisSet, baselineBasis, ppm, t, fitRange);
            jac = jac.';
            
            % estimating the sigma based on the residual
            [indMin, indMax] = FitObject.ppmToIndex(ppm, fitRange);
            sigm   = std(obj.Model.fit.residual(indMin:indMax));
            
            % find zero indices
            zeroInd = find(~jac);
            
%             % OPTION 1 - get rid of zero elements
%             jac = nonzeros(jac).';
%             fisher = (1./(sigm^2)) .* jac.'*jac;
%             invFisher = pinv(fisher);
%             crlbs = diag(sqrtm(invFisher));
%                         
%             insertNumberInVector = @(a, x, n)cat(1,  x(1:n), a, x(n+1:end));
%             
%             for pp = 1:length(zeroInd)
%                 crlbs = insertNumberInVector(0, crlbs, zeroInd);
%             end
            
            
            % OPTION 2 - replace zeros with epsilon
            jac(zeroInd) = eps;
            jac = real(jac);
            fisher = jac.'*jac;
            invFisher = pinv(fisher);
            crlbs = diag(invFisher);
            
            CRLB = FitObject.x2pars(crlbs, nBasisFcts);
            
            %Relative CRLBs in percent
            relativeCRLB = 100 * CRLB.metAmpl./parsOut.metAmpl ./1024;

            % Save table with basis function names and relative CRLB
            % Only use basis functions that are included
            obj.Model.CRLB = table(basisSet.names(:, logical(basisSet.includeInFit))', relativeCRLB);

        end
        
        
        function [init] = getStartingValuesFromInitialFit(obj, init, parameter)
            
            % Collect the number of basis functions included in the fit
            nBasisFcts = sum(obj.BasisSets.includeInFit);
            
            if strcmp(obj.Options.parametrizations.(parameter).fun, 'free')
                switch parameter
                    case {'ph0', 'ph1', 'gaussLB', 'metAmpl', 'baseAmpl'}
                        % these were all determined during the initial fit
                        init.(parameter) = obj.Model.init.parsOut.(parameter)';
                    case {'freqShift', 'lorentzLB'}
                        % duplicate the initial value by the number of
                        % metabolite basis functions
                        init.(parameter) = repmat(obj.Model.init.parsOut.(parameter), [1, nBasisFcts]);
                end
                
            else
                error('Continue coding here when you start to use a parametrization');
            end
            
        end
        
        
        function [init, lb, ub] = initializeParameters(obj, init, lb, ub, parameter)
            
            % Collect number of basis functions and baseline spline elements
            nBasisFcts = sum(obj.BasisSets.includeInFit);
            nBaselineComps = size(obj.BaselineBasis, 2);
            
            if strcmp(obj.Options.parametrizations.(parameter).fun, 'free')
                switch parameter
                    case {'ph0', 'ph1', 'gaussLB'}
                        % one parameter per spectrum
                        nParamsPerSpec = 1;
                    case {'metAmpl', 'freqShift', 'lorentzLB'}
                        % one parameter per metabolite per spectrum
                        nParamsPerSpec = nBasisFcts;
                    case {'baseAmpl'}
                        % one parameter per spline function per spectrum
                        nParamsPerSpec = nBaselineComps;
                end
                init.(parameter) = repmat(obj.Options.parametrizations.(parameter).init, [size(obj.Data.fids,2), nParamsPerSpec]);
                lb.(parameter)   = repmat(obj.Options.parametrizations.(parameter).lb,   [size(obj.Data.fids,2), nParamsPerSpec]);
                ub.(parameter)   = repmat(obj.Options.parametrizations.(parameter).ub,   [size(obj.Data.fids,2), nParamsPerSpec]);
            else
                error('Continue coding here when you start to use a parametrization');
            end
            
        end
        
        
        
        % Visualization stuff below
        
        % Plot a spectrum
        function plotSpectra(obj, plotRange)
            % default to the provided fit range
            if nargin < 2
                plotRange = obj.Options.optimFreqFitRange;
            end
                
            figure;
            plot(obj.Data.ppm, real(fftshift(fft(obj.Data.fids,[],1),1)), 'k');
            set(gca, 'XDir', 'reverse', 'XLim', plotRange);
            
        end
        
        % Plot a complete fit
        function plotFit(obj, plotRange)
            % default to the provided fit range
            if nargin < 2
                plotRange = obj.Options.optimFreqFitRange;
            end
                
            % calculate plots
            ppm         = obj.Data.ppm;
            data        = obj.Model.fit.data;
            fit         = obj.Model.fit.fit;
            residual    = obj.Model.fit.residual;
            baseline    = obj.Model.fit.baseline;
            metabs      = obj.Model.fit.metabs;
            figure;
            hold on;
            plot(ppm, real(data), 'k');
            plot(ppm, real(fit), 'r');
            plot(ppm, real(residual) + max(real(data)), 'k');
            for rr = 1:size(metabs,2)
                plot(ppm, real(metabs(:,rr) + baseline), 'g');
            end
            plot(ppm, real(baseline'), 'b');
            hold off;
            
            set(gca, 'XDir', 'reverse', 'XLim', plotRange);
            xlabel('chemical shift (ppm');
            
        end
        
        % Plot an initial fit
        function plotInitFit(obj, plotRange)
            % default to the provided fit range
            if nargin < 2
                plotRange = obj.Options.optimFreqFitRange;
            end
                
            % calculate plots
            ppm         = obj.Data.ppm;
            data        = obj.Model.init.fit.data;
            fit         = obj.Model.init.fit.fit;
            residual    = obj.Model.init.fit.residual;
            baseline    = obj.Model.init.fit.baseline;
            metabs      = obj.Model.init.fit.metabs;
            
            figure;
            hold on;
            plot(ppm, real(data), 'k');
            plot(ppm, real(fit), 'r');
            plot(ppm, real(residual) + max(real(data)), 'k');
            for rr = 1:size(metabs,2)
                plot(ppm, real(metabs(:,rr) + baseline), 'g');
            end
            plot(ppm, real(baseline), 'b');
            hold off;
            
            set(gca, 'XDir', 'reverse', 'XLim', plotRange);
            xlabel('chemical shift (ppm');
            
        end
        
        % Plot the basis set
        function plotBasisSet(obj, plotRange)
            
            % default to the provided fit range
            if nargin < 2
                plotRange = obj.Options.optimFreqFitRange;
            end
                
            % fft time-domain basis set
            fdBasisSpecs = real(fftshift(fft(obj.BasisSets.fids,[],1),1));
            % calculate stagger
            stag = median(max(abs(real(fdBasisSpecs))));
            % calculate text position
            ppm = obj.Data.ppm;
            xText = round((max(plotRange) - min(plotRange))*0.1 + min(plotRange));
                        
            figure;
            hold on;

            for rr = 1:length(obj.BasisSets.names)
                % plot basis functions not included in the fit in red
                if obj.BasisSets.includeInFit(rr) == 1
                    colorLine = 'k';
                else
                    colorLine = 'r';
                end
                plot(ppm, fdBasisSpecs(:,rr) - stag*(rr-1), colorLine);
                text(xText, -stag*(rr-1)+0.5*stag, obj.BasisSets.names{rr}, 'Color', colorLine, 'Clipping', 'on');
                
                
            end
            hold off;
            set(gca, 'XDir', 'reverse', 'XLim', plotRange);
            xlabel('chemical shift (ppm');
            
        end
        
        % Plot the baseline basis
        function plotBaselineBasis(obj, plotRange)
            
            % default to the provided fit range
            if nargin < 2
                plotRange = obj.Options.optimFreqFitRange;
            end
                
            figure;
            plot(obj.Data.ppm, real(obj.BaselineBasis), 'k');
            set(gca, 'XDir', 'reverse', 'XLim', plotRange);
            xlabel('chemical shift (ppm');
            
        end
        
        
        % Rescale a basis set function
        function rescaleBasisSetFunction(obj, name, factor)
            
            if nargin < 3
                rescaleToMean = 1;
            else
                rescaleToMean = 0;
            end
            
            idxName = find(strcmpi(name, obj.BasisSets.names));
            if isempty(idxName)
                error('Could not find basis function %s in the basis set.', name);
            else
                if rescaleToMean
                    % Calculate the norm of the basis function to be scaled
                    currentNorm = norm(obj.BasisSets.fids(:,idxName));
                    % Calculate the mean of the norm of all other basis
                    % functions
                    meanNorm = mean(norm(obj.BasisSets.fids(:,~strcmpi(name, obj.BasisSets.names))));
                    % Calculate conversion factor
                    factor = meanNorm/currentNorm;
                else
                    % do nothing
                end
                obj.BasisSets.fids(:,idxName) = obj.BasisSets.fids(:,idxName) .* factor;
            end
            
        end
        
        
        % Rescale an entire basis set
        % This helps with cases where the input basis set is just so poorly
        % scaled that it messes with the initial fit
        function rescaleBasisSetToOne(obj, factor)
            
            % If no input is provided, the basis set is scaled such that
            % the maximum in the frequency domain is one
            if nargin < 2
                rescaleToOne = 1;
            else
                rescaleToOne = 0;
            end
            
            % Find maximum of the basis set array (frequency domain)
            fdBasisSpecs = fftshift(fft(obj.BasisSets.fids,[],1),1);
            currentMax = max(real(fdBasisSpecs), [], 'all');
            
            if rescaleToOne
                factor = 1/currentMax;
            end
            
            obj.BasisSets.fids = obj.BasisSets.fids .* factor;
            
        end
        
        % Removes TMS peak from the basis set
        function removeTMSFromBasisSet(obj)
            
            % If no input is provided, the basis set is scaled such that
            % the maximum in the frequency domain is one
            if nargin < 2
                rescaleToOne = 1;
            else
                rescaleToOne = 0;
            end
            
            % Find maximum of the basis set array (frequency domain)
            fdBasisSpecs = fftshift(fft(obj.BasisSets.fids,[],1),1);
            currentMax = max(real(fdBasisSpecs), [], 'all');
            
            if rescaleToOne
                factor = 1/currentMax;
            end
            
            obj.BasisSets.fids = obj.BasisSets.fids .* factor;
            
        end
            

            
    end
    
    
    % Static methods, helper functions
    methods (Static)
        
        % Forward models, loss functions, Jacobian and gradient functions
        % below
        function sse = lossFunction(x, data, basisSet, baselineBasis, ppm, t, fitRange)
            
            [indMin, indMax] = FitObject.ppmToIndex(ppm, fitRange);
            
            % apply the forward model
            data        = fftshift(fft(data, [], 1),1);
            prediction  = FitObject.forwardModel(x, basisSet, baselineBasis, ppm, t);
            diffVec     = data(indMin:indMax) - prediction(indMin:indMax);
            
            % OPTION 1
            sse         = real(sum(diffVec .* conj(diffVec)));
            
            % OPTION 2
            % sse         = sum(real(diffVec).^2);
%             % normalize by 1/(sd(noise))^2
%             noiseEstFrac = round(0.8*size(data,1));
%             noise       = detrend(real(data(noiseEstFrac:end)));
%             sse         = sse .* 1/(std(noise).^2);
            
            
        end
        
        function [Y, baseline, metabs] = forwardModel(x, basisSet, baselineBasis, ppm, t)
            
            % Define the default 1-D forward model first
            fids = basisSet.fids;
            nBasisFcts = sum(basisSet.includeInFit);
            
            inputParams     = FitObject.x2pars(x, nBasisFcts);
            gaussLB     = inputParams.gaussLB;
            lorentzLB   = inputParams.lorentzLB;
            freqShift   = inputParams.freqShift;
            metAmpl     = inputParams.metAmpl;
            baseAmpl    = inputParams.baseAmpl;
            ph0         = inputParams.ph0;
            ph1         = inputParams.ph1;
            
            timeDomainMultiplier = zeros(size(fids));
            for ll = 1:nBasisFcts
                timeDomainMultiplier(:,ll) = exp(-(1i*freqShift(ll) + lorentzLB(ll) + gaussLB.^2.*t).*t).';    
            end
            
            Fl = timeDomainMultiplier .* fids;
            specs = fftshift(fft(Fl, [], 1),1);
            
            mets = metAmpl' * specs.';
            baseline = baseAmpl' * baselineBasis.';
            
            T_ph = exp(-1j .* (ph0 + ph1.*ppm));
            Y = T_ph .* (mets + baseline);
            Y = Y.';
            
            metabs = repmat(T_ph.', [1 nBasisFcts]) .* specs .* metAmpl';
            baseline = T_ph .* baseline;
        end
        

        
        function [grad, jac] = forwardGradient(x, data, basisSet, baselineBasis, ppm, t, fitRange)
            
            % get fit range
            [indMin, indMax] = FitObject.ppmToIndex(ppm, fitRange);
            
            % get prediction
            prediction  = FitObject.forwardModel(x, basisSet, baselineBasis, ppm, t);
            
            % Construct the Jacobian matrix of partial derivatives
            fids = basisSet.fids;
            nBasisFcts = size(fids,2);
            nBaselineComps = size(baselineBasis, 2);
            
            inputParams = FitObject.x2pars(x, nBasisFcts);
            gaussLB = inputParams.gaussLB;
            lorentzLB = inputParams.lorentzLB;
            freqShift = inputParams.freqShift;
            metAmpl = inputParams.metAmpl;
            baseAmpl = inputParams.baseAmpl;
            ph0 = inputParams.ph0;
            ph1 = inputParams.ph1;
            
            % Apply the time-domain parameters (Gauss, Lorentz, frequency
            % shift)
            tdMultiplier = zeros(size(fids));
            for ll = 1:nBasisFcts
                tdMultiplier(:,ll) = exp(-(1j*freqShift(ll) + lorentzLB(ll) + gaussLB.^2.*t).*t).';    
            end
            tdTransformedBasis = tdMultiplier .* fids;
            
            % Create a couple of useful vectors and matrices containing the
            % phase term. (We rather want to do the repmat here to avoid
            % making the terms even uglier later)
            T_ph = exp(-1j .* (ph0 + ph1.*ppm).');
            T_ph_basis      = repmat(T_ph, [1, nBasisFcts]);
            T_ph_baseline   = repmat(T_ph, [1, nBaselineComps]);
            
            % Create matrices containing the time vector and its square.
            % We need these as time appears in the derivative of the
            % time-domain signal with respect to the time-domain parameters
            T_t     = repmat(t.', [1, nBasisFcts]);
            T_tt    = repmat(t.'.*t.', [1, nBasisFcts]);
            
            % These are derivatives according to the chain rule
            fdtBasis = fftshift(fft(T_t .* tdTransformedBasis, [], 1),1);
            fdttBasis = fftshift(fft(T_tt .* tdTransformedBasis, [], 1),1);
            
            % partial derivatives
            dYdph0          = (-1j) .* prediction;
            dYdph1          = (-1j) .* ppm.' .* prediction; 
            dYdgaussLB      = -2 .* gaussLB .* T_ph .* fdttBasis * metAmpl; 
            dYdlorentzLB    = T_ph_basis .* (-1)  .* fdtBasis .* metAmpl';
            dYdfreqShift    = T_ph_basis .* (-1j) .* fdtBasis .* metAmpl';
            dYdmetAmpl      = T_ph_basis .* fftshift(fft(tdTransformedBasis, [], 1),1);
            dYdbaseAmpl     = T_ph_baseline .* baselineBasis;
            
            % reduce to fit range
            dYdph0          = dYdph0(indMin:indMax,:);
            dYdph1          = dYdph1(indMin:indMax,:);
            dYdgaussLB      = dYdgaussLB(indMin:indMax,:);
            dYdlorentzLB    = dYdlorentzLB(indMin:indMax,:);
            dYdfreqShift    = dYdfreqShift(indMin:indMax,:);
            dYdmetAmpl      = dYdmetAmpl(indMin:indMax,:);
            dYdbaseAmpl     = dYdbaseAmpl(indMin:indMax,:);
            
            % reduce data and prediction to fit range;
            prediction  = prediction(indMin:indMax);
            data        = fftshift(fft(data,[],1),1);
            data        = data(indMin:indMax);
            
            
            % concatenate vertically to create the matrix of partial
            % derivatives for each data point (the Jacobian of the forward
            % model)
            jac = cat(2, dYdph0, dYdph1, dYdgaussLB, dYdlorentzLB, dYdfreqShift, dYdmetAmpl, dYdbaseAmpl);
            
%             % normalize by 1/(sd(noise))^2
%             noiseEstFrac = round(0.8*size(data,1));
%             noise       = detrend(real(data(noiseEstFrac:end)));
%             jac         = jac .* 1/(std(noise).^2);
            
            % OPTION 1
            % apply the product rule to the loss function:
            % let u = (data-pred) and v = conj(data-pred)
            % d(uv)/dx = vdu/dx + udv/dx            (ddata/dx = 0)
            %          = v(-dpred/dx) + u(-conj(dpred/dx)
            %          = conj(data-pred)*(-dpred/dx) +
            %          (data-pred)*(-conj*dpred/dx)
            grad = real(sum((data - prediction).*(-conj(jac)) + (-jac .* conj(data-prediction)),1)).';
            
            % OPTION 2
            % apply the chain rule to the loss function:
            %grad = real(sum(2.*(-jac).*(data - prediction))).';

        end
        
        
        function sse = initialFitLossFunction(x, tdData, basisSet, baselineBasis, ppm, t, fitRange)
            
            % get time-domain data, fft, define fit range
            fdData = fftshift(fft(tdData,[],1),1);
            [indMin, indMax] = FitObject.ppmToIndex(ppm, fitRange);
            
            % get basis sets, find number of basis functions
            fids = basisSet.fids;
            nBasisFcts = size(fids,2);
            
            % apply non-linear parameters to basis set and fft
            ph0         = x(1);
            gaussLB     = exp(x(2));
            lorentzLB   = exp(x(3));
            freqShift   = x(4);
            specs       = FitObject.transformBasis(fids, gaussLB, lorentzLB, freqShift, t);
            
            % append basis set with the spline baseline basis functions
            fullBasis = cat(2, specs, baselineBasis);
            % apply zero-order phase
            fullBasis = exp(-1j*ph0) .* fullBasis;
            
            % crop to fit range
            fdData = fdData(indMin:indMax,:);
            fullBasis = fullBasis(indMin:indMax,:);
            
            % concatenate real and imaginary
            fdData = cat(1, real(fdData), imag(fdData));
            fullBasis = cat(1, real(fullBasis), imag(fullBasis));
            
            % solve the linear equation system
            beta = real(pinv(fullBasis) * fdData);
            
            % project metabolite coefficients to >0
            beta(beta(1:nBasisFcts) < 0) = 0;
            
            % make spectrum
            prediction = fullBasis * beta;
            
            % apply the forward model
            diffVec     = fdData - prediction;
            sse         = mean(abs(diffVec).*2);
            
        end
        
        function specs = transformBasis(fids, gaussLB, lorentzLB, freqShift, t)
            timeDomainMultiplier = zeros(size(fids));
            nBasisFcts = size(fids,2);
            for ll = 1:nBasisFcts
                timeDomainMultiplier(:,ll) = exp(-(1i*freqShift + lorentzLB + gaussLB.^2.*t).*t).';
            end
            
            Fl = timeDomainMultiplier .* fids;
            specs = fftshift(fft(Fl, [], 1), 1);
        end
            
        function [f,g,h] = fminunc_wrapper(x,F,G,H)
            % [f,g,h] = fminunc_wrapper( x, F, G, H )
            % for use with Matlab's "fminunc"
            f = F(x);
            if nargin > 2 && nargout > 1
                g = G(x);
            end
            if nargin > 3 && nargout > 2
                h = H(x);
            end
        end
        
        function [indMin, indMax] = ppmToIndex(ppmVector, ppmRange)
            
            % Returns indices of the ppm vector corresponding to a given
            % 2-element range
            [~,indMin] = min(abs(ppmVector - ppmRange(1)));
            [~,indMax] = min(abs(ppmVector - ppmRange(2)));
            
        end
        
        function x = pars2x(paramStruct)
            
            % Converts a parameter struct into a 1-D x vector that can be
            % passed on to solvers
            x = [paramStruct.ph0, paramStruct.ph1, paramStruct.gaussLB, ...
                 paramStruct.lorentzLB, paramStruct.freqShift, ... 
                 paramStruct.metAmpl, paramStruct.baseAmpl]';
            
        end
        
        function paramStruct = x2pars(x, nBasisFcts)
            
            % Converts a 1-D x vector into a parameter struct
            paramStruct.ph0 = x(1);
            paramStruct.ph1 = x(2);
            paramStruct.gaussLB = x(3);
            paramStruct.lorentzLB = x(4:3+nBasisFcts);
            paramStruct.freqShift = x(4+nBasisFcts:3+2*nBasisFcts);
            paramStruct.metAmpl = x(4+2*nBasisFcts:3+3*nBasisFcts);
            paramStruct.baseAmpl = x(4+3*nBasisFcts:end);

            
        end
        
        function CRLB = getDiagonalElementsOfCRLB(x, nBasisFcts)
            
            % Converts a 1-D x vector into a parameter struct
            CRLB.ph0 = diag(x(1,1));
            CRLB.ph1 = diag(x(2,2));
            CRLB.gaussLB = diag(x(3,3));
            CRLB.lorentzLB = diag(x(4:3+nBasisFcts, 4:3+nBasisFcts));
            CRLB.freqShift = diag(x(4+nBasisFcts:3+2*nBasisFcts, 4+nBasisFcts:3+2*nBasisFcts));
            CRLB.metAmpl = diag(x(4+2*nBasisFcts:3+3*nBasisFcts, 4+2*nBasisFcts:3+3*nBasisFcts));
            CRLB.baseAmpl = diag(x(4+3*nBasisFcts:end, 4+3*nBasisFcts:end));
            
            
        end
        
        
        
    end
    

end