% op_phaseFirstPoint.m
% Georg Oeltzschner, Johns Hopkins University 2021
% 
% USAGE:
% [out] = op_phaseFirstPoint(in);
% 
% DESCRIPTION:
% Phases all FIDs by the first point.
% 
% INPUTS:
% in        = Input data structure. 
%
% OUTPUTS:
% out       = Output data structure


function [out] = op_phaseFirstPoint(in)

if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

% find the first points
firstPoints = in.fids(1,:,:,:);

% extend along t dimension
conjFirstPoints = repmat(conj(firstPoints), [in.sz(1), 1]);

% apply 
fids = in.fids .* conjFirstPoints;
    
%re-calculate Specs using fft
specs=fftshift(fft(fids,[],in.dims.t),in.dims.t);


%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
    
end
