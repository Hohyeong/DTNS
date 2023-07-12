function slBusOut = NavSatStatus(msgIn, slBusOut, varargin)
%#codegen
%   Copyright 2021 The MathWorks, Inc.
    slBusOut.Status = int8(msgIn.Status);
    slBusOut.Service = uint16(msgIn.Service);
end
