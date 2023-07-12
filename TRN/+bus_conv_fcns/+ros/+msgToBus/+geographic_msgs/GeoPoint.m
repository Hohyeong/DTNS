function slBusOut = GeoPoint(msgIn, slBusOut, varargin)
%#codegen
%   Copyright 2021 The MathWorks, Inc.
    slBusOut.Latitude = double(msgIn.Latitude);
    slBusOut.Longitude = double(msgIn.Longitude);
    slBusOut.Altitude = double(msgIn.Altitude);
end
