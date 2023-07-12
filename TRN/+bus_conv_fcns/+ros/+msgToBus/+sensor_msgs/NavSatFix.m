function slBusOut = NavSatFix(msgIn, slBusOut, varargin)
%#codegen
%   Copyright 2021 The MathWorks, Inc.
    currentlength = length(slBusOut.Header);
    for iter=1:currentlength
        slBusOut.Header(iter) = bus_conv_fcns.ros.msgToBus.std_msgs.Header(msgIn.Header(iter),slBusOut(1).Header(iter),varargin{:});
    end
    slBusOut.Header = bus_conv_fcns.ros.msgToBus.std_msgs.Header(msgIn.Header,slBusOut(1).Header,varargin{:});
    currentlength = length(slBusOut.Status);
    for iter=1:currentlength
        slBusOut.Status(iter) = bus_conv_fcns.ros.msgToBus.sensor_msgs.NavSatStatus(msgIn.Status(iter),slBusOut(1).Status(iter),varargin{:});
    end
    slBusOut.Status = bus_conv_fcns.ros.msgToBus.sensor_msgs.NavSatStatus(msgIn.Status,slBusOut(1).Status,varargin{:});
    slBusOut.Latitude = double(msgIn.Latitude);
    slBusOut.Longitude = double(msgIn.Longitude);
    slBusOut.Altitude = double(msgIn.Altitude);
                    currentlength = length(slBusOut.PositionCovariance);
                    slBusOut.PositionCovariance = double(msgIn.PositionCovariance(1:currentlength));
    slBusOut.PositionCovarianceType = uint8(msgIn.PositionCovarianceType);
end
