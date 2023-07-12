function slBusOut = HomePosition(msgIn, slBusOut, varargin)
%#codegen
%   Copyright 2021 The MathWorks, Inc.
    currentlength = length(slBusOut.Header);
    for iter=1:currentlength
        slBusOut.Header(iter) = bus_conv_fcns.ros.msgToBus.std_msgs.Header(msgIn.Header(iter),slBusOut(1).Header(iter),varargin{:});
    end
    slBusOut.Header = bus_conv_fcns.ros.msgToBus.std_msgs.Header(msgIn.Header,slBusOut(1).Header,varargin{:});
    currentlength = length(slBusOut.Geo);
    for iter=1:currentlength
        slBusOut.Geo(iter) = bus_conv_fcns.ros.msgToBus.geographic_msgs.GeoPoint(msgIn.Geo(iter),slBusOut(1).Geo(iter),varargin{:});
    end
    slBusOut.Geo = bus_conv_fcns.ros.msgToBus.geographic_msgs.GeoPoint(msgIn.Geo,slBusOut(1).Geo,varargin{:});
    currentlength = length(slBusOut.Position);
    for iter=1:currentlength
        slBusOut.Position(iter) = bus_conv_fcns.ros.msgToBus.geometry_msgs.Point(msgIn.Position(iter),slBusOut(1).Position(iter),varargin{:});
    end
    slBusOut.Position = bus_conv_fcns.ros.msgToBus.geometry_msgs.Point(msgIn.Position,slBusOut(1).Position,varargin{:});
    currentlength = length(slBusOut.Orientation);
    for iter=1:currentlength
        slBusOut.Orientation(iter) = bus_conv_fcns.ros.msgToBus.geometry_msgs.Quaternion(msgIn.Orientation(iter),slBusOut(1).Orientation(iter),varargin{:});
    end
    slBusOut.Orientation = bus_conv_fcns.ros.msgToBus.geometry_msgs.Quaternion(msgIn.Orientation,slBusOut(1).Orientation,varargin{:});
    currentlength = length(slBusOut.Approach);
    for iter=1:currentlength
        slBusOut.Approach(iter) = bus_conv_fcns.ros.msgToBus.geometry_msgs.Vector3(msgIn.Approach(iter),slBusOut(1).Approach(iter),varargin{:});
    end
    slBusOut.Approach = bus_conv_fcns.ros.msgToBus.geometry_msgs.Vector3(msgIn.Approach,slBusOut(1).Approach,varargin{:});
end
