%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEM
% methods:  DEM(DEMData, DBinfo)
%           GetElevation(lat, lon)
%           GetIndexFromLatLon(obj, lat, lon)
% * Maintainer : Sungjoong Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef DEM < handle
    properties
        NCols
        NRows
        MinX
        MaxX
        MinY
        MaxY
        DB
        ResX
        ResY
        MaxH
    end
    methods
        function obj = DEM(DEMData, DBinfo)
            % DEMData : Matrix 형태의 원본 DB data
            % DBinfo  : inforDB 구조체(혹은 inforMap)
            % MapOrDB : DB or Map flag
            
            obj.DB = DEMData;
            obj.NCols = size(DEMData,2); % 3601 for 1-deg DTED2 data;
            obj.NRows = size(DEMData,1); % 3601 for 1-deg DTED2 data;
            obj.MaxH = max(max(DEMData));
            
            obj.ResX = DBinfo.GridDegDistMap; % = 1/3600 for DTED2;
            obj.ResY = DBinfo.GridDegDistMap; % = 1/3600 for DTED2;
            
            obj.MinX = DBinfo.initGridDegLonMap; % 125.0; [125, 130);
            obj.MaxY = DBinfo.initGridDegLatMap; % 40.0; (35, 40];
            obj.MaxX = obj.MinX + obj.ResX*(obj.NCols-1);
            obj.MinY = obj.MaxY - obj.ResY*(obj.NRows-1);

        end
        
        function height = GetElevation(obj, lat, lon)
            % Query terrain elevation. 
            % Use inerpolation of adjacent terrain elevation data.
            % Args:
            % obj       | self
            % lat       | latitude [deg] 
            % lon       | longitude [deg]
            % Returns:
            % height    | terrain elevation at query point [m]
            % Example:
            % DB.GetElevation(latitude(deg), longitude(deg));
            
            % Get smllest index of box that contains query point
            [idxLat, idxLon] = GetIndexFromLatLon(obj, lat, lon);
            
            % Remainder length
            lenLon  = lon - (obj.MinX + idxLon * obj.ResX);
            lenLat  = lat - (obj.MinY + idxLat * obj.ResY);
            
            % Weighting factor
            f       = lenLon / obj.ResX;
            g       = lenLat / obj.ResY;
            
            % Terrain elevation at corners
            h11     = obj.DB(idxLat  , idxLon  );
            h21     = obj.DB(idxLat+1, idxLon  );
            h12     = obj.DB(idxLat  , idxLon+1);
            h22     = obj.DB(idxLat+1, idxLon+1);
            
            hLow    = (1-f) * double(h11) + f * double(h12);
            hHigh   = (1-f) * double(h21) + f * double(h22);
            
            height  = (1-g) * hLow + g * hHigh;            
                        
        end
        
        function [idxLat, idxLon] = GetIndexFromLatLon(obj, lat, lon)
            % Query the smallest lat/lon indices of the grid which contains
            % Query point
            % Args:
            % obj       | self
            % lat       | latitude [deg] 
            % lon       | longitude [deg]
            % Returns:
            % idxLat    | smallest lat index of the grid
            % idxLon    | smallest lon index of the grid
            
            idxLon  = fix( (lon - obj.MinX) / obj.ResX ) + 1;
            idxLat  = fix( (lat - obj.MinY) / obj.ResY ) + 1;

        end
        
        
        function height = CalcHeightBL(obj, Lon, Lat)
            % DB.CalcHeightBL(longitude(deg), latitude(deg)); 로 호출
            % 사용금지 !!!!!
                      
            if obj.MinX <= Lon && Lon < obj.MaxX && obj.MinY < Lat && Lat <= obj.MaxY
                IndexX = fix((Lon - obj.MinX) / obj.ResX);
                IndexY = fix((obj.MaxY - Lat) / obj.ResY);

%                 LengthX = Lon - (obj.MinX + IndexX     * obj.ResX);
%                 LengthY = Lat - (obj.MaxY - (IndexY+1) * obj.ResY);

%                 f = LengthX / obj.ResX;
%                 g = LengthY / obj.ResY;
% 
%                 H1 = (obj.DB(IndexX+1, IndexY+2));
%                 H2 = (obj.DB(IndexX+1, IndexY+1));
%                 H3 = (obj.DB(IndexX+2, IndexY+1));
%                 H4 = (obj.DB(IndexX+2, IndexY+2));
% 
%                 R1 = (1 - f) * double(H1) + f * double(H4);
%                 R2 = (1 - f) * double(H2) + f * double(H3);
% 
%                 height = (1 - g) * R1 + g * R2;
                
                LengthX = (Lon - obj.MinX) - IndexX * obj.ResX;
                LengthY = (obj.MaxY - Lat) - IndexY * obj.ResY;
                
                f = LengthX / obj.ResX;
                g = LengthY / obj.ResY;

                % Matlab call convention : M(idx_y{row}, idx_x{column})
                H1 = (obj.DB(IndexY+1, IndexX+1));
                H2 = (obj.DB(IndexY+2, IndexX+1));
                H3 = (obj.DB(IndexY+2, IndexX+2));
                H4 = (obj.DB(IndexY+1, IndexX+2));

                R1 = (1 - f) * double(H1) + f * double(H4);
                R2 = (1 - f) * double(H2) + f * double(H3);

                height = (1 - g) * R1 + g * R2;
                
            else
                % Query coordinate is out of data bound
                warning('Elevation query at ( Lat : %.5f deg, Lon : %.5f deg) is out of database bound', Lat, Lon);
                height = -1;
            end
        end
    end
end