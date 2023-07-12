// Zone 1 영역 지형고도 조회용
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <io.h>
#include "MAP.h"

#ifdef _WIN32
#include <WinSock2.h>
#elif
#include <arpa/inet.h>
#endif

char fileName[BUF_SIZE];
MAPDATA	_MapData[2];
MAP_POST_HANDLE _map_postHandle[5];
unsigned int _map_postHandleFront = 0;

//------------------------------------------------------------------------------------------------------
// Function : InitialMap
//		지형 FILE에 저장된 디렉토리 경로를 설정
//
// Arguments : 
//		char* path  - 디렉토리 경로
//		int size - 문자열 길이
//------------------------------------------------------------------------------------------------------
int InitialMap(char* path, int size)

{	strncpy(_MapData[0].mapPath, path, size);
	_MapData[0].mapPath[size] = '\0';

	strncpy(_MapData[1].mapPath, path, size);
	_MapData[1].mapPath[size] = '\0';

	_map_postHandleFront = 0;
	
	return _access(_MapData[0].mapPath, 0);
}

//------------------------------------------------------------------------------------------------------
// Function : UpdateMapCache
//		지형 및 DVOF Cache를 업데이트 한다.
//
// Arguments : 
//		double lat  - 위도
//		double lon  - 경도
//------------------------------------------------------------------------------------------------------
int UpdateCache(double lat, double lon)

{	int result  = 0;
	if (IsUpdateBoundary(0, lat, lon) == 0)
	
{		result |= UpdateMapCache(lat, lon);
	}

	if (IsUpdateBoundary(1, lat, lon) == 0)
	
{		result |= UpdateObsCache(lat, lon);
	}

	return result;
}

//------------------------------------------------------------------------------------------------------
// Function : UpdateMapCache
//		지형 Cache를 업데이트 한다.
//
// Arguments : 
//		double lat  - 위도
//		double lon  - 경도
//------------------------------------------------------------------------------------------------------
int UpdateMapCache(double lat, double lon)

{	int latDeg = (int)(lat / CELL_INTERVAL) / 5;
	int latMinute = (int)((int)(lat*100) / (int)(CELL_INTERVAL*100)) % 5 * (int)(CELL_INTERVAL * 60);
	int lonDeg = (int)(lon / CELL_INTERVAL) / 5;
	int lonMinute = (int)((int)(lon * 100) / (int)(CELL_INTERVAL * 100)) % 5 * (int)(CELL_INTERVAL * 60);

	ReadMaps(lat, lon);

	double lat2 = latDeg + latMinute*(1.0 / 60.0);
	double lon2 = lonDeg + lonMinute*(1.0 / 60.0);
	SetUpdateBoundary(0, lat2, lon2);
	SetScanBoundary(0);

	_MapData[0].mapCache.resolution = RESOLUTION;

	return SUCCESS;
}

int UpdateObsCache(double lat, double lon)

{	int latDeg = (int)(lat / CELL_INTERVAL) / 5;
	int latMinute = (int)((int)(lat * 100) / (int)(CELL_INTERVAL * 100)) % 5 * (int)(CELL_INTERVAL * 60);
	int lonDeg = (int)(lon / CELL_INTERVAL) / 5;
	int lonMinute = (int)((int)(lon * 100) / (int)(CELL_INTERVAL * 100)) % 5 * (int)(CELL_INTERVAL * 60);
	
	// clean the obstacle
	memset(_MapData[1].mapCache.map, 0x00, sizeof(_MapData[1].mapCache.map));
	
	double lat2 = latDeg + latMinute*(1.0 / 60.0);
	double lon2 = lonDeg + lonMinute*(1.0 / 60.0);
	SetUpdateBoundary(1, lat2, lon2);
	SetScanBoundary(1);

	_MapData[1].mapCache.resolution = RESOLUTION;

	return SUCCESS;
}


int ReadMaps(double lat, double lon)

{	int result = SUCCESS;

	// 부동소수점 오류로 소수점 8자리에서 반올림
	double roundLat = lat + 0.00000001;
	double roundlon = lon + 0.00000001;

	int baseLatDeg = (int)(roundLat / CELL_INTERVAL) / DEVIDE_CNT;
	int baseLatMin = (int)floor((roundLat - ((double)baseLatDeg)) * 60 / (double)CELLUNIT) * CELLUNIT;

	int baseLonDeg = (int)(roundlon / CELL_INTERVAL) / DEVIDE_CNT;
	int baseLonMin = (int)floor((roundlon - ((double)baseLonDeg)) * 60 / (double)CELLUNIT) * CELLUNIT;

	int mid = CELLCNT / 2;
	for (int y = -mid; y <= mid; y++)
	
{		for (int x = -mid; x <= mid; x++)
		
{			// Longitude
			int lonDeg = baseLonDeg;
			int lonMin = baseLonMin + (x * CELLUNIT);
			lonDeg += lonMin / 60;
			lonMin %= 60;

			if (lonDeg > 0 && lonMin < 0)
			
{				lonDeg--;
				lonMin += 60;
			}
			else if (lonDeg < 0 && lonMin > 0)
			
{				lonDeg++;
				lonMin -= 60;
			}

			// Latitude
			int latDeg = baseLatDeg;
			int latMin = baseLatMin + (y * CELLUNIT);
			latDeg += latMin / 60;
			latMin %= 60;

			if (latDeg > 0 && latMin < 0)
			
{				latDeg--;
				latMin += 60;
			}
			else if (latDeg < 0 && latMin > 0)
			
{				latDeg++;
				latMin -= 60;
			}
		

			// Read a Sheet
			if (ReadCell(x + mid, y + mid, latDeg, latMin, lonDeg, lonMin) == FAIL)
{				result = FAIL;
			}
		}
	}

	return result;
}


int ReadCell(int x, int y, int latDeg, int latMin, int lonDeg, int lonMin)


{	sprintf(fileName, "%s/%c%02d_%02d_%c%03d_%02d.map",
		_MapData[0].mapPath,
		(latDeg >= 0.0) ? 'N' : 'S',
		abs(latDeg),
		abs(latMin),
		(lonDeg >= 0.0) ? 'E' : 'W',
		abs(lonDeg),
		abs(lonMin));

	printf("%s\n", fileName);

	//  Read file
	FILE* fd = fopen(fileName, "rb");
	if (fd == NULL)
	

{		// 맵 오프인 경우, 알고리즘에서 MAP OFF 판단이 가능하도록
		// 해당 고도 데이터를 -1로 설정
		ProcessMapOff(0, x, y);

		return NOT_EXIST_MAP;
	}

	MAP_HEADER header;
	fread(&header, 1, sizeof(header), fd);
	for (int i = 0; i < 7; i++)
	
{		header.raw[i] = ntohs(header.raw[i]);
	}

	short temp[MAPCACHE_SIZE];
	for (int lonIdx = 0; lonIdx < header.data.numOfLon-1; lonIdx++)
	
{		int size = fread(temp, 1, header.data.numOfLat * sizeof(short), fd);
		for (int latIdx = 0; latIdx < header.data.numOfLat-1; latIdx++)
		
{			int lon_idx = x*CELLCACHE_SIZE + lonIdx;
			int lat_idx = y*CELLCACHE_SIZE + latIdx;
			//printf("%d, %d\n", lon_idx, lat_idx);
			_MapData[0].mapCache.map[lon_idx][lat_idx] = ntohs(temp[latIdx]);
		}
	}

	fclose(fd);
	return 0;

}


void ProcessMapOff(int type, int x, int y)
{	for (int lonIdx = 0; lonIdx < CELLCACHE_SIZE; lonIdx++)
	
{		for (int latIdx = 0; latIdx < CELLCACHE_SIZE; latIdx++)
		
{			int lon_idx = x * CELLCACHE_SIZE + lonIdx;
			int lat_idx = y * CELLCACHE_SIZE + latIdx;
			
			_MapData[type].mapCache.map[lon_idx][lat_idx] = -1.0;
		}
	}
}

//------------------------------------------------------------------------------------------------------
// Description :
//		
//------------------------------------------------------------------------------------------------------
void SetUpdateBoundary(int type, double lat, double lon)

{	
	_MapData[type].mapCache.updateBoundary.leftBottom.lat = lat;
	_MapData[type].mapCache.updateBoundary.leftBottom.lon = lon;

	_MapData[type].mapCache.updateBoundary.leftTop.lat = (lat + CELL_INTERVAL);
	_MapData[type].mapCache.updateBoundary.leftTop.lon = lon;

	_MapData[type].mapCache.updateBoundary.rightTop.lat = (lat + CELL_INTERVAL);
	_MapData[type].mapCache.updateBoundary.rightTop.lon = (lon + CELL_INTERVAL);

	_MapData[type].mapCache.updateBoundary.rightBottom.lat = lat;
	_MapData[type].mapCache.updateBoundary.rightBottom.lon = (lon + CELL_INTERVAL);
}


//------------------------------------------------------------------------------------------------------
// Description :
//		
//------------------------------------------------------------------------------------------------------
void SetScanBoundary(int type)

{	_MapData[type].mapCache.scanBoundary.leftBottom.lat = _MapData[type].mapCache.updateBoundary.leftBottom.lat - CELL_INTERVAL;
	_MapData[type].mapCache.scanBoundary.leftBottom.lon = _MapData[type].mapCache.updateBoundary.leftBottom.lon - CELL_INTERVAL;

	_MapData[type].mapCache.scanBoundary.leftTop.lat = _MapData[type].mapCache.updateBoundary.leftTop.lat + CELL_INTERVAL;
	_MapData[type].mapCache.scanBoundary.leftTop.lon = _MapData[type].mapCache.updateBoundary.leftTop.lon - CELL_INTERVAL;

	_MapData[type].mapCache.scanBoundary.rightTop.lat = _MapData[type].mapCache.updateBoundary.rightTop.lat + CELL_INTERVAL;
	_MapData[type].mapCache.scanBoundary.rightTop.lon = _MapData[type].mapCache.updateBoundary.rightTop.lon + CELL_INTERVAL;

	_MapData[type].mapCache.scanBoundary.rightBottom.lat = _MapData[type].mapCache.updateBoundary.rightBottom.lat - CELL_INTERVAL;
	_MapData[type].mapCache.scanBoundary.rightBottom.lon = _MapData[type].mapCache.updateBoundary.rightBottom.lon + CELL_INTERVAL;
}

//------------------------------------------------------------------------------------------------------
// Description :
//		항공기 이동에 따른 지형DB 에 업데이트 유무를 확인
//
// Arguments : 
//		double lat  - 위도
//		double lon  - 경도
//
// Retrun :
//		0 - 업데이트 필요
//		1 - 현재 캐시되어 있음.
//------------------------------------------------------------------------------------------------------
int IsUpdateBoundary(int type, double lat, double lon)

{	if ( lat >=  _MapData[type].mapCache.updateBoundary.leftBottom.lat
		&& lon >= _MapData[type].mapCache.updateBoundary.leftBottom.lon
		&& lat <= _MapData[type].mapCache.updateBoundary.rightTop.lat
		&& lon <= _MapData[type].mapCache.updateBoundary.rightTop.lon)
	
{		return 1;
	}

	return 0;
}



//------------------------------------------------------------------------------------------------------
// Description :
//		고도정보 조회
//------------------------------------------------------------------------------------------------------
double MAPGetTerrainHeight( double lat, double lon)

{	int leftLonIndex = IndexOfLon(0, lon);
	int bottomLatIndex = IndexOfLat(0, lat);

	double h_DB = -1.0;
	int map_off = 0;

	if ((bottomLatIndex - 1) >= 0
		&& bottomLatIndex < MAPCACHE_SIZE
		&& leftLonIndex >= 0
		&& (leftLonIndex + 1) < MAPCACHE_SIZE)
	
{		short P11 = _MapData[0].mapCache.map[leftLonIndex][bottomLatIndex];			// Bottom Left 		lat lower
		short P12 = _MapData[0].mapCache.map[leftLonIndex][bottomLatIndex + 1];		// top left 
		short P21 = _MapData[0].mapCache.map[leftLonIndex + 1][bottomLatIndex];		// bottom right				
		short P22 = _MapData[0].mapCache.map[leftLonIndex + 1][bottomLatIndex + 1];	// top right 

		if (P11 < 0 || P11 > 9000)
		
{			map_off = 1;
		}

		if (P12 < 0 || P12 > 9000)
		
{			map_off = 1;
		}

		if (P21 < 0 || P21 > 9000)
		
{			map_off = 1;
		}

		if (P22 < 0 || P22 > 9000)
		
{			map_off = 1;
		}

		if (map_off == 0)
{			double leftLon = _MapData[0].mapCache.scanBoundary.leftBottom.lon + (leftLonIndex * _MapData[0].mapCache.resolution * ARCSEC2DEG);
			double bottomLat = _MapData[0].mapCache.scanBoundary.leftBottom.lat + (bottomLatIndex  * _MapData[0].mapCache.resolution * ARCSEC2DEG);

			double x = (lon - leftLon)*DEG2RAD;
			double y = (lat - bottomLat)*DEG2RAD;
			double d = _MapData[0].mapCache.resolution * ARCSEC2DEG * DEG2RAD;

			long double P_low = (d - x) / d * P11 + (x / d) * P21;
			long double P_high = (d - x) / d * P12 + (x / d) * P22;
			h_DB = (d - y) / d * P_low + (y / d) * P_high;
		}
	}
	else
{		h_DB = - 1.0;
	}

	return h_DB;
}


int IndexOfLat(int type, double lat)

{	double d = _MapData[type].mapCache.resolution * ARCSEC2DEG * DEG2RAD;

	int index = (int)((((lat - _MapData[type].mapCache.scanBoundary.leftBottom.lat) * DEG2RAD) / d));
	if (index < 0 || index >= MAPCACHE_SIZE)
	
{		return -1;
	}

	return index;
}

int IndexOfLon(int type, double lon)

{	double d = _MapData[type].mapCache.resolution * ARCSEC2DEG * DEG2RAD;

	int index = (int)((((lon - _MapData[type].mapCache.scanBoundary.leftBottom.lon)* DEG2RAD) / d));
	if (index < 0 || index >= MAPCACHE_SIZE)
	
{		return -1;
	}

	return index;
}

int Resolution(int type)

{	int resolution = 3;
	if (_MapData[type].mapCache.resolution != 0)
	
{		resolution = _MapData[type].mapCache.resolution;
	}
	return resolution;
}

//------------------------------------------------------------------------------------------------------
// Description :
//		Post 방식 지형/지물 조회를 위한 기준위치 지정
//
// Arguments : 
//		int type - 지형/지물 유형 (0 : 지형, 1: 지물)
//		double lat  - 위도
//		double lon  - 경도
//
// Retrun :
//		핸들
//------------------------------------------------------------------------------------------------------
int MAPGetPostHandle(int type, double lat, double lon)

{	int handle = -1;
	int leftLonIndex = IndexOfLon(type, lon);
	int bottomLatIndex = IndexOfLat(type, lat);
	if ((bottomLatIndex - 1) >= 0
		&& bottomLatIndex < MAPCACHE_SIZE
		&& leftLonIndex >= 0
		&& (leftLonIndex + 1) < MAPCACHE_SIZE)

{		// handle
		handle = _map_postHandleFront;

		// set parameters
		_map_postHandle[handle].type = type;
		_map_postHandle[handle].x_offset = leftLonIndex;
		_map_postHandle[handle].y_offset = bottomLatIndex;

		// moved a front point of the queue
		_map_postHandleFront = (_map_postHandleFront + 1) % 5;
	}

	return handle;	
}


//------------------------------------------------------------------------------------------------------
// Description :
//		Post 방식 지형/지물 조회
//
// Arguments : 
//		int handle - MAPGetPostHandle()함수로 부터 얻은 기준위치에 대한 핸들
//		int x_offset  - 경도 방향 오프셋
//		int y_offset  - 위도 방향 오프셋
//
// Retrun :
//		지형/지물 고도
//------------------------------------------------------------------------------------------------------
short MAPGetPostData(int handle, int x_offset, int y_offset)

{	double h_DB = -1.0;
	int map_off = 0;
	short height = -1.0;

	if (handle >= 0 && handle < 5)

{		int type = _map_postHandle[handle].type;
		int x = _map_postHandle[handle].x_offset + x_offset;
		int y = _map_postHandle[handle].y_offset + y_offset;

		if (type ==0 || type == 1)
{			if (x >= 0
				&& x < MAPCACHE_SIZE
				&& y >= 0
				&& y < MAPCACHE_SIZE)

{				height = _MapData[type].mapCache.map[x][y];
			}
		}
	}
	
	return height;
}