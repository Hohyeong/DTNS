#ifndef _MAPHANDLER_HEADER
#define _MAPHANDLER_HEADER

#define CELLCACHE_SIZE	(720)
#define CELLCNT			(3)
#define MAPCACHE_SIZE	(CELLCACHE_SIZE*CELLCNT)

#define DEVIDE_CNT		(5)
#define CELL_LAT_CNT	(180*5)
#define CELL_LON_CNT	(360*5)
#define CELL_X_CENTER	(CELL_LON_CNT/2-1)
#define CELL_Y_CENTER	(CELL_LAT_CNT/2-1)
#define CELL_INTERVAL	(0.2)
#define RESOLUTION		(1)
#define ARCSEC2DEG      (0.000277777777777777)
#define DEG2RAD         (0.0174532925199433)
#define RAD2DEG         (57.2957795131)

#define CELLUNIT		(60/DEVIDE_CNT)

//---------------------------------------------------------------------------------------------
// MAP DATA AREA
//---------------------------------------------------------------------------------------------
#define BUF_SIZE		(1024)
#define SUCCESS			(0)
#define FAIL			(-1)

#define NO_LOADED		(0x01)
#define NOT_EXIST_MAP	(0x02)
#define MAPUPDATE_ERROR (0x04)
#define DVOFUPDATE_ERRO (0x08)


#ifdef	__cplusplus
extern "C" {
#endif

#pragma pack(push, 1)

typedef struct
{
	double lat;
	double lon;
} COORDINATE;

typedef struct
{
	COORDINATE leftBottom;
	COORDINATE leftTop;
	COORDINATE rightTop;
	COORDINATE rightBottom;
} BOUNDARY;

typedef struct
{
	int resolution;
	BOUNDARY updateBoundary;
	BOUNDARY scanBoundary;

	//------------------------------------------------------------------------------------------------------
	//                         most left longitude          
	//                         ↓
	//                        -------------
	// most bottom latiude → | ? | ? | ? |  [0] 
	//                        | ? | ? | ? |  [1]
	//                        | ? | ? | ? |  [2]
	//                        -------------
	//                         [0] [1] [2]
	//------------------------------------------------------------------------------------------------------
	short map[MAPCACHE_SIZE][MAPCACHE_SIZE];

} MAP_CACHE;

typedef union
{
	unsigned short raw[7];
	struct
	{
		/* 0-2bit : resolution, 4bit : east, 5 bit : north */
		unsigned short flag;
		unsigned short lonDeg;
		unsigned short lonMinute;
		unsigned short latDeg;
		unsigned short latMinute;
		unsigned short numOfLat;
		unsigned short numOfLon;
	} data;
} MAP_HEADER;

typedef struct
{
	MAP_CACHE	mapCache;
	//int			mapError;
	char		mapPath[BUF_SIZE];
} MAPDATA;

typedef struct{
	int type;	// 0: Terrain , 1: Obstacle
	int x_offset;
	int y_offset;
} MAP_POST_HANDLE;

#pragma pack(pop) 


int UpdateCache(double lat, double lon);
int UpdateMapCache(double lat, double lon);
int UpdateObsCache(double lat, double lon);
void SetUpdateBoundary(int type, double lat, double lon);
void SetScanBoundary();
int IsUpdateBoundary(int type, double lat, double lon);
int IndexOfLat(int type, double lat);
int IndexOfLon(int type, double lon);
int Resolution(int type);
int ReadMaps(double lat, double lon);
int ReadCell(int x, int y, int latDeg, int latMin, int lonDeg, int lonMin);
void ProcessMapOff(int type, int x, int y);

//-------------------------------------------------------------
// Simulink기반 시험에 필요한 함수
int InitialMap(char* path, int size);
int UpdateCache(double lat, double lon);
//-------------------------------------------------------------
// 알고리즘 내에서 사용 가능 함수
double MAPGetTerrainHeight(double lat, double lon);
int MAPGetPostHandle(int type, double lat, double lon);
short MAPGetPostData(int handle, int x_offset, int y_offset);
//-------------------------------------------------------------

extern MAPDATA		_MapData[2];

#ifdef	__cplusplus
}
#endif

#endif