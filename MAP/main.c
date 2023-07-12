#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MAP.h"


// MAP폴더 위치
// MAP영역(위도34 경도126 ~ 위도38 경도130)
#define FOLDER_PATH "C:\\Users\\JHC\\Documents\\DTNS\\MAP\\MAP_DB"

// DVOF 조회기능 비활성화 (모두 0으로 출력)
// Only DTED Lv2 지형조회용
int main(){

	// 초기 위치
	double lat = 37.2;
	double lon = 128.1;
		
	// 맵 폴더 위치 초기화
	int result = InitialMap(FOLDER_PATH, strlen(FOLDER_PATH));
	if (result != 0){
		printf("MAP Folder Path is wrong!\n");
		return -1;
	}

	// MAP 조회 테스트
	for (int row = 0; row < 1000; row++)
	{
		// 맵(지형/지물) 캐시 업데이트 (2160 x 2160)
		UpdateCache(lat, lon);

		// Interpolation 방식 지형 조회
		// Interpolation 방식은 지형 조회에만 있음.
		double height = MAPGetTerrainHeight(lat, lon);
		printf("Terrain Height(%lf,%lf) : %lfm\n", lat, lon, height);

		// Post방식 맵 조회 (지형)
		int handle1 = MAPGetPostHandle(0, lat, lon);
		short height1 = MAPGetPostData(handle1, 0, 0);
		short height2 = MAPGetPostData(handle1, 1, 1);
		short height3 = MAPGetPostData(handle1, -1, -1);
		short height4 = MAPGetPostData(handle1, 600, -600);

		// Post방식 맵 조회 (지물) 
		int handle2 = MAPGetPostHandle(1, lat, lon);
		short height5 = MAPGetPostData(handle2, 0, 0);
		short height6 = MAPGetPostData(handle2, 600, -600);
		
		lat -= 0.001;
		lon -= 0.001;
	}

	return 0;
}