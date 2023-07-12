#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MAP.h"


// MAP���� ��ġ
// MAP����(����34 �浵126 ~ ����38 �浵130)
#define FOLDER_PATH "C:\\Users\\JHC\\Documents\\DTNS\\MAP\\MAP_DB"

// DVOF ��ȸ��� ��Ȱ��ȭ (��� 0���� ���)
// Only DTED Lv2 ������ȸ��
int main(){

	// �ʱ� ��ġ
	double lat = 37.2;
	double lon = 128.1;
		
	// �� ���� ��ġ �ʱ�ȭ
	int result = InitialMap(FOLDER_PATH, strlen(FOLDER_PATH));
	if (result != 0){
		printf("MAP Folder Path is wrong!\n");
		return -1;
	}

	// MAP ��ȸ �׽�Ʈ
	for (int row = 0; row < 1000; row++)
	{
		// ��(����/����) ĳ�� ������Ʈ (2160 x 2160)
		UpdateCache(lat, lon);

		// Interpolation ��� ���� ��ȸ
		// Interpolation ����� ���� ��ȸ���� ����.
		double height = MAPGetTerrainHeight(lat, lon);
		printf("Terrain Height(%lf,%lf) : %lfm\n", lat, lon, height);

		// Post��� �� ��ȸ (����)
		int handle1 = MAPGetPostHandle(0, lat, lon);
		short height1 = MAPGetPostData(handle1, 0, 0);
		short height2 = MAPGetPostData(handle1, 1, 1);
		short height3 = MAPGetPostData(handle1, -1, -1);
		short height4 = MAPGetPostData(handle1, 600, -600);

		// Post��� �� ��ȸ (����) 
		int handle2 = MAPGetPostHandle(1, lat, lon);
		short height5 = MAPGetPostData(handle2, 0, 0);
		short height6 = MAPGetPostData(handle2, 600, -600);
		
		lat -= 0.001;
		lon -= 0.001;
	}

	return 0;
}