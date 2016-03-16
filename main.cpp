#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MainDirection.h"

#define TIME_TEST(times, call) do {		\
		int i;							\
		clock_t t0 = clock();			\
		for (i = 0; i < times; i++ ) {	\
			call;						\
		}								\
		clock_t t1 = clock();			\
		printf("average running time: %f\n", (float)(t1 - t0) / times / CLOCKS_PER_SEC);	\
	} while(0)

#define DATA_RATIO 1

void loadLog(char* path, DataVec* dv)
{
	const int K = 1;
	int x, y, z;
	int (*p)[3];
	char line[128];
	FILE* file = fopen(path, "r");
	for (p = dv->data, dv->len = 0; p < dv->data + DATA_MAXLEN; ++p, ++(dv->len)) {
		if (!fgets(line, sizeof(line), file)) break;
		if (sscanf(line, "%*f %d %d %d", &x, &y, &z) != 3) break;
		(*p)[0] = x * DATA_RATIO;
		(*p)[1] = y * DATA_RATIO;
		(*p)[2] = z * DATA_RATIO;
	}
	fclose(file);
}

char* filename(char* path)
{
	char* p = NULL;
	if (!p) p = strrchr(path, '\\');
	if (!p) p = strrchr(path, '/');
	if (!p) return path;
	return p + 1;
}

void show(char* path)
{
	DataVec dv;
	StVals st;
	float eigenVal[3];
	float eigenVec[3][3];

	// load
	loadLog(path, &dv);

	// statistic
	calcStatistics(&dv, &st);
	printf("Avg:" "\t{ %ld, %ld, %ld }\n", st.avg[0], st.avg[1], st.avg[2]);
	printf("Cov:" "\t{{ %ld, %ld, %ld },\n"
		          "\t { %ld, %ld, %ld },\n"
				  "\t { %ld, %ld, %ld }}\n", 
				  st.cov[0][0], st.cov[0][1], st.cov[0][2],
				  st.cov[1][0], st.cov[1][1], st.cov[1][2],
				  st.cov[2][0], st.cov[2][1], st.cov[2][2]);

	// 
	eigenSystem(st.cov, eigenVal, eigenVec);
	printf("lambda:" "\t{ %f, %f, %f }\n", eigenVal[0], eigenVal[1], eigenVal[2]);
	printf("vector:" "\t{{ %f, %f, %f },\n"
					 "\t { %f, %f, %f },\n"
					 "\t { %f, %f, %f }}\n", 
					 eigenVec[0][0], eigenVec[0][1], eigenVec[0][2],
					 eigenVec[1][0], eigenVec[1][1], eigenVec[1][2],
					 eigenVec[2][0], eigenVec[2][1], eigenVec[2][2]);

	char plt[1024];
	sprintf(plt, "%s.plt", path);
	FILE* file = fopen(plt, "w");
	fprintf(file,
		    "set view equal xyz\n"
		    "set arrow from %ld,%ld,%ld rto %f,%f,%f lc rgb 'red'\n"
			"set arrow from %ld,%ld,%ld rto %f,%f,%f lc rgb 'green'\n"
			"set arrow from %ld,%ld,%ld rto %f,%f,%f lc rgb 'blue'\n"
			"splot '%s' u 2:3:4 w d lc rgb 'black'",
			st.avg[0], st.avg[1], st.avg[2],
			eigenVec[0][0] * 128, eigenVec[0][1] * 128, eigenVec[0][2] * 128,
			st.avg[0], st.avg[1], st.avg[2],
			eigenVec[1][0] * 128, eigenVec[1][1] * 128, eigenVec[1][2] * 128,
			st.avg[0], st.avg[1], st.avg[2],
			eigenVec[2][0] * 128, eigenVec[2][1] * 128, eigenVec[2][2] * 128,
			filename(path));
	fclose(file);
}

void deal(char* prefix)
{
	char buff[1024];
	sprintf(buff, "dir /B \"%s*.log\" > files.txt", prefix);
	system(buff);
	FILE* file = fopen("files.txt", "r");
	int len = sprintf(buff, "%s", prefix);
	while (fgets(buff+len, sizeof(buff) - len, file)) {
		buff[strlen(buff) - 1] = '\0';
		printf("\nDeal with %s\n", buff);
		show(buff);
	}
	fclose(file);
}

#ifndef TEST
int main(int argc, char* argv[])
{
	//char* path = argc > 1 ? argv[1] : "501053.log";
	//TIME_TEST(1, show(path));
	deal("Data\\");
	system("pause");
}
#endif
