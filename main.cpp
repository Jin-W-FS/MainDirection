#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MainDirection.h"
#include "timeit.h"

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
	fixed eigenVal[3];
	fixed eigenVec[3][3];

	// load
	loadLog(path, &dv);

	// statistic
	calcStatistics(&dv, &st);
	printf("Avg:" "\t{ %.3f, %.3f, %.3f }\n", 
		fixed2float(st.avg[0]), fixed2float(st.avg[1]), fixed2float(st.avg[2]));
	printf("Cov:" "\t{{ %.3f, %.3f, %.3f },\n"
		          "\t { %.3f, %.3f, %.3f },\n"
				  "\t { %.3f, %.3f, %.3f }}\n", 
				  fixed2float(st.cov[0][0]), fixed2float(st.cov[0][1]), fixed2float(st.cov[0][2]),
				  fixed2float(st.cov[1][0]), fixed2float(st.cov[1][1]), fixed2float(st.cov[1][2]),
				  fixed2float(st.cov[2][0]), fixed2float(st.cov[2][1]), fixed2float(st.cov[2][2]));

	 
	eigenSystem(st.cov, eigenVal, eigenVec);
	printf("lambda:" "\t{ %.3f, %.3f, %.3f }\n",
		fixed2float(eigenVal[0]), fixed2float(eigenVal[1]), fixed2float(eigenVal[2]));
	printf("vector:" "\t{{ %.3f, %.3f, %.3f },\n"
					 "\t { %.3f, %.3f, %.3f },\n"
					 "\t { %.3f, %.3f, %.3f }}\n", 
					 fixed2float(eigenVec[0][0]), fixed2float(eigenVec[0][1]), fixed2float(eigenVec[0][2]),
					 fixed2float(eigenVec[1][0]), fixed2float(eigenVec[1][1]), fixed2float(eigenVec[1][2]),
					 fixed2float(eigenVec[2][0]), fixed2float(eigenVec[2][1]), fixed2float(eigenVec[2][2]));

	char plt[1024];
	sprintf(plt, "%s.plt", path);
	FILE* file = fopen(plt, "w");
	fprintf(file,
		    "set view equal xyz\n"
		    "set arrow from %.3f,%.3f,%.3f rto %.3f,%.3f,%.3f lc rgb 'red'\n"
			"set arrow from %.3f,%.3f,%.3f rto %.3f,%.3f,%.3f lc rgb 'green'\n"
			"set arrow from %.3f,%.3f,%.3f rto %.3f,%.3f,%.3f lc rgb 'blue'\n"
			"splot '%s' u 2:3:4 w d lc rgb 'black'",
			fixed2float(st.avg[0]), fixed2float(st.avg[1]), fixed2float(st.avg[2]),
			fixed2float(eigenVec[0][0]) * 128, fixed2float(eigenVec[0][1]) * 128, fixed2float(eigenVec[0][2]) * 128,
			fixed2float(st.avg[0]), fixed2float(st.avg[1]), fixed2float(st.avg[2]),
			fixed2float(eigenVec[1][0]) * 128, fixed2float(eigenVec[1][1]) * 128, fixed2float(eigenVec[1][2]) * 128,
			fixed2float(st.avg[0]), fixed2float(st.avg[1]), fixed2float(st.avg[2]),
			fixed2float(eigenVec[2][0]) * 128, fixed2float(eigenVec[2][1]) * 128, fixed2float(eigenVec[2][2]) * 128,
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
	//char* path = argc > 1 ? argv[1] : "Data\\501053.log";
	//TIMEIT(1, show(path));
	deal("Data\\");
	system("pause");
}
#endif
