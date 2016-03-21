#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MainDirection.h"
#include "timeit.h"

void loadLog(char* path, DataVec* dv)
{
	int (*p)[3];
	char line[128];
	FILE* file = fopen(path, "r");
	for (p = dv->data, dv->len = 0; p < dv->data + DATA_MAXLEN; ++p, ++(dv->len)) {
		if (!fgets(line, sizeof(line), file)) break;
		if (sscanf(line, "%*f %d %d %d", &(*p)[0], &(*p)[1], &(*p)[2]) != 3) break;
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

#define st2float(x) ((float)(x)/(float)(1<<STVAL_SHIFT))
void show(char* path)
{
	DataVec dv;
	StVals st;
	real realCov[3][3];
	real eigenVal[3];
	real eigenVec[3][3];

	// load
	loadLog(path, &dv);

	// statistic
	calcStatistics(&dv, &st);
	printf("Avg:" "\t{ %.1f, %.1f, %.1f }\n",
				  st2float(st.avg[0]), st2float(st.avg[1]), st2float(st.avg[2]));
	printf("Cov:" "\t{{ %.1f, %.1f, %.1f },\n"
		          "\t { %.1f, %.1f, %.1f },\n"
				  "\t { %.1f, %.1f, %.1f }}\n", 
				  st2float(st.cov[0][0]), st2float(st.cov[0][1]), st2float(st.cov[0][2]),
				  st2float(st.cov[1][0]), st2float(st.cov[1][1]), st2float(st.cov[1][2]),
				  st2float(st.cov[2][0]), st2float(st.cov[2][1]), st2float(st.cov[2][2]));

	stcov2real(&st, realCov);
	eigenSystem(realCov, eigenVal, eigenVec);
	printf("lambda:" "\t{ %.3f, %.3f, %.3f }\n",
		real2float(eigenVal[0]), real2float(eigenVal[1]), real2float(eigenVal[2]));
	printf("vector:" "\t{{ %.3f, %.3f, %.3f },\n"
					 "\t { %.3f, %.3f, %.3f },\n"
					 "\t { %.3f, %.3f, %.3f }}\n", 
					 real2float(eigenVec[0][0]), real2float(eigenVec[0][1]), real2float(eigenVec[0][2]),
					 real2float(eigenVec[1][0]), real2float(eigenVec[1][1]), real2float(eigenVec[1][2]),
					 real2float(eigenVec[2][0]), real2float(eigenVec[2][1]), real2float(eigenVec[2][2]));

	char plt[1024];
	sprintf(plt, "%s.plt", path);
	FILE* file = fopen(plt, "w");
	fprintf(file,
		    "set view equal xyz\n"
		    "set arrow from %.3f,%.3f,%.3f rto %.3f,%.3f,%.3f lc rgb 'red'\n"
			"set arrow from %.3f,%.3f,%.3f rto %.3f,%.3f,%.3f lc rgb 'green'\n"
			"set arrow from %.3f,%.3f,%.3f rto %.3f,%.3f,%.3f lc rgb 'blue'\n"
			"splot '%s' u 2:3:4 w d lc rgb 'black'",
			st2float(st.avg[0]), st2float(st.avg[1]), st2float(st.avg[2]),
			real2float(eigenVec[0][0]) * 128, real2float(eigenVec[0][1]) * 128, real2float(eigenVec[0][2]) * 128,
			st2float(st.avg[0]), st2float(st.avg[1]), st2float(st.avg[2]),
			real2float(eigenVec[1][0]) * 128, real2float(eigenVec[1][1]) * 128, real2float(eigenVec[1][2]) * 128,
			st2float(st.avg[0]), st2float(st.avg[1]), st2float(st.avg[2]),
			real2float(eigenVec[2][0]) * 128, real2float(eigenVec[2][1]) * 128, real2float(eigenVec[2][2]) * 128,
			filename(path));
	fclose(file);
}

void deal(char* prereal)
{
	char buff[1024];
	sprintf(buff, "dir /B \"%s*.log\" > files.txt", prereal);
	system(buff);
	FILE* file = fopen("files.txt", "r");
	int len = sprintf(buff, "%s", prereal);
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
#if 1
	char* path = argc > 1 ? argv[1] : "Data\\501053.log";
	TIMEIT(1, show(path));
#else
	deal("Data\\");
#endif
	system("pause");
}
#endif
