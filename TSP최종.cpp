/* 2015-11-25 ~ 2015-11-26 김영해
인공지능 유전자 알고리즘 - 상호 교환 돌연변이 사용
순회방문원문제
도시 1에서 시작하여 모든 도시를 방문 후
도시 1로 귀환하는 가장 짧은 경로 구하기 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <share.h>

#define DATANO		30		// 데이터 수
#define POOLSIZE	30		// 풀 사이즈
#define LASTG		600		// 중단 세대
#define MRATE		0.01	// 돌연변이 확률
#define STARTCITY	1		// 시작 및 돌아올 도시 번호
#define YES			1
#define NO			0
#define FILEPATH	"..\\result_GA.txt"	// 파일 경로, 이 프로젝트 폴더에 생성한다.

//char pathbuf[100];			// Bestfit의 경우 최종 경로를 출력하기 위한 buffer
FILE *fp;						// 결과를 파일로 저장하기 위해 사용할 파일 포인터

double calcDistance(double x1, double x2, double y1, double y2);
double evalfit(int gene[]);					// 적응도 계산
void mating(int pool[POOLSIZE][DATANO]);	// 교차
void mutation(int pool[POOLSIZE][DATANO]);	// 돌연변이
void printp(int pool[POOLSIZE][DATANO], int generation);	// 결과 출력
void initpool(int pool[POOLSIZE][DATANO]);	// 초기 집단 생성 

double rndn(double l);								// n 미만의 난수 생성 double
int rndn2(int l);									// n 미만의 난수 생성 int
int rndnInit(int cityPool[], int canSelectSize);	// 초기 집단 생성시 불규칙하게 도시 선택
int select(double roulette[POOLSIZE], double totalfitness);	// 부모 선택
void copypool(int pool[POOLSIZE][DATANO], int nextpool[POOLSIZE][DATANO]);	// 다음 세대의 풀 복사
void crossing(int m[], int p[], int c1[], int c2[]);	// 유전자 교차

double q[][2] = {
	3.9, 2.1, 9.7, 7.1, 8.2, 6.8, 2.9, 0.2, 7.2, 0.7,
	6.5, 0.3, 3.1, 7.4, 2.7, 2.2, 1.1, 3.4, 0.6, 8.8,
	7.4, 8.6, 9.4, 5.4, 2.7, 7.6, 4.1, 8.5, 6.7, 1.1,
	2.0, 3.9, 1.7, 4.4, 0.1, 9.9, 9.9, 3.9, 4.6, 1.6,
	8.4, 0.9, 0.7, 9.1, 2.2, 3.4, 6.5, 2.7, 5.6, 5.4,
	7.8, 3.3, 3.6, 6.3, 8.9, 7.1, 3.8, 2.9, 4.7, 5.2
};

double poolAver[LASTG];				// 세대들의 평균 적합도 저장 전역배열
double poolElite[LASTG];				// 세대들의 엘리트 저장 전역 배열
double poolSuper = 0;						// 세대들의 슈퍼엘리트 저장 전역 배열

int main() {
	int pool[POOLSIZE][DATANO];				// 유전자 풀
	int generation;							// 현재 세대수
	clock_t start = clock();				// 실행 시작 시간

	srand((unsigned int)time(NULL));		// 난수 초기화

	if ((fp = _fsopen(FILEPATH, "w+", _SH_DENYNO)) == NULL)
		perror("File open error");

	initpool(pool);							// 초기 집단 생성
	fprintf(fp, "pool[초기집단] : \n");
	for (int i = 0; i < POOLSIZE; i++) {	// 초기집단 출력
		for (int j = 0; j < DATANO; j++) {
			fprintf(fp, "%d ", pool[i][j]);
		}
		fprintf(fp, "\n\n");
	}

	for (generation = 0; generation < LASTG; ++generation) {	// 중단 세대까지 반복
		fprintf(fp, "%d세대\n", generation);
		mating(pool);		// 교차
		mutation(pool);		// 돌연변이
		printp(pool, generation);		// 결과 출력
	}

	fprintf(fp, "세대\t\t평균\t\t엘리트\n");
	for (generation = 0; generation < LASTG; ++generation) {	// 정리하여 결과 출력
		fprintf(fp, "%d세대\t %lf \t", generation, poolAver[generation]);
		fprintf(fp, "%lf\n", poolElite[generation]);
	}
	fprintf(fp, "SuperElite : %lf\n", poolSuper);



	fprintf(fp, "Execution time: %.3lfs\n", (double)(clock() - start) / CLOCKS_PER_SEC);

	fclose(fp);

	return 0;
}

int select(double roulette[POOLSIZE], double totalfitness) {
	int i;				// 반복 제어 문자
	double ball;		// 공
	double acc = 0;		// 적응도 적산값

	for (i = 0; i < POOLSIZE; ++i) {
		acc += totalfitness / roulette[i];	// 거리가 가까울수록 적합하기 때문에 역수를 취함
	}

	ball = rndn(acc);	// 역수의 합들
	acc = 0;
	for (i = 0; i < POOLSIZE; ++i) {
		acc += totalfitness / roulette[i];
		if (acc > ball)	// 대응하는 유전자
			break;
	}

	return i;
}

void mating(int pool[POOLSIZE][DATANO]) {
	int i;							// 반복 제어 함수
	double totalfitness = 0;		// 적응도 합계
	int nextpool[POOLSIZE][DATANO];	// 자식 세대
	double roulette[POOLSIZE];		// 적응도 저장
	int mama, papa;					// 부모 유전자 번호

	for (i = 0; i < POOLSIZE; ++i) {	// 룰렛 작성
		roulette[i] = evalfit(pool[i]);	// 적응도 계산

		totalfitness += roulette[i];	// 적응도 합계 계산
	}

	for (i = 0; i < POOLSIZE / 2; ++i) {// 선택과 교차를 반복
		do {	// 부모 선택
			mama = select(roulette, totalfitness);
			papa = select(roulette, totalfitness);
		} while (mama == papa);// 부모 중복 제거
		crossing(pool[mama], pool[papa], nextpool[i * 2], nextpool[i * 2 + 1]);	// 두 유전자 교차
	}
	copypool(pool, nextpool);	// 다음 세대 풀 복사

}

void crossing(int m[], int p[], int c1[], int c2[]) {
	int i, j;		// 반복 제어 변수
	int cp1, cp2;	// 교차할 점 2개
	int tmp;		// 교환 임시 변수
	c1[DATANO - 1] = STARTCITY;	// 마지막에는 무조건 STARTCITY
	c2[DATANO - 1] = STARTCITY;

	cp1 = rndn2(DATANO - 1);	// 교차점 결정 마지막 도시는 무조건 STARTCITY여야 하므로 마지막 제외 0 ~ 18
	cp2 = rndn2(DATANO - 1);

	if (cp1 > cp2) {			// cp2가 더 크도록 순서 교정
		tmp = cp1;
		cp1 = cp2;
		cp2 = tmp;
	}

	for (i = cp1; i <= cp2; i++) {	// 선택된 부분 교차하여 자식에게 줌
		c1[i] = p[i];
		c2[i] = m[i];
	}

	int isHaving;	// 가지고 있는 인자인가?
	int INDEX;		// 가지고 있지 않은 인자를 자식에게 넣을 인덱스

	INDEX = 0;		// 자식의 첫 번째에 넣는다
	for (i = 0; i < DATANO - 1; i++) {	// c1에게 mama가 없는 인자를 채워줌
		for (j = cp1; j <= cp2; j++) {
			if (c1[j] == m[i]) {
				isHaving = YES;			// 자식이 가지고 있다면 필요 없다.
				break;
			}
			else if (c1[j] != m[i])		// 확인해 본 곳까지의 자식이 인자를 가지고 있지 않으면
				isHaving = NO;
		}

		if (isHaving == NO) {
			if (INDEX < cp1) {			// cp1 앞쪽에 가지고 있지 않은 인자를 넣음
				c1[INDEX] = m[i];
				INDEX++;
				if (INDEX == cp1)
					INDEX = cp2 + 1;
			}
			else if (INDEX == cp1) {	// cp1 앞쪽이 다찼다면 cp2뒤로 옮김
				INDEX = cp2 + 1;
				c1[INDEX] = m[i];
				INDEX++;
			}
			else if (INDEX > cp2) {		// cp2 뒤쪽에 가지고 있지 않는 인자를 넣음
				c1[INDEX] = m[i];
				INDEX++;
			}
		}

	}

	INDEX = 0;		// 자식의 첫번째에 넣는다
	for (i = 0; i < DATANO - 1; i++) {	// c2에게 papa가 없는 인자를 채워줌
		for (j = cp1; j <= cp2; j++) {
			if (c2[j] == p[i]) {
				isHaving = YES;			// 자식이 가지고 있다면 필요 없다.
				break;
			}
			else if (c2[j] != p[i])		// 확인해 본 곳까지의 자식이 인자를 가지고 있지 않으면
				isHaving = NO;
		}

		if (isHaving == NO) {
			if (INDEX < cp1) {			// cp1 앞쪽에 가지고 있지 않은 인자를 넣음
				c2[INDEX] = p[i];
				INDEX++;
				if (INDEX == cp1)
					INDEX = cp2 + 1;
			}
			else if (INDEX == cp1) {	// cp1 앞쪽이 다찼다면 cp2뒤로 옮김
				INDEX = cp2 + 1;
				c2[INDEX] = p[i];
				INDEX++;
			}
			else if (INDEX > cp2) {		// cp2 뒤쪽에 가지고 있지 않는 인자를 넣음
				c2[INDEX] = p[i];
				INDEX++;
			}
		}
	}
}


void copypool(int pool[POOLSIZE][DATANO], int nextpool[POOLSIZE][DATANO]) {
	int i, j;//반복 제어 함수

	for (i = 0; i < POOLSIZE; i++) {
		for (j = 0; j < DATANO; j++) {
			pool[i][j] = nextpool[i][j];
		}
	}

}

double calcDistance(double x1, double x2, double y1, double y2) {
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

double evalfit(int g[]) {
	int i;				// 반복 제어 변수
	double fitness = 0;	// 적응도 도중값

	for (i = 0; i < DATANO; i++) {
		if (i == 0)
			// 처음 방문 도시와 STARTCITY과 거리
			fitness += calcDistance(q[g[i] - 1][0], q[STARTCITY - 1][0], q[g[i] - 1][1], q[STARTCITY - 1][1]);
		else
			// 다음 방문 도시와 현재 도시의 거리
			fitness += calcDistance(q[g[i] - 1][0], q[g[i - 1] - 1][0], q[g[i] - 1][1], q[g[i - 1] - 1][1]);
	}

	return fitness;
}

void printp(int pool[POOLSIZE][DATANO], int generation) {
	int i, j;				// 반복 제어 변수
	double fitness;			// 적응도
	double totalfitness = 0;// 적응도 합계
	int elite;
	double bestfit = 0;		// 엘리트 유전자 처리용 변수

	for (i = 0; i < POOLSIZE; i++) {
		for (j = 0; j < DATANO; j++) {
			fprintf(fp, "%1d ", pool[i][j]);
		}
		fitness = evalfit(pool[i]);
		fprintf(fp, "\t%lf\n", fitness);

		if (bestfit == 0) {	//bestfit 초기 값 설정
			bestfit = fitness;
			elite = i;
		}
		else if (fitness < bestfit) {	// 엘리트 해
			bestfit = fitness;
			elite = i;
		}

		totalfitness += fitness;
	}

	fprintf(fp, "\nElite: %d\tBestfit: %lf\n", elite, bestfit);			// 엘리트 해의 적응도 출력
	fprintf(fp, "Average Fitness: %lf\n\n", totalfitness / POOLSIZE);	// 평균 적응도 출력

																		/*정리된 결과 전역변수들에 저장*/
	poolAver[generation] = totalfitness / POOLSIZE;		//해당 세대의 평균 적응도 저장
	poolElite[generation] = bestfit;					//해당 세대의 엘리트 적응도 저장
	if (poolSuper == 0)		//super 엘리트 초기값 설정
		poolSuper = bestfit;
	else if (poolSuper > bestfit)
		poolSuper = bestfit;							//전체 세대에서 찾은 슈퍼엘리트 저장
}

void initpool(int pool[POOLSIZE][DATANO]) {
	int i, j, k;			// 반복 제어변수
	int cityPool[DATANO];	// 선택 가능한 도시의 집합
	int canSelecSize;		// 선택 가능한 도시들의 갯수

	for (i = 0; i < POOLSIZE; ++i) {
		canSelecSize = DATANO - 1;		// STARTCITY를 선택하지 않기 위해 - 1

		for (k = 0; k < DATANO; k++) {	// 선택 가능 도시의 집합 초기화
			if (k == 0)
				cityPool[k] = STARTCITY;// 선택되지 않는 배열의 0번 인자에 시작 도시를 넣음
			cityPool[k] = k + 1;
		}

		for (j = 0; j < DATANO; ++j) {
			pool[i][j] = rndnInit(cityPool, canSelecSize);	// 마지막을 제외한 도시 순서를 랜덤하게 설정
			canSelecSize--;									// 선택된 도시 제 선택 되지 않도록 범위에서 제외
		}
	}
}

double rndn(double l) {
	double rndno;	// 생성한 난수

	while ((rndno = ((double)rand() / RAND_MAX) * l) == l);// 인자값 미만의 난수생성 

	return rndno;
}

int rndn2(int l) {
	int rndno;		// 생성한 난수

	while ((rndno = ((double)rand() / RAND_MAX) * l) == l);// 인자값 미만의 난수생성 

	return rndno;
}

int rndnInit(int cityPool[], int canSelectSize) {
	int rndCity;		// 생성한 난수
	int tmp;			// 배열의 위치를 바꾸기 위한 임시 변수
	int selectedCity;	// 선택된 도시, 

	if (canSelectSize == 0) {						// 만약 도시가 한개 남았다면
		selectedCity = STARTCITY - 1;
		return rndCity = cityPool[selectedCity];	// STARTCITY 선택
	}
	else
		selectedCity = rndn2(canSelectSize) + 1;	// STARTCITY을 선택 하지 않기위해 + 1

	rndCity = cityPool[selectedCity];
	if ((canSelectSize) != selectedCity) {			// 선택된 도시가 제선택 되지 않도록 배열의 끝으로 이동
		tmp = cityPool[canSelectSize];
		cityPool[canSelectSize] = cityPool[selectedCity];
		cityPool[selectedCity] = tmp;
	}
	return rndCity;									// 선택된 도시 반환
}

void mutation(int pool[POOLSIZE][DATANO]) {
	int i, j, k;								// 상호 교환될 도시순서의 인덱스
	int tmp;									// 교환 임시 저장 변수

	for (i = 0; i < POOLSIZE; ++i) {
		if (rndn(100) / 100.0 < MRATE) {	// 돌연변이가 일어난다면
			do {
				j = rndn2(DATANO - 1);			// 마지막 도시는 무조건 STARTCITY이여야 하므로 0~18값만 가짐
				k = rndn2(DATANO - 1);
			} while (i == j);

			tmp = pool[i][k];					// 임의 선택된 두 도시의 순서 교체
			pool[i][k] = pool[i][j];
			pool[i][j] = tmp;
		}
	}
}