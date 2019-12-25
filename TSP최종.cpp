/* 2015-11-25 ~ 2015-11-26 �迵��
�ΰ����� ������ �˰��� - ��ȣ ��ȯ �������� ���
��ȸ�湮������
���� 1���� �����Ͽ� ��� ���ø� �湮 ��
���� 1�� ��ȯ�ϴ� ���� ª�� ��� ���ϱ� */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <share.h>

#define DATANO		30		// ������ ��
#define POOLSIZE	30		// Ǯ ������
#define LASTG		600		// �ߴ� ����
#define MRATE		0.01	// �������� Ȯ��
#define STARTCITY	1		// ���� �� ���ƿ� ���� ��ȣ
#define YES			1
#define NO			0
#define FILEPATH	"..\\result_GA.txt"	// ���� ���, �� ������Ʈ ������ �����Ѵ�.

//char pathbuf[100];			// Bestfit�� ��� ���� ��θ� ����ϱ� ���� buffer
FILE *fp;						// ����� ���Ϸ� �����ϱ� ���� ����� ���� ������

double calcDistance(double x1, double x2, double y1, double y2);
double evalfit(int gene[]);					// ������ ���
void mating(int pool[POOLSIZE][DATANO]);	// ����
void mutation(int pool[POOLSIZE][DATANO]);	// ��������
void printp(int pool[POOLSIZE][DATANO], int generation);	// ��� ���
void initpool(int pool[POOLSIZE][DATANO]);	// �ʱ� ���� ���� 

double rndn(double l);								// n �̸��� ���� ���� double
int rndn2(int l);									// n �̸��� ���� ���� int
int rndnInit(int cityPool[], int canSelectSize);	// �ʱ� ���� ������ �ұ�Ģ�ϰ� ���� ����
int select(double roulette[POOLSIZE], double totalfitness);	// �θ� ����
void copypool(int pool[POOLSIZE][DATANO], int nextpool[POOLSIZE][DATANO]);	// ���� ������ Ǯ ����
void crossing(int m[], int p[], int c1[], int c2[]);	// ������ ����

double q[][2] = {
	3.9, 2.1, 9.7, 7.1, 8.2, 6.8, 2.9, 0.2, 7.2, 0.7,
	6.5, 0.3, 3.1, 7.4, 2.7, 2.2, 1.1, 3.4, 0.6, 8.8,
	7.4, 8.6, 9.4, 5.4, 2.7, 7.6, 4.1, 8.5, 6.7, 1.1,
	2.0, 3.9, 1.7, 4.4, 0.1, 9.9, 9.9, 3.9, 4.6, 1.6,
	8.4, 0.9, 0.7, 9.1, 2.2, 3.4, 6.5, 2.7, 5.6, 5.4,
	7.8, 3.3, 3.6, 6.3, 8.9, 7.1, 3.8, 2.9, 4.7, 5.2
};

double poolAver[LASTG];				// ������� ��� ���յ� ���� �����迭
double poolElite[LASTG];				// ������� ����Ʈ ���� ���� �迭
double poolSuper = 0;						// ������� ���ۿ���Ʈ ���� ���� �迭

int main() {
	int pool[POOLSIZE][DATANO];				// ������ Ǯ
	int generation;							// ���� �����
	clock_t start = clock();				// ���� ���� �ð�

	srand((unsigned int)time(NULL));		// ���� �ʱ�ȭ

	if ((fp = _fsopen(FILEPATH, "w+", _SH_DENYNO)) == NULL)
		perror("File open error");

	initpool(pool);							// �ʱ� ���� ����
	fprintf(fp, "pool[�ʱ�����] : \n");
	for (int i = 0; i < POOLSIZE; i++) {	// �ʱ����� ���
		for (int j = 0; j < DATANO; j++) {
			fprintf(fp, "%d ", pool[i][j]);
		}
		fprintf(fp, "\n\n");
	}

	for (generation = 0; generation < LASTG; ++generation) {	// �ߴ� ������� �ݺ�
		fprintf(fp, "%d����\n", generation);
		mating(pool);		// ����
		mutation(pool);		// ��������
		printp(pool, generation);		// ��� ���
	}

	fprintf(fp, "����\t\t���\t\t����Ʈ\n");
	for (generation = 0; generation < LASTG; ++generation) {	// �����Ͽ� ��� ���
		fprintf(fp, "%d����\t %lf \t", generation, poolAver[generation]);
		fprintf(fp, "%lf\n", poolElite[generation]);
	}
	fprintf(fp, "SuperElite : %lf\n", poolSuper);



	fprintf(fp, "Execution time: %.3lfs\n", (double)(clock() - start) / CLOCKS_PER_SEC);

	fclose(fp);

	return 0;
}

int select(double roulette[POOLSIZE], double totalfitness) {
	int i;				// �ݺ� ���� ����
	double ball;		// ��
	double acc = 0;		// ������ ���갪

	for (i = 0; i < POOLSIZE; ++i) {
		acc += totalfitness / roulette[i];	// �Ÿ��� �������� �����ϱ� ������ ������ ����
	}

	ball = rndn(acc);	// ������ �յ�
	acc = 0;
	for (i = 0; i < POOLSIZE; ++i) {
		acc += totalfitness / roulette[i];
		if (acc > ball)	// �����ϴ� ������
			break;
	}

	return i;
}

void mating(int pool[POOLSIZE][DATANO]) {
	int i;							// �ݺ� ���� �Լ�
	double totalfitness = 0;		// ������ �հ�
	int nextpool[POOLSIZE][DATANO];	// �ڽ� ����
	double roulette[POOLSIZE];		// ������ ����
	int mama, papa;					// �θ� ������ ��ȣ

	for (i = 0; i < POOLSIZE; ++i) {	// �귿 �ۼ�
		roulette[i] = evalfit(pool[i]);	// ������ ���

		totalfitness += roulette[i];	// ������ �հ� ���
	}

	for (i = 0; i < POOLSIZE / 2; ++i) {// ���ð� ������ �ݺ�
		do {	// �θ� ����
			mama = select(roulette, totalfitness);
			papa = select(roulette, totalfitness);
		} while (mama == papa);// �θ� �ߺ� ����
		crossing(pool[mama], pool[papa], nextpool[i * 2], nextpool[i * 2 + 1]);	// �� ������ ����
	}
	copypool(pool, nextpool);	// ���� ���� Ǯ ����

}

void crossing(int m[], int p[], int c1[], int c2[]) {
	int i, j;		// �ݺ� ���� ����
	int cp1, cp2;	// ������ �� 2��
	int tmp;		// ��ȯ �ӽ� ����
	c1[DATANO - 1] = STARTCITY;	// ���������� ������ STARTCITY
	c2[DATANO - 1] = STARTCITY;

	cp1 = rndn2(DATANO - 1);	// ������ ���� ������ ���ô� ������ STARTCITY���� �ϹǷ� ������ ���� 0 ~ 18
	cp2 = rndn2(DATANO - 1);

	if (cp1 > cp2) {			// cp2�� �� ũ���� ���� ����
		tmp = cp1;
		cp1 = cp2;
		cp2 = tmp;
	}

	for (i = cp1; i <= cp2; i++) {	// ���õ� �κ� �����Ͽ� �ڽĿ��� ��
		c1[i] = p[i];
		c2[i] = m[i];
	}

	int isHaving;	// ������ �ִ� �����ΰ�?
	int INDEX;		// ������ ���� ���� ���ڸ� �ڽĿ��� ���� �ε���

	INDEX = 0;		// �ڽ��� ù ��°�� �ִ´�
	for (i = 0; i < DATANO - 1; i++) {	// c1���� mama�� ���� ���ڸ� ä����
		for (j = cp1; j <= cp2; j++) {
			if (c1[j] == m[i]) {
				isHaving = YES;			// �ڽ��� ������ �ִٸ� �ʿ� ����.
				break;
			}
			else if (c1[j] != m[i])		// Ȯ���� �� �������� �ڽ��� ���ڸ� ������ ���� ������
				isHaving = NO;
		}

		if (isHaving == NO) {
			if (INDEX < cp1) {			// cp1 ���ʿ� ������ ���� ���� ���ڸ� ����
				c1[INDEX] = m[i];
				INDEX++;
				if (INDEX == cp1)
					INDEX = cp2 + 1;
			}
			else if (INDEX == cp1) {	// cp1 ������ ��á�ٸ� cp2�ڷ� �ű�
				INDEX = cp2 + 1;
				c1[INDEX] = m[i];
				INDEX++;
			}
			else if (INDEX > cp2) {		// cp2 ���ʿ� ������ ���� �ʴ� ���ڸ� ����
				c1[INDEX] = m[i];
				INDEX++;
			}
		}

	}

	INDEX = 0;		// �ڽ��� ù��°�� �ִ´�
	for (i = 0; i < DATANO - 1; i++) {	// c2���� papa�� ���� ���ڸ� ä����
		for (j = cp1; j <= cp2; j++) {
			if (c2[j] == p[i]) {
				isHaving = YES;			// �ڽ��� ������ �ִٸ� �ʿ� ����.
				break;
			}
			else if (c2[j] != p[i])		// Ȯ���� �� �������� �ڽ��� ���ڸ� ������ ���� ������
				isHaving = NO;
		}

		if (isHaving == NO) {
			if (INDEX < cp1) {			// cp1 ���ʿ� ������ ���� ���� ���ڸ� ����
				c2[INDEX] = p[i];
				INDEX++;
				if (INDEX == cp1)
					INDEX = cp2 + 1;
			}
			else if (INDEX == cp1) {	// cp1 ������ ��á�ٸ� cp2�ڷ� �ű�
				INDEX = cp2 + 1;
				c2[INDEX] = p[i];
				INDEX++;
			}
			else if (INDEX > cp2) {		// cp2 ���ʿ� ������ ���� �ʴ� ���ڸ� ����
				c2[INDEX] = p[i];
				INDEX++;
			}
		}
	}
}


void copypool(int pool[POOLSIZE][DATANO], int nextpool[POOLSIZE][DATANO]) {
	int i, j;//�ݺ� ���� �Լ�

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
	int i;				// �ݺ� ���� ����
	double fitness = 0;	// ������ ���߰�

	for (i = 0; i < DATANO; i++) {
		if (i == 0)
			// ó�� �湮 ���ÿ� STARTCITY�� �Ÿ�
			fitness += calcDistance(q[g[i] - 1][0], q[STARTCITY - 1][0], q[g[i] - 1][1], q[STARTCITY - 1][1]);
		else
			// ���� �湮 ���ÿ� ���� ������ �Ÿ�
			fitness += calcDistance(q[g[i] - 1][0], q[g[i - 1] - 1][0], q[g[i] - 1][1], q[g[i - 1] - 1][1]);
	}

	return fitness;
}

void printp(int pool[POOLSIZE][DATANO], int generation) {
	int i, j;				// �ݺ� ���� ����
	double fitness;			// ������
	double totalfitness = 0;// ������ �հ�
	int elite;
	double bestfit = 0;		// ����Ʈ ������ ó���� ����

	for (i = 0; i < POOLSIZE; i++) {
		for (j = 0; j < DATANO; j++) {
			fprintf(fp, "%1d ", pool[i][j]);
		}
		fitness = evalfit(pool[i]);
		fprintf(fp, "\t%lf\n", fitness);

		if (bestfit == 0) {	//bestfit �ʱ� �� ����
			bestfit = fitness;
			elite = i;
		}
		else if (fitness < bestfit) {	// ����Ʈ ��
			bestfit = fitness;
			elite = i;
		}

		totalfitness += fitness;
	}

	fprintf(fp, "\nElite: %d\tBestfit: %lf\n", elite, bestfit);			// ����Ʈ ���� ������ ���
	fprintf(fp, "Average Fitness: %lf\n\n", totalfitness / POOLSIZE);	// ��� ������ ���

																		/*������ ��� ���������鿡 ����*/
	poolAver[generation] = totalfitness / POOLSIZE;		//�ش� ������ ��� ������ ����
	poolElite[generation] = bestfit;					//�ش� ������ ����Ʈ ������ ����
	if (poolSuper == 0)		//super ����Ʈ �ʱⰪ ����
		poolSuper = bestfit;
	else if (poolSuper > bestfit)
		poolSuper = bestfit;							//��ü ���뿡�� ã�� ���ۿ���Ʈ ����
}

void initpool(int pool[POOLSIZE][DATANO]) {
	int i, j, k;			// �ݺ� �����
	int cityPool[DATANO];	// ���� ������ ������ ����
	int canSelecSize;		// ���� ������ ���õ��� ����

	for (i = 0; i < POOLSIZE; ++i) {
		canSelecSize = DATANO - 1;		// STARTCITY�� �������� �ʱ� ���� - 1

		for (k = 0; k < DATANO; k++) {	// ���� ���� ������ ���� �ʱ�ȭ
			if (k == 0)
				cityPool[k] = STARTCITY;// ���õ��� �ʴ� �迭�� 0�� ���ڿ� ���� ���ø� ����
			cityPool[k] = k + 1;
		}

		for (j = 0; j < DATANO; ++j) {
			pool[i][j] = rndnInit(cityPool, canSelecSize);	// �������� ������ ���� ������ �����ϰ� ����
			canSelecSize--;									// ���õ� ���� �� ���� ���� �ʵ��� �������� ����
		}
	}
}

double rndn(double l) {
	double rndno;	// ������ ����

	while ((rndno = ((double)rand() / RAND_MAX) * l) == l);// ���ڰ� �̸��� �������� 

	return rndno;
}

int rndn2(int l) {
	int rndno;		// ������ ����

	while ((rndno = ((double)rand() / RAND_MAX) * l) == l);// ���ڰ� �̸��� �������� 

	return rndno;
}

int rndnInit(int cityPool[], int canSelectSize) {
	int rndCity;		// ������ ����
	int tmp;			// �迭�� ��ġ�� �ٲٱ� ���� �ӽ� ����
	int selectedCity;	// ���õ� ����, 

	if (canSelectSize == 0) {						// ���� ���ð� �Ѱ� ���Ҵٸ�
		selectedCity = STARTCITY - 1;
		return rndCity = cityPool[selectedCity];	// STARTCITY ����
	}
	else
		selectedCity = rndn2(canSelectSize) + 1;	// STARTCITY�� ���� ���� �ʱ����� + 1

	rndCity = cityPool[selectedCity];
	if ((canSelectSize) != selectedCity) {			// ���õ� ���ð� ������ ���� �ʵ��� �迭�� ������ �̵�
		tmp = cityPool[canSelectSize];
		cityPool[canSelectSize] = cityPool[selectedCity];
		cityPool[selectedCity] = tmp;
	}
	return rndCity;									// ���õ� ���� ��ȯ
}

void mutation(int pool[POOLSIZE][DATANO]) {
	int i, j, k;								// ��ȣ ��ȯ�� ���ü����� �ε���
	int tmp;									// ��ȯ �ӽ� ���� ����

	for (i = 0; i < POOLSIZE; ++i) {
		if (rndn(100) / 100.0 < MRATE) {	// �������̰� �Ͼ�ٸ�
			do {
				j = rndn2(DATANO - 1);			// ������ ���ô� ������ STARTCITY�̿��� �ϹǷ� 0~18���� ����
				k = rndn2(DATANO - 1);
			} while (i == j);

			tmp = pool[i][k];					// ���� ���õ� �� ������ ���� ��ü
			pool[i][k] = pool[i][j];
			pool[i][j] = tmp;
		}
	}
}