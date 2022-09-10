/*
 * NOTE: use this command to build:
 * gcc -o filters filters.c -lm -lfftw3
 */

#include <stdio.h>
#include <malloc.h>
#include <complex.h>
#include <math.h>

/*
*****************************************************************************
* DEFINES
*****************************************************************************
*/

#define FIR_OLDER 21
#define CLEAN 0
#define MAX 10000
#define BUFFER_SIZE 1023 // 2^10 -1
#define M_PI 3.1415

/*
*****************************************************************************
* VARIABLES
*****************************************************************************
*/

enum filterType
{
	LOW_PASS,
	HIGH_PASS,
	BAND_PASS,
	BAND_STOP
};
enum windowType
{
	RECTANGULAR,
	BATLETT,
	HANNING,
	HAMMING,
	BLACKMAN
};

/*
*****************************************************************************
* FUNCTIONS PROTOTYPE
*****************************************************************************
*/

double *create1TransSinc(int windowLength, double transFreq, double sampFreq, enum filterType type);

double *createWindow(double *in, double *out, int windowLength, enum windowType type);

double *createKaiserWindow(double *in, double *out, int windowLength, double beta);
double modZeroBessel(double x);

int outputFFT(char *filename, double *window, int windowLength, double sampFreq);

double fir_filter(double rawData, double *firCoef, int currentIndex, double *circularBuff, int windowLength);

/*
*****************************************************************************
* FUNCTIONS IMPLEMENT
*****************************************************************************
*/

// Do cac trong so la doi xunng nen ta chi dung mot mang co kich thuoc (windowLength -1)/2 de luu.
double *create1TransSinc(int windowLength, double transFreq, double sampFreq, enum filterType type)
{
	// WindowLeng phai la so le
	if (windowLength % 2 == 0)
	{
		fprintf(stderr, "create1TransSinc: The windowLength can not is odd");
		return NULL;
	}

	// Allocate memory to store coefficient.
	int arraySize = (windowLength + 1) / 2;
	double *coefficient = (double *)malloc((arraySize + 1) * sizeof(double));
	if (coefficient == NULL)
	{
		fprintf(stderr, "create1TransSinc: Cannot allocate memory for coefficient");
		return NULL;
	}

	// Check type type
	if (type != LOW_PASS && type != HIGH_PASS)
	{
		fprintf(stderr, "create1TransSinc: The type must be LOW_PASS or HIGH_PASS");
		return NULL;
	}

	// Caculate coefficent of filter
	int k;
	double val = 2.0 * (double)transFreq / (double)sampFreq;

	if (type == HIGH_PASS)
		val = 1.0 - val;

	coefficient[0] = val;

	for (k = 1; k < arraySize; k++)
	{
		val = sin(M_PI * k * coefficient[0]) / (M_PI * k);
		if (type == HIGH_PASS)
			val = -val;
		coefficient[k] = val;
	}

	return coefficient;
}

double *createWindow(double *in, double *out, int windowLength, enum windowType type)
{
	// WindowLeng phai la so le
	if (windowLength % 2 == 0)
	{
		fprintf(stderr, "create1TransSinc: The windowLength can not is odd");
		return NULL;
	}

	// Allocate memory to store window.
	int arraySize = (windowLength + 1) / 2;
	if (out == NULL)
	{
		out = (double *)malloc((arraySize + 1) * sizeof(double));
		if (out == NULL)
		{
			fprintf(stderr, "create1TransSinc: Cannot allocate memory for coefficient");
			return NULL;
		}
	}

	// Luu y: Trong cong thuc vi du n = [0,M], trong nay bien chay k = [-M/2,M/2].
	// Bien doi qua lai 2 bien bang cat dat k = n - M/2 va thay vao cong thuc.

	int M = windowLength - 1;
	int halfLength = (M - 1) / 2;
	int k;
	switch (type)
	{
	case RECTANGULAR:
		for (k = 0; k < arraySize; k++)
		{
			out[k] = 1.0;
		}
		break;

	case BATLETT:
		for (k = 0; k < arraySize; k++)
		{
			out[k] = 1.0 - 2.0 * k / M;
		}
		break;

	case HANNING:
		for (k = 0; k < arraySize; k++)
		{
			out[k] = 0.5 + 0.5 * cos(M_PI * 2 * k / M);
		}
		break;

	case HAMMING:
		for (k = 0; k < arraySize; k++)
		{
			out[k] = 0.54 + 0.46 * cos(M_PI * 2 * k / M);
		}
		break;

	case BLACKMAN:
		for (k = 0; k < arraySize; k++)
		{
			out[k] = 0.42 + 0.5 * cos(M_PI * 2 * k / M) + 0.08 * cos(4 * M_PI * k / M);
		}
		break;
	default:
		fprintf(stderr, "create1TransSinc: Bad type of window");
		return NULL;
	}

	// If input has been given, multiply with out
	if (in != NULL)
	{
		for (k = 0; k < arraySize; k++)
		{
			out[k] *= in[k];
		}
	}
	return out;
}

/*
 * Function	: createKaiserWindow
 * Propose	: calculate param for Keise windows
 * Param 	: Transition Width (transWidth) is given in Hz
			  Sampling Frequency (sampFreq) is given in Hz
			  Window Length (windowLength) will be set
 * Return	: Pointer to array store coefficient of windows
 */
void calulateKaiserParam(int *M, double *beta, double ripple, double transWidth, double transFreq)
{
	double A = -20.0 * log10(ripple);
	double tw = 2 * M_PI * transWidth / transFreq;

	if (A <= 21)
	{
		*M = ceil(5.79 / tw);
		*beta = 0;
	}
	if (A > 21 && A <= 50)
	{
		*M = ceil((A - 7.95) / (2.285 * tw));
		*beta = 0.5842 * pow((A - 21), 0.4) + 0.07886 * (A - 21);
	}
	else
	{
		*M = ceil((A - 7.95) / (2.285 * tw));
		*beta = 0.1102 * (A - 8.7);
	}
}

double modZeroBessel(double x)
{
	int i;

	double x_2 = x / 2;
	double num = 1;
	double fact = 1;
	double result = 1;

	for (i = 1; i < 20; i++)
	{
		num *= x_2 * x_2;
		fact *= i;
		result += num / (fact * fact);
		//	printf("%f %f %f\n", num, fact, result);
	}

	return result;
}

double *createKaiserWindow(double *in, double *out, int windowLength, double beta)
{
	int arraySize = (windowLength - 1) / 2;
	int M = windowLength - 1;
	int k;

	if (out == NULL)
	{
		out = (double *)malloc(arraySize * sizeof(double));
		if (out == NULL)
		{
			fprintf(stderr, "createKaiserWindow: Can not allocate memory for window\n");
			return NULL;
		}
	}

	// caculate param
	for (k = 0; k < arraySize; k++)
	{
		out[k] = modZeroBessel(beta * sqrt(1 - (4 * k * k) / (M * M))) / modZeroBessel(beta);
	}

	if (in != NULL)
	{
		for (k = 0; k < arraySize; k++)
			out[k] *= in[k];
	}

	return out;
}

/*
 * Function	: fir_filter
 * Propose	: Calculate data before use FFR filter
 * Param 	: rawData, firCoef, bReset
 * Result 	: data filted.
 */
double fir_filter(double rawData, double *firCoef, int currentIndex, double *circularBuff, int windowLength)
{
	double dataFilted = 0.0f;
	int M = (windowLength - 1) / 2;

	dataFilted += firCoef[0] * circularBuff[(currentIndex - M + 1) & BUFFER_SIZE]; // (currentIndex + M)&BUFF_SIZE la phan tu trung tam
	int k;
	for (k = 1; k <= M; k++)
	{
		dataFilted += firCoef[k] * (circularBuff[(currentIndex - M + 1 - k) & BUFFER_SIZE] + circularBuff[(currentIndex - M + 1 + k) & BUFFER_SIZE]);
	}
	return dataFilted;
}

int main()
{
	/* LOW_PASS param*/
	double sampFreq = 8000.0;
	double transFreq = 100.0;
	int windowLength = 21;

	/* Kaiser Param*/
	double ripple = 0.01;
	double transWidth = 100.0;
	int windowLengthKaser;
	double beta;

	/* Caculate Coeffice FIR param*/
	int halfLength = (windowLength - 1) / 2;
	double *lpf = create1TransSinc(windowLength, transFreq, sampFreq, LOW_PASS);
	double *lpf_hamming = createWindow(lpf, NULL, windowLength, HAMMING);

	/* test */
	// int k;
	// for (k = 0; k <= halfLength; k++)
	// {
	// 	printf("%6.5f \t %6.5f\n", lpf[k], lpf_hamming[k]);
	// }

	/* Test data*/
	double t;
	int count = 0;
	double circularBuff[BUFFER_SIZE + 1] = {0.0f};
	double dataFilted[10000];
	double dataFilted2[10000];
	int currentIndex = 0;

	for (t = 0.0; t < 0.1; t += 1.0 / sampFreq)
	{
		circularBuff[currentIndex] = 50 * sin(2 * M_PI * 100 * t) + rand() % 15 * cos(2 * M_PI * 200 * t) + rand() % 10 * cos(2 * M_PI * 300 * t + 3 * M_PI / 4.0);
		dataFilted[count] = fir_filter(circularBuff[currentIndex], lpf, currentIndex, circularBuff, windowLength);
		dataFilted2[count] = fir_filter(circularBuff[currentIndex], lpf_hamming, currentIndex, circularBuff, windowLength);

		// test
		printf("%6.4f \t %6.4f \t %6.4f %f6.4\n", t, circularBuff[currentIndex], dataFilted[count], dataFilted2[count]);

		count++;
		currentIndex++;
		currentIndex &= BUFFER_SIZE;
	}

	return 0;
}