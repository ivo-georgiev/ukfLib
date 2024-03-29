#define _GNU_SOURCE 1
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "ukfCfg.h"
#include "ukfCfg1.h"

#define cfg0 (int)0
#define cfg1 (int)1

/*---------------------------------------------*/
/*           Function Prototype                */
/*---------------------------------------------*/
void show_matrix_obj(Matrix_t A);
void show_matrix(double *A, int n, int m);
int ukf_test(void);
int mtxlib_test(void);
int mtxlib_test2(void);
void report_compiler(void);
FILE *stream_log;

int verbose = 0;

int main(int argc, char *argv[])
{
	char *fname = NULL;
	int error = 0;
	int opt, index;

	while ((opt = getopt(argc, argv, "vo::")) != -1)
		switch(opt)
		{
			case 'v': verbose++; break;
			case 'o':
					  fname = (optarg)?strdupa(optarg):strdupa(argv[optind++]);
					  break;
			default:
					  fprintf(stderr, "Usage: %s [-v] [-ofile]\n", argv[0]);
					  exit(EXIT_FAILURE);
		}
	for (index = optind; index < argc; index++)
		printf ("Non-option argument %s\n", argv[index]);
	if (!fname)
		fname = strdupa("results.txt");

	if ((strlen(fname) < 1) || (stream_log = fopen(fname, "wt")) == 0)
	{
		perror(fname);
		stream_log = stderr;
	}
	if (verbose)
	{
		printf("verbose = %d\n", verbose);
		printf("output: %s\n", fname);
		report_compiler();
	}

	// generic matrix operation test
	error |= mtxlib_test();

	error |= mtxlib_test2();

	// UKF test start here
	error |= ukf_test();

	fclose(stream_log);

	return 0;
}

void show_matrix_obj(Matrix_t A)
{
	int i, j;

	for (i = 0; i < A.nrow; i++)
	{
		for (j = 0; j < A.ncol; j++)
		{
			fprintf(stream_log, "%2.14f ", A.val[A.ncol * i + j]);
		}
		fprintf(stream_log, "\n");
	}
	fprintf(stream_log, "\n");
}

void show_matrix(double *A, int n, int m)
{
	int i, j;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			fprintf(stream_log, "%2.14f ", A[m * i + j]);
		}
		fprintf(stream_log, "\n");
	}
	fprintf(stream_log, "\n");
}
/**
 *
 * Initialize and test UKF C implementation against expected result. Filter is tested in the loop from 15 steps.
 * Total root square error is accumulated in the same loop for each state in order to show deviation from reference matlab solution.
 */
int ukf_test(void)
{
	_Bool tfInitCfg0 = 0;
	_Bool tfInitCfg1 = 0;
	UKF_t ukfIo[2];
	uint32_t simLoop;

	// UKF filter measurement input(data log is generated in matlab and used for UKF simulation for 15 iteration)
	static const double yt[2][15] = { { 0, 16.085992708563385, 12.714829185978214, 14.528500994457660, 19.105561355310275, 23.252820029388918,
		29.282949862903255, 36.270058819651275, 44.244884173240955, 47.394243121124411, 55.988905459180458,
		61.667450941562109, 68.624980301613647, 76.337963872393104, 82.611325690835159 }, // y1 test
				 { 0, 16.750821420874981, 14.277640835870006, 16.320754051600520, 20.560460303503849, 24.827446289454556, 31.290961393448615, 36.853553457560210,
					 42.157283183453522, 49.382835230961490, 57.516319669684677, 65.664496283509095, 71.428712755732704, 79.241720894223079, 84.902760328915676 } };

	// UKF filter expected system states calculated with matlab script for 15 iterations
	static const double x_exp[15][4] = {
/*          x1                   x2                   x3                   x4*/
		{ 4.901482729572258, 4.576939885855807, 49.990342921246459, 49.958134463327802 },
		{ 10.103304943868373, 9.409135720815829, 50.226544716205318, 49.750795004242228 },
		{ 15.132069573131298, 14.138974122835807, 50.429540890147599, 49.191128327737864 },
		{ 20.322823824348411, 19.096919763991380, 50.836860772439010, 49.189580207886742 },
		{ 24.940146120267713, 23.399758647105461, 49.577386595072561, 47.383813382660449 },
		{ 30.021901202161214, 27.882145089050120, 49.977568320123794, 46.551744562547626 },
		{ 34.844137036519108, 32.753891693435087, 49.474006205027358, 47.190693993547214 },
		{ 39.048783329419251, 38.499098203031146, 47.606199725375902, 50.001113730363919 },
		{ 43.883085498256158, 42.383331307689538, 47.209657232695072, 46.747611757031784 },
		{ 49.479941190207498, 47.255980559687778, 49.911944505272395, 47.887236233284476 },
		{ 55.928745858553086, 51.180270357916882, 53.472542964944132, 46.510558543249353 },
		{ 61.636426955126616, 55.275415649334157, 54.052126522632797, 45.262815392265203 },
		{ 67.755622369652016, 59.602096868661732, 55.881393486598796, 45.326766104509289 },
		{ 73.045763444967164, 63.838187852739992, 54.782159791340007, 44.291415099856643 },
		{ 80.489525793047093, 66.908477563332085, 58.973616985147245, 42.638148924845950 }
	};

	// UKF initialization: CFG0
	tfInitCfg0 = ukf_init(&ukfIo[cfg0], &UkfMatrixCfg0);

	if (tfInitCfg0 == 0)
	{
		double err[4] = { 0, 0, 0, 0 };
		double absErrAccum[4] = { 0, 0, 0, 0 };

		// UKF simulation CFG0: BEGIN
		for (simLoop = 1; simLoop < 15; simLoop++)
		{
			double *const py_cfg0 = ukfIo[cfg0].input.y.val;

			// UKF:CFG0 apply/load system measurements in working array for current iteration.
			py_cfg0[0] = yt[0][simLoop];
			py_cfg0[1] = yt[1][simLoop];

			// UKF:CFG0 periodic task call
			(void)ukf_step(&ukfIo[cfg0]);

			err[0] = fabs(ukfIo[cfg0].update.x.val[0] - x_exp[simLoop - 1][0]);
			err[1] = fabs(ukfIo[cfg0].update.x.val[1] - x_exp[simLoop - 1][1]);
			err[2] = fabs(ukfIo[cfg0].update.x.val[2] - x_exp[simLoop - 1][2]);
			err[3] = fabs(ukfIo[cfg0].update.x.val[3] - x_exp[simLoop - 1][3]);

			fprintf(stream_log, "Loop: %d |system states : ukf.m | system states : est | system states : impl. diff \n", (int)simLoop);
			fprintf(stream_log, "          %2.14f        %2.14f       %2.14f\n", x_exp[simLoop - 1][0], ukfIo[cfg0].update.x.val[0], err[0]);
			fprintf(stream_log, "          %2.14f        %2.14f       %2.14f\n", x_exp[simLoop - 1][1], ukfIo[cfg0].update.x.val[1], err[1]);
			fprintf(stream_log, "          %2.14f        %2.14f       %2.14f\n", x_exp[simLoop - 1][2], ukfIo[cfg0].update.x.val[2], err[2]);
			fprintf(stream_log, "          %2.14f        %2.14f       %2.14f\n", x_exp[simLoop - 1][3], ukfIo[cfg0].update.x.val[3], err[3]);

			// accumulate the differennce between reference matlab implementation and results from C code execution
			absErrAccum[0] += err[0];
			absErrAccum[1] += err[1];
			absErrAccum[2] += err[2];
			absErrAccum[3] += err[3];
		}
		fprintf(stream_log, "Accumulated error: CFG0 \n");
		fprintf(stream_log, "%2.16f  \n%2.16f  \n%2.16f  \n%2.16f \n", absErrAccum[0], absErrAccum[1], absErrAccum[2], absErrAccum[3]);

		// UKF simulation CFG0: END
	}
	else
	{
		// initialization fail
	}

	// UKF initialization: CFG1(free pendulum)
	tfInitCfg1 = ukf_init(&ukfIo[cfg1], &UkfMatrixCfg1);

	if (tfInitCfg1 == 0)
	{
		static double tetha = 0.5;   // initial conditions for angle
		static double tetha_dot = 0; // initial conditions for angle speed
		const double B = 0.05;		  // kg*s/m
		const double l = 0.613;
		const double m = 0.5;
		const double g = 9.81;
		static const double T0 = 0.0001;

		// UKF simulation: BEGIN
		for (simLoop = 0; simLoop < 70; simLoop++)
		{
			double *const py_cfg1 = ukfIo[cfg1].input.y.val;
			double err[2] = { 0, 0 };

			// UKF:CFG1 apply/load system measurements in working array for current iteration

			tetha = tetha + T0 * tetha_dot;
			tetha_dot = tetha_dot - ((T0 * B * tetha_dot) / m) - ((T0 * g) / l) * sin(tetha);

			py_cfg1[0] = tetha;

			// UKF:CFG0 periodic task call
			(void)ukf_step(&ukfIo[cfg1]);

			err[0] = fabs(ukfIo[cfg1].update.x.val[0] - tetha);
			err[1] = fabs(ukfIo[cfg1].update.x.val[1] - tetha_dot);

			fprintf(stream_log, "Loop: %d |system states : real | system states : est | system states : err \n", (int)simLoop);
			fprintf(stream_log, "          %2.14f       %2.14f      %2.14f\n", tetha, ukfIo[1].update.x.val[0], err[0]);
			fprintf(stream_log, "          %2.14f      %2.14f     %2.14f\n", tetha_dot, ukfIo[1].update.x.val[1], err[1]);
		}
		// UKF simulation: END
	}
	else
	{
		// initialization fail
		// TBD
		return -1;
	}
	return 0;
}
/**
 *
 * Test some generic matrix operations from mtxLib.c
 *
 */
int mtxlib_test(void)
{
	/*------------------------------------------------*/
	/*             Sample Matrix not related with UKF */
	/*------------------------------------------------*/
	static double symMtx[5][5] = { { 3.009964742796415, -1.009481719676593, 1.161774056277037, -0.291928007896218, -0.775215823251770 },
		{ -1.009481719676593, 3.009964742796415, -1.009481719676593, 1.161774056277037, -0.291928007896218 },
		{ 1.161774056277037, -1.009481719676593, 3.009964742796415, -1.009481719676593, 1.161774056277037 },
		{ -0.291928007896218, 1.161774056277037, -1.009481719676593, 3.009964742796415, -1.009481719676593 },
		{ -0.775215823251770, -0.291928007896218, 1.161774056277037, -1.009481719676593, 3.009964742796415 } };

	static double symMtxChol[5][5] = { { 1.734924996302842, -0.581858997840148, 0.669639355449256, -0.168265491890613, -0.446829589119858 },
		{ 0, 1.634443284249678, -0.379239855780691, 0.650904975441148, -0.337680621986338 }, { 0, 0, 1.554903536627710, -0.418003689501540, 0.857240820764834 },
		{ 0, 0, 0, 1.543776893059448, -0.328117294491480 }, { 0, 0, 0, 0, 1.361527478565284 } };

	static double Identity_5x5[4][4] = {
		{ 1.0, 0, 0, 0, },
		{ 0, 1.0, 0, 0, },
		{ 0, 0, 1.0, 0, },
		{ 0, 0, 0, 1.0 } };

	static double TestMatrix_0_4x4[4][4] = { { 3.0, 5.0, -1.0, -4 }, { 1.0, 4.0, -0.7, -3 }, { 0, -2.0, 0, 1 }, { -2.0, 6.0, 0, 0.3 } };

#if 0
	static double TestMatrix_1_3x3[3][3] =
	{{10.5, 2.17, 3.03},
		{ 0.44, 0.59, 6.89},
		{ 7.56, 8.17, 9.21}};

	static double TestMatrix_2_3x3[3][3] =
	{ {1.11, 29.3, 31.2},
		{45.3, 5.17, 6.11},
		{7.61, 88.0, 9.34}};
#endif

	static double TestMatrix_1_2x3[2][3] = { { 1.11, 29.3, 31.2 }, // size 3x3
		{ 45.3, 5.17, 6.11 } };

	static double TestMatrixDest_3x2[3][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } };
	Matrix_t myFactMatrix;
	Matrix_t myTestMatx = { 0, 0, 0, NULL };
	Matrix_t myChol = { 0, 0, 0, NULL };
	Matrix_t Im = { 0, 0, 0, NULL };
	Matrix_t oTestMatrix_0_4x4 = { 0, 0, 0, NULL };
	Matrix_t oTestMatrixDest_3x2 = { 0, 0, 0, NULL };

	mtx_init(&myTestMatx, &TestMatrix_1_2x3[0][0], NROWS(TestMatrix_1_2x3), NCOL(TestMatrix_1_2x3), COLXROW(TestMatrix_1_2x3));
	mtx_init(&myChol, &symMtxChol[0][0], NROWS(symMtxChol), NCOL(symMtxChol), COLXROW(symMtxChol));
	mtx_init(&Im, &Identity_5x5[0][0], NROWS(Identity_5x5), NCOL(Identity_5x5), COLXROW(Identity_5x5));
	mtx_init(&oTestMatrix_0_4x4, &TestMatrix_0_4x4[0][0], NROWS(TestMatrix_0_4x4), NCOL(TestMatrix_0_4x4), COLXROW(TestMatrix_0_4x4));
	mtx_init(&oTestMatrixDest_3x2, &TestMatrixDest_3x2[0][0], NROWS(TestMatrixDest_3x2), NCOL(TestMatrixDest_3x2), COLXROW(TestMatrixDest_3x2));
	// show_matrix(&TestMatrix_1_2x3[0][0],2,3);
	// show_matrix_obj(myTestMatx);

	(void)mtx_init(&myFactMatrix, &symMtx[0][0], NROWS(symMtx), NCOL(symMtx), COLXROW(symMtx));
	// show_matrix(&symMtx[0][0],5,5);

	(void)mtx_chol_lower(&myFactMatrix);
	if (verbose > 1) show_matrix_obj(myFactMatrix);

	/*show_matrix_obj(myChol);
	  mtx_transp_square(&myChol);
	  show_matrix_obj(myChol);
	  */

	if (verbose > 1) show_matrix_obj(Im);
	if (verbose > 1) show_matrix_obj(oTestMatrix_0_4x4);
	mtx_inv(&oTestMatrix_0_4x4, &Im);
	if (verbose > 1) show_matrix_obj(Im);

	mtx_identity(&Im);
	if (verbose > 1) show_matrix_obj(Im);

	// show_matrix_obj(myTestMatx);
	// mtx_transp_dest(&myTestMatx,&oTestMatrixDest_3x2);
	// show_matrix_obj(oTestMatrixDest_3x2);*/
	return 0; // OK
}

/**
 * prints the size of the floating point types
 *
 */
void report_compiler(void)
{
	fprintf(stderr, "sizeof float = %d bits\nsizeof double = %d bits\nsizeof long double = %d bits\n", 8 * __SIZEOF_FLOAT__, 8 * __SIZEOF_DOUBLE__,
			8 * __SIZEOF_LONG_DOUBLE__);
}

void test_mtx_cpy(Matrix_t const *const a, Matrix_t *const b);
double distance1(Matrix_t *a, Matrix_t *b);
double distance2(Matrix_t *a, Matrix_t *b);
double distance_inf(Matrix_t *a, Matrix_t *b);

#include "mtxGen.c"

double m_temp[12][12];
Matrix_t o_temp={12*12,12,12,(double*) m_temp};
double m_temp2[12][12];
Matrix_t o_temp2={12*12,12,12,(double*) m_temp2};

int mtxlib_test2(void)
{
	int i;
	int n = (int)sizeof(allSym)/(int)sizeof(*allSym);
	int res = 0;

	for (i=0;i<n;i++)
	{
		// allSym[i]; loCholSym[i]; invSym[i];
		res |= !(allSym[i]->nelem == loCholSym[i]->nelem && invSym[i]->nelem == loCholSym[i]->nelem);
		o_temp.nelem = allSym[i]->nelem;
		o_temp.nrow = allSym[i]->nrow;
		o_temp.ncol = allSym[i]->ncol;
		mtx_identity(&o_temp);
		test_mtx_cpy(allSym[i], &o_temp2);
		res |= (MTX_OPERATION_OK != mtx_inv(&o_temp2, &o_temp));
		double d1 = distance1(invSym[i], &o_temp);
		double d2 = distance2(invSym[i], &o_temp);
		double d3 = distance_inf(invSym[i], &o_temp);
		// fprintf(stderr, "%3d d = %G\n", i, d);
		fprintf(stream_log, "inv: %3d d1 = %G\td2 = %G\td3 = %G\n", i, d1, d2, d3);
	}

	for (i=0;i<n;i++)
	{
		test_mtx_cpy(allSym[i], &o_temp);

		res |= (MTX_OPERATION_OK != mtx_chol_lower(&o_temp));
		double d1 = distance1(loCholSym[i], &o_temp);
		double d2 = distance2(loCholSym[i], &o_temp);
		double d3 = distance_inf(loCholSym[i], &o_temp);
		fprintf(stream_log, "chL: %3d d1 = %G\td2 = %G\td3 = %G\n", i, d1, d2, d3);
		if (verbose > 1)
		{
			show_matrix_obj(*loCholSym[i]);
			show_matrix_obj(o_temp);
		}
	}

	for (i=0;i<n;i++)
	{
		test_mtx_cpy(allSym[i], &o_temp);

		res |= (MTX_OPERATION_OK != mtx_chol_upper(&o_temp));
		double d1 = distance1(upCholSym[i], &o_temp);
		double d2 = distance2(upCholSym[i], &o_temp);
		double d3 = distance_inf(upCholSym[i], &o_temp);
		fprintf(stream_log, "chU: %3d d1 = %G\td2 = %G\td3 = %G\n", i, d1, d2, d3);
		if (verbose > 1)
		{
			show_matrix_obj(*upCholSym[i]);
			show_matrix_obj(o_temp);
		}
	}
	n = sizeof(allMagic)/sizeof(*allMagic);

	for (i=0;i<n;i++)
	{
		test_mtx_cpy(allMagic[i], &o_temp);

		res |= (MTX_OPERATION_OK != mtx_transp_square(&o_temp));
		double d1 = distance1(transpMagic[i], &o_temp);
		double d2 = distance2(transpMagic[i], &o_temp);
		double d3 = distance_inf(transpMagic[i], &o_temp);
		fprintf(stream_log, "trn: %3d d1 = %G\td2 = %G\td3 = %G\n", i, d1, d2, d3);
		if (verbose > 1)
		{
			show_matrix_obj(*upCholSym[i]);
			show_matrix_obj(o_temp);
		}
	}

	return res;
}

void test_mtx_cpy(Matrix_t const *const a, Matrix_t *const b)
{
	int i, j;

	b->ncol = a->ncol;
	b->nrow = a->nrow;
	b->nelem = a->nelem;

	for (i=0;i<a->nrow;i++)
	{
		int o = a->ncol * i;
		for (j=0;j<a->ncol;j++)
		{
			b->val[o + j] = a->val[o + j];
		}
	}
}

double distance1(Matrix_t *a, Matrix_t *b)
{
	double d = 0;
	int i, j;

	for (i=0;i<a->nrow;i++)
	{
		int o = a->ncol * i;
		for (j=0;j<a->ncol;j++)
		{
			double aa = a->val[o + j];
			double bb = b->val[o + j];

			d += fabs(aa - bb);
		}
	}
	return d;
}


double distance2(Matrix_t *a, Matrix_t *b)
{
	double d = 0;
	int i, j;

	for (i=0;i<a->nrow;i++)
	{
		int o = a->ncol * i;
		for (j=0;j<a->ncol;j++)
		{
			double aa = a->val[o + j];
			double bb = b->val[o + j];
			double t = (aa - bb);
			d += t*t;
		}
	}
	return sqrt(d);
}


double distance_inf(Matrix_t *a, Matrix_t *b)
{
	double d = 0;
	int i, j;

	for (i=0;i<a->nrow;i++)
	{
		int o = a->ncol * i;
		for (j=0;j<a->ncol;j++)
		{
			double aa = a->val[o + j];
			double bb = b->val[o + j];
			double t = fabs(aa - bb);
			d = fmax(d,t);
		}
	}
	return sqrt(d);
}

