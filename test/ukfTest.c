/**
 * GNU General Public License v3
 * Copyright (c) 2017 ivo-georgiev
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <string.h>

#include "ukfCfg.h"
#include "ukfCfg1.h"

#define cfg0 (int)0
#define cfg1 (int)1

/*---------------------------------------------*/
/*           Function Prototype                */
/*---------------------------------------------*/
int ukf_test(void);
FILE *stream_log;

int main(int argc, char *argv[])
{
    // Prevent compiler warnings for unused parameters
    (void)argc;
    (void)argv;

    int error = 0;

    // Directly open the results file using a safe, standard C approach
    stream_log = fopen("results.txt", "wt");
    if (stream_log == NULL)
    {
        perror("Failed to open results.txt");
        stream_log = stderr;
    }

    // Inform user that the standalone UKF validation framework is starting
    printf("Starting UKF validation via MATLAB-generated flight profiles...\n");

    // Execute the core Unscented Kalman Filter test scenario
    error = ukf_test();

    // Safely close the logging file stream if it was successfully opened
    if (stream_log != stderr)
    {
        fclose(stream_log);
    }

    return error;
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





