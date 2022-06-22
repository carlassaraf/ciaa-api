/*
 * fft.h
 *
 *  Created on: Jun 22, 2022
 *      Author: fabri
 */

#ifndef FFT_H_
#define FFT_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

/* Custom constants */
typedef enum {
	FFT_FORWARD,
	FFT_REVERSE
} fft_dir_t;

/* Mathematical constants*/
#define PI			3.14159265358979323846
#define TWO_PI 		2 * PI
#define FOUR_PI 	4 * PI
#define SIX_PI 		6 * PI

#define sq(x) 	((x)*(x))

/* Windowing type */
typedef enum {
	FFT_WIN_TYP_RECTANGLE, 			/* Rectangle (Box car) */
	FFT_WIN_TYP_HAMMING, 			/* Hamming */
	FFT_WIN_TYP_HANN, 				/* Hann */
	FFT_WIN_TYP_TRIANGLE,			/* Triangle (Bartlett) */
	FFT_WIN_TYP_NUTTALL,			/* Nuttall */
	FFT_WIN_TYP_BLACKMAN,			/* Blackman */
	FFT_WIN_TYP_BLACKMAN_NUTTALL,	/* Blackman nuttall */
	FFT_WIN_TYP_BLACKMAN_HARRIS,	/* Blackman harris*/
	FFT_WIN_TYP_FLT_TOP,			/* Flat top */
	FFT_WIN_TYP_WELCH 				/* Welch */
} fft_windowing_type_t;

/* Prototype functions */
void fft_init(float *vReal, float *vImag, uint16_t samples, float samplingFrequency);
void fft_windowing(fft_windowing_type_t windowType, uint8_t dir);
void fft_compute(uint8_t dir);
void fft_complex_to_magnitude(void);
void fft_major_peak(float* mag_out, float* freq_out, float magFact);

#endif /* FFT_H_ */
