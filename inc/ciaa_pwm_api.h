/*
 * ciaa_pwm_api.h
 *
 *  Created on: Jun 5, 2022
 *      Author: fabri
 */

#ifndef CIAA_PWM_API_H_
#define CIAA_PWM_API_H_

#include "chip.h"
#include "ciaa_gpio_api.h"

#define	CIAA_PWM	LPC_SCT

typedef struct {
	pin_t *pin;
	uint8_t index;
	uint8_t pwm;
} pwm_output_t;

typedef struct {
	uint32_t frequency;
	uint8_t clkdiv;
} pwm_config_t;

typedef enum {
	CTOUT_0, CTOUT_1, CTOUT_2, CTOUT_3,
	CTOUT_4, CTOUT_5, CTOUT_6, CTOUT_7,
	CTOUT_8, CTOUT_9, CTOUT_10, CTOUT_11,
	CTOUT_12, CTOUT_13, CTOUT_14, CTOUT_15
} ctout_t;

extern pwm_output_t PWM0;
extern pwm_output_t PWM1;
extern pwm_output_t PWM2;
extern pwm_output_t PWM3;
extern pwm_output_t PWM4;
extern pwm_output_t PWM5;
extern pwm_output_t PWM6;

void pwm_init(void);
void pwm_config_init(pwm_config_t config);
void pwm_set_frequency(uint32_t frequency);
void pwm_set_output_pin(pwm_output_t output);

static inline pwm_config_t pwm_get_default_config(void) {
	/* Declare default configuration: 1KHz frequency, Clock Div 1 */
	pwm_config_t config = { 1000, 1 };
	return config;
}

static inline void pwm_start(void) { CIAA_PWM->CTRL_U &= ~(SCT_CTRL_HALT_L | SCT_CTRL_HALT_H); }

static inline void pwm_stop(void) {
	/* Stop the timer */
	CIAA_PWM->CTRL_U |= (SCT_CTRL_HALT_L | SCT_CTRL_HALT_H);
	/* Clear counter */
	CIAA_PWM->CTRL_U |= (SCT_CTRL_CLRCTR_L | SCT_CTRL_CLRCTR_H);
}

static inline void pwm_set_start(bool enabled) {
	if(enabled) { pwm_start(); }
	else { pwm_stop(); }
}

static inline void pwm_set_match(uint8_t match, uint32_t value) { CIAA_PWM->MATCHREL[match].U = value; }

static inline void pwm_set_duty_ticks(pwm_output_t output, uint32_t ticks) { pwm_set_match(output.pwm, ticks); }

static inline void pwm_set_wrap(uint32_t value) { CIAA_PWM->MATCHREL[SCT_MATCH_0].U = value; }

#endif /* CIAA_PWM_API_H_ */
