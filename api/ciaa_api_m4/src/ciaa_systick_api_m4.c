/*
 * ciaa_systick_api_m4.c
 *
 *  Created on: Jun 7, 2022
 *      Author: fabri
 */

#include "ciaa_systick_api_m4.h"

/* Interrupt handler function pointer */
void (*systick_handler)(void) = { NULL };
/* Tracks time in ms since the SysTick was enabled */
absolute_time_t absoluteTimeMs = 0;

/*
 * 	@brief	Initialize SysTick with given period
 *
 * 	@param	us: SysTick period in microseconds
 *
 * 	@return	None
 */
void systick_init(uint32_t us) {
	/* Calculate the number of ticks to match */
	uint32_t ticks = SystemCoreClock * (us / 1E6);
    /* Configure SysTick */
	SysTick_Config(ticks);
}

/*
 *	@brief	SysTick interrupt handler
 *
 *	@param	None
 *
 *	@return	None
 */
void SysTick_Handler(void) {
	/* Increment absolute time counter by 1 every milliseccond */
	absoluteTimeMs++;
	/* Check if there is any handler and call it */
	if(*systick_handler) {
		(*systick_handler)();
	}
}
