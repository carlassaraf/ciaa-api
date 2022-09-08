/*
 * ciaa_ipc_api.c
 *
 *  Created on: Aug 27, 2022
 *      Author: carlassaraf
 */

#include "ciaa_ipc_api.h"

/* Pointer to shared memory space */
static ipc_queue_t *q_ipc = (ipc_queue_t*) SHARED_MEM_IPC;

/*
 * 	@brief	Function to initialize the IPC message queue
 *
 * 	@param	data: pointer to the data to be sent
 * 	@param	size: message size
 * 	@param	count: number of items to be held in the queue
 *
 * 	@return	None, will not return if there is error in given arguments
 */
void ipc_queue_init(void *data, int size, int count) {
	/* Pointer to shared queue */
	ipc_queue_t *qrd;
	/* Sanity Check */
	if (!size || !count || !data) { while (1); }
	/* Check if size is a power of 2 */
	if (count & (count - 1)) { while (1); }
	/* Point to the shared RAM */
	qrd = q_ipc;
	/* Clear shared RAM */
	memset(qrd, 0, sizeof(*qrd));
	/* Initialize shared queue */
	qrd->count = count;
	qrd->size = size;
	qrd->data = data;
	qrd->valid = QUEUE_MAGIC_VALID;
}

/*
 * 	@brief	Function to push a message into queue with timeout
 *
 * 	@param	data: pointer to data to be pushed
 * 	@param	tout: non-zero value - timeout value in milliseconds,
 *                zero value - no blocking,
 *                negative value - blocking
 *
 * 	@return	queue_insert on success,
 * 			queue_full or queue_error on failure,
 *          queue_timeout when there is a timeout
 */
ipc_status_t ipc_push_tout(void *data, int tout) {
	/* Pointer to shared queue */
	ipc_queue_t *qwr = q_ipc;
	/* Check if queue is initialized */
	if (!ipc_queue_is_valid(qwr)) { return queue_error; }
	/* Check if queue is full */
	if ((tout == 0) && ipc_queue_is_full(qwr)) { return queue_full; }
	/* Wait for queue to be empty */
	while (ipc_queue_is_full(qwr)) {
		if (ipc_event_handler()) { return queue_timeout; }
	}
	/* Push item to queue */
	qwr->data = data;
	memcpy(qwr->data + ((qwr->head & (qwr->count - 1)) * qwr->size), data, qwr->size);
	qwr->head++;
	/* Insert was a success */
	return queue_insert;
}

/*
 * 	@brief	Function to read a message from queue with timeout
 *
 * 	@param	data: pointer to store popped data
 * 	@param	tout: non-zero value - timeout value in milliseconds,
 *                zero value - no blocking,
 *                negative value - blocking
 *
 * 	@return	queue_valid on success,
 * 			queue_empty or queue_error on failure,
 *          queue_timeout when there is a timeout
 */
ipc_status_t ipc_pop_tout(void *data, int tout) {
	/* Pointer to shared queue */
	ipc_queue_t *qrd = q_ipc;
	/* Check if queue is initialized */
	if (!ipc_queue_is_valid(qrd)) { return queue_error; }
	/* Check if queue is empty */
	if ((tout == 0) && ipc_queue_is_empty(qrd)) { return queue_empty; }
	/* Wait for queue to be full */
	while (ipc_queue_is_empty(qrd)) {
		if (ipc_event_handler()) { return queue_timeout; }
	}
	/* Pop the queue item */
	memcpy(data, qrd->data + ((qrd->tail & (qrd->count - 1)) * qrd->size), qrd->size);
	qrd->tail++;
	data = qrd->data;
	/* Pop was a success */
	return queue_valid;
}
