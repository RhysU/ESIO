/*
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * Adapted from the GNU Scientific Library
 * Copyright (C) 2010 The PECOS Development Team
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __ESIO_ERROR_H__
#define __ESIO_ERROR_H__

#include <stdio.h>

/** @file
 * Provides standardized error numbers and error handling routines.  Error
 * reporting follows the design and conventions used in the <a
 * href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a> (GSL) <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Error-Handling.html">
 * error handling routines</a>.  Much of ESIO's error code is a direct copy
 * of GSL's API and source code.  Notable exceptions are the MPI error handling
 * macros which are an improved copy of ideas found in <a
 * href="http://www.mcs.anl.gov/petsc/">PETSc</a>.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Standardized error status codes used throughout ESIO.
 * Where possible these codes are numerically equivalent to
 * <a href="http://www.gnu.org/software/gsl/manual/html_node/Error-Codes.html">
 * GSL's error codes</a>.
 */
enum esio_error_status {
    ESIO_FAILURE  = -1, /**< Failure */
    ESIO_SUCCESS  =  0, /**< Success */
    ESIO_EFAULT   =  3, /**< Invalid pointer */
    ESIO_EINVAL   =  4, /**< Invalid argument supplied by user */
    ESIO_EFAILED  =  5, /**< Generic failure */
    ESIO_ESANITY  =  7, /**< Sanity check failed - shouldn't happen */
    ESIO_ENOMEM   =  8, /**< Memory allocation failed */
    ESIO_EZERODIV = 12  /**< Tried to divide by zero */
};

/**
 * Calls the error handler last set using esio_set_error_handler
 * when invoked.  This is the entry point to the error handling system.
 *
 * The default behavior is to log the error to the stream specified using
 * esio_set_stream. The functions esio_set_stream,
 * esio_set_stream_handler, and esio_set_error_handler can be used to
 * modify this behavior.
 *
 * @param reason Reason for the error.
 * @param file File in which the error was reported.
 * @param line Line at which the error was reported.
 * @param esio_errno Error code to report.  Should be one of
 *      esio_error_status if at all possible.
 *
 * @see Most clients should not use this function directly; instead use one of
 *      the convenience macros: ESIO_ERROR, ESIO_ERROR_VAL,
 *      ESIO_ERROR_VOID, ESIO_ERROR_NULL
 */
void esio_error(const char * reason,
                const char * file,
                int line,
                int esio_errno);

/**
 * Print an error message to the current error stream.
 * If a esio_stream_handler_t has been specified, it is used.
 * If a stream has been set using esio_set_stream, it is used.
 * Lastly, the routine prints the error message to standard error.
 *
 * @param label Label used to identify the error.
 * @param file File in which the error was reported.
 * @param line Line at which the error was reported.
 * @param reason Reason for the error.
 */
void esio_stream_printf(const char *label,
                        const char *file,
                        int line,
                        const char *reason);

/**
 * Look up a human-readable error message for the given error status.
 *
 * @param esio_errno Error code to look up.
 *
 * @return A message suitable for use in logging or error messages.
 */
const char * esio_strerror(const int esio_errno);

/**
 * Defines the function prototype necessary for an error handler.
 * Error handlers should be reentrant safe if possible.
 *
 * @param reason Reason for the error.
 * @param file File in which the error was reported.
 * @param line Line at which the error was reported.
 * @param esio_errno Error code to report.
 *
 * @see esio_set_error_handler
 */
typedef void esio_error_handler_t(const char * reason,
                                  const char * file,
                                  int line,
                                  int esio_errno);

/**
 * Defines the function prototype necessary for a stream handler.
 * Stream handlers should be reentrant safe if possible.
 *
 * @param reason Reason for the error.
 * @param file File in which the error was reported.
 * @param line Line at which the error was reported.
 * @param esio_errno Error code to report.
 *
 * @see esio_set_stream_handler
 */
typedef void esio_stream_handler_t(const char * label,
                                   const char * file,
                                   int line,
                                   const char * reason);

/**
 * Sets the current error handler for the process.
 * Invoked by esio_error when an error occurs.
 *
 * @param new_handler New error handler to use.
 *
 * @return the previous error handler in use.
 */
esio_error_handler_t *
esio_set_error_handler(esio_error_handler_t * new_handler);

/**
 * An error handler implementation that disables all error reporting.
 * Primarily intended for use in test environments.
 *
 * @return the previous error handler in use.
 */
esio_error_handler_t *
esio_set_error_handler_off(void);

/**
 * Sets the current stream handler for the process.
 * Used by the default error handling behavior, and possibly by
 * other custom error handling routines.
 *
 * @param new_handler New stream handler to use.
 *
 * @return the previous stream handler in use.
 */
esio_stream_handler_t *
esio_set_stream_handler(esio_stream_handler_t * new_handler);

/**
 * Set the default stream for error message display.  Default
 * behavior is to use stderr.
 *
 * @param new_stream New stream to use.
 *
 * @return the previous stream in use.
 */
FILE * esio_set_stream(FILE * new_stream);

/**
 * Invokes esio_error and returns the value \c esio_errno.
 * Automatically provides file and line information.
 *
 * @param reason Message to report.
 * @param esio_errno Error status to report and returned from the current
 *      function.
 */
#define ESIO_ERROR(reason, esio_errno) \
       do { \
       esio_error (reason, __FILE__, __LINE__, esio_errno) ; \
       return esio_errno ; \
       } while (0)

/**
 * Invokes esio_error using \c esio_errno and returns the value \c
 * value.  Automatically provides file and line information.
 *
 * @param reason Message to report.
 * @param esio_errno Error status to report.
 * @param value Value to return from the current function.
 */
#define ESIO_ERROR_VAL(reason, esio_errno, value) \
       do { \
       esio_error (reason, __FILE__, __LINE__, esio_errno) ; \
       return value ; \
       } while (0)

/**
 * Invokes esio_error using \c esio_errno and returns from the current
 * function.  Automatically provides file and line information.
 *
 * @param reason Message to report.
 * @param esio_errno Error status to report.
 */
#define ESIO_ERROR_VOID(reason, esio_errno) \
       do { \
       esio_error (reason, __FILE__, __LINE__, esio_errno) ; \
       return ; \
       } while (0)

/**
 * Invokes esio_error using \c esio_errno and returns NULL
 * from the current function.  Useful for out-of-memory conditions.
 * Automatically provides file and line information.
 *
 * @param reason Message to report.
 * @param esio_errno Error status to report.
 */
#define ESIO_ERROR_NULL(reason, esio_errno) \
        ESIO_ERROR_VAL(reason, esio_errno, 0)

/**
 * Invokes esio_error using \c esio_errno but \em does \em not return
 * from the current function.  Automatically provides file and line
 * information.
 *
 * @param reason Message to report.
 * @param esio_errno Error status to report.
 */
#define ESIO_ERROR_REPORT(reason, esio_errno) \
       do { \
       esio_error (reason, __FILE__, __LINE__, esio_errno) ; \
       } while (0)

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/* Internal helper macro for implementing ESIO_MPICHKx macros */
#define ESIO_MPICHKx_TEMPLATE(esio_error_macro,stmt) \
    do { \
        const int _chk_stat = (stmt); \
        if (_chk_stat != MPI_SUCCESS) { \
            char _chk_reason[255]; \
            char *_chk_mpistring = NULL; \
            int _chk_len; \
            const int _chk_string_stat \
                = MPI_Error_string(_chk_stat,_chk_mpistring,&_chk_len); \
            snprintf(_chk_reason, sizeof(_chk_reason)/sizeof(_chk_reason[0]), \
                    "Encountered MPI error code %d: %s", _chk_stat, \
                    (_chk_string_stat == MPI_SUCCESS) \
                    ? _chk_mpistring : "UNKNOWN"); \
            esio_error_macro(_chk_reason, ESIO_EFAILED); \
        } \
    } while(0)
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * Executes \c stmt once handling any resulting MPI error per \c
 * ESIO_ERROR.  Any relevant message is looked up using \c MPI_Error_string
 * and reported.  \c ESIO_EFAILED is the return value provided to \c
 * ESIO_ERROR.
 *
 * @param stmt Statement, presumably an MPI call, to be executed.
 * @note <tt>esio/mpi.h</tt> must be <tt>include</tt>d for the macro
 *       expansion to compile correctly.
 * @warning This macro expands to a not insignificant amount of code.
 *          It should not be used in performance critical regions.
 * @see PETSc's CHKERRQ for the original inspiration for this macro.
 */
#define ESIO_MPICHKQ(stmt) \
    ESIO_MPICHKx_TEMPLATE(ESIO_ERROR,stmt)

/**
 * Executes \c stmt once handling any resulting MPI error per \c
 * ESIO_ERROR_NULL.  Any relevant message is looked up using \c
 * MPI_Error_string and reported.
 *
 * @param stmt Statement, presumably an MPI call, to be executed.
 * @note <tt>esio/mpi.h</tt> must be <tt>include</tt>d for the macro
 *       expansion to compile correctly.
 * @warning This macro expands to a not insignificant amount of code.
 *          It should not be used in performance critical regions.
 * @see PETSc's CHKERRQ for the original inspiration for this macro.
 */
#define ESIO_MPICHKN(stmt) \
    ESIO_MPICHKx_TEMPLATE(ESIO_ERROR_NULL,stmt)

/**
 * Executes \c stmt once handling any resulting MPI error per \c
 * ESIO_ERROR_VOID.  Any relevant message is looked up using \c
 * MPI_Error_string and reported.
 *
 * @param stmt Statement, presumably an MPI call, to be executed.
 * @note <tt>esio/mpi.h</tt> must be <tt>include</tt>d for the macro
 *       expansion to compile correctly.
 * @warning This macro expands to a not insignificant amount of code.
 *          It should not be used in performance critical regions.
 * @see PETSc's CHKERRV for the original inspiration for this macro.
 */
#define ESIO_MPICHKV(stmt) \
    ESIO_MPICHKx_TEMPLATE(ESIO_ERROR_VOID,stmt)

/**
 * Executes \c stmt once handling any resulting MPI error per \c
 * ESIO_ERROR_REPORT. The current function \em continues \em executing. Any
 * relevant message is looked up using \c MPI_Error_string and reported.
 *
 * @param stmt Statement, presumably an MPI call, to be executed.
 * @note <tt>esio/mpi.h</tt> must be <tt>include</tt>d for the macro
 *       expansion to compile correctly.
 * @warning This macro expands to a not insignificant amount of code.
 *          It should not be used in performance critical regions.
 * @see PETSc's CHKERRV for the original inspiration for this macro.
 */
#define ESIO_MPICHKR(stmt) \
    ESIO_MPICHKx_TEMPLATE(ESIO_ERROR_REPORT,stmt)

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __ESIO_ERROR_H__ */
