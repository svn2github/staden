#ifndef _TEXT_OUTPUT_H_
#define _TEXT_OUTPUT_H_

/*
 * Usage: verror(priority, format, args...);
 * NB: don't pass more than 8K per call
 */
#define ERR_WARN 0
#define ERR_FATAL 1
void verror(int priority, char *name, char *fmt, ...);

/*
 * Usage: vmessage(format, args...);
 * NB: don't pass more than 8K per call
 */
void vmessage(char *fmt, ...);

/*
 * Adds a new header to the text output window.
 */
void vfuncheader(char *fmt, ...);

/*
 * As vfuncheader, but only outputting when necessary.
 */
void vfuncgroup(int group, char *fmt, ...);

void vfuncparams(char *fmt, ...);

void UpdateTextOutput(void);

void start_message(void);
void end_message(char *parent);

/*
 * Outputs data to a persistent log file. The log file is always appended to,
 * but only contains headers and error messages.
 *
 * 'fn' is the filename to log to, or NULL if we use whichever one is currently
 * open. Not initialising the function will not cause any problems.
 *
 * Giving a filename of "" will close the log file.
 */
void log_file(char *fn, char *message);

/* 
 * Controls whether vmessage output should also be written to the log file 
 * (in addition to vfuncheader and verror messages). 
 * A value of 0 means do not log. Any other values implies logging. 
 */ 
void log_vmessage(int log);

#endif
