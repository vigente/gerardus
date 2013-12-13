/* $RCSfile: reporter.h,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:57:39 $ */
#ifndef REPORTER_FILE
#define REPORTER_FILE

typedef int (*qsreport_string_fct)(void* dest, const char *s); 

typedef struct qsstring_reporter
{
    qsreport_string_fct report_fct; 
    void *dest; 
} qsstring_reporter;

extern int ILL_fprintf(void *dest, const char *s); 

void ILLstring_reporter_init(qsstring_reporter *reporter, 
                             qsreport_string_fct fct,
                             void* dest);

void ILLstring_reporter_copy(qsstring_reporter *dest,
                             qsstring_reporter *src);

#define ILLstring_report(s, reporter)  \
	((reporter)->report_fct((reporter)->dest, s) < 0)
	                 /* used from with ILL fct to report progress */

/* REPORTER_FILE */
#endif  
