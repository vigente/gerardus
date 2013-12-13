/*  $RCSfile: qsopt.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef __QS_QSOPT_H
#define __QS_QSOPT_H

#include <stdio.h>

#ifdef WIN32 

#ifdef QSLIB_EXPORTS
#define QSLIB_INTERFACE __declspec(dllexport)
#else
//#define QSLIB_INTERFACE __declspec(dllimport)
#define QSLIB_INTERFACE 
#endif

#else 
#define QSLIB_INTERFACE extern 
#endif 

#ifdef WIN32 
typedef struct QSLIB_INTERFACE qsdata *QSprob;
typedef struct QSLIB_INTERFACE qsbasis *QSbas;
#else 
typedef struct qsdata *QSprob;
typedef struct qsbasis *QSbas;
#endif 

/****************************************************************************/
/*                                                                          */
/*                 PARAMETERS TO SPECIFY OBJECTIVE SENSE                    */
/*                                                                          */
/****************************************************************************/

#define QS_MIN       (1)
#define QS_MAX       (-1)

/****************************************************************************/
/*                                                                          */
/*                 PARAMETERS THAT CAN BE SET BY setparam                   */
/*                                                                          */
/****************************************************************************/


#define QS_PARAM_PRIMAL_PRICING    0
#define QS_PARAM_DUAL_PRICING      2
#define QS_PARAM_SIMPLEX_DISPLAY   4
#define QS_PARAM_SIMPLEX_MAX_ITERATIONS 5
#define QS_PARAM_SIMPLEX_MAX_TIME  6
#define QS_PARAM_SIMPLEX_SCALING   7


/****************************************************************************/
/*                                                                          */
/*                     VALUES FOR PRICING PARAMETERS                        */
/*                                                                          */
/****************************************************************************/

#define QS_PRICE_PDANTZIG 1
#define QS_PRICE_PDEVEX 2
#define QS_PRICE_PSTEEP 3
#define QS_PRICE_PMULTPARTIAL 4

#define QS_PRICE_DDANTZIG 6
#define QS_PRICE_DSTEEP 7
#define QS_PRICE_DMULTPARTIAL 8
#define QS_PRICE_DDEVEX 9


/****************************************************************************/
/*                                                                          */
/*                         VALUES FOR BASIS STATUS                          */
/*                                                                          */
/****************************************************************************/


#define QS_COL_BSTAT_LOWER     '0'
#define QS_COL_BSTAT_BASIC     '1'
#define QS_COL_BSTAT_UPPER     '2'
#define QS_COL_BSTAT_FREE      '3'

#define QS_ROW_BSTAT_LOWER     '0'
#define QS_ROW_BSTAT_BASIC     '1'
#define QS_ROW_BSTAT_UPPER     '2'


/****************************************************************************/
/*                                                                          */
/*         Return Status for QSopt_primal, QSopt_dual, QSget_status         */
/*                                                                          */
/****************************************************************************/

#define QS_LP_OPTIMAL           1
#define QS_LP_INFEASIBLE        2
#define QS_LP_UNBOUNDED         3
#define QS_LP_ITER_LIMIT        4
#define QS_LP_TIME_LIMIT        5
#define QS_LP_UNSOLVED          6
#define QS_LP_ABORTED		7
#define QS_LP_MODIFIED        100

/*
#define QS_LP_PRIMAL_FEASIBLE   11
#define QS_LP_PRIMAL_INFEASIBLE 12
#define QS_LP_PRIMAL_UNBOUNDED  13
#define QS_LP_DUAL_FEASIBLE     14
#define QS_LP_DUAL_INFEASIBLE   15
#define QS_LP_DUAL_UNBOUNDED    16
*/


/****************************************************************************/
/*                                                                          */
/*                         QSopt Constants                                  */
/*                                                                          */
/****************************************************************************/

#define QS_MAXDOUBLE (1e30)


/****************************************************************************/
/*                                                                          */
/*                      QSopt Library Functions                             */
/*                                                                          */
/****************************************************************************/
#ifdef  __cplusplus 
extern "C" {
#endif

#ifdef WIN32
/* 
 *  in WINDOWS we make 
 *     solver_main/reader_main part of DLL
 */
QSLIB_INTERFACE int solver_main(int argc, char **argv); 
QSLIB_INTERFACE int reader_main(int argc, char **argv); 
#endif 

QSLIB_INTERFACE void
    QSfree (void *ptr),
    QSfree_prob (QSprob p),
    QSfree_basis (QSbas B);

QSLIB_INTERFACE int
    QSopt_primal (QSprob p, int *status),
    QSopt_dual (QSprob p, int *status),
    QSopt_pivotin_col (QSprob p, int ccnt, int *clist),
    QSopt_pivotin_row (QSprob p, int rcnt, int *rlist),
    QSopt_strongbranch (QSprob p, int ncand, int *candidatelist,
        double *xlist, double *down_vals, double *up_vals, int iterations,
        double objbound),
    QSchange_objsense (QSprob p, int newsense),
    QSget_objsense (QSprob p, int *newsense),
    QSnew_col (QSprob p, double obj, double lower, double upper,
        const char *name),
    QSadd_cols (QSprob p, int num, int *cmatcnt, int *cmatbeg, int *cmatind,
        double *cmatval, double *obj, double *lower, double *upper,
        const char **names),
    QSadd_col (QSprob p, int cnt, int *cmatind, double *cmatval, double obj,
        double lower, double upper, const char *name),
    QSnew_row (QSprob p, double rhs, char sense, const char *name),
    QSadd_rows (QSprob p, int num, int *rmatcnt, int *rmatbeg,
        int *rmatind, double *rmatval, double *rhs, char *sense,
        const char **names),
    QSadd_row (QSprob p, int cnt, int *rmatind, double *rmatval, double rhs,
        char sense, const char *name),
    QSdelete_rows (QSprob p, int num, int *dellist),
    QSdelete_row (QSprob p, int rowindex),
    QSdelete_setrows (QSprob p, int *flags),
    QSdelete_named_row (QSprob p, const char *rowname),
    QSdelete_named_rows_list (QSprob p, int num, const char **rownames),
    QSdelete_cols (QSprob p, int num, int *dellist),
    QSdelete_col (QSprob p, int colindex),
    QSdelete_setcols (QSprob p, int *flags),
    QSdelete_named_column (QSprob p, const char *colname),
    QSdelete_named_columns_list (QSprob p, int num, const char **colnames),
    QSchange_senses (QSprob p, int num, int *rowlist, char *sense),
    QSchange_sense (QSprob p, int rowindex, char sense),
    QSchange_coef (QSprob p, int rowindex, int colindex, double coef),
    QSchange_objcoef (QSprob p, int indx, double coef),
    QSchange_rhscoef (QSprob p, int indx, double coef),
    QSchange_bounds (QSprob p, int num, int *collist, char *lu,
        double *bounds),
    QSchange_bound (QSprob p, int indx, char lu, double bound),
    QSload_basis (QSprob p, QSbas B),
    QSread_and_load_basis (QSprob p, const char *filename),
    QSload_basis_array (QSprob p, char *cstat, char *rstat),
    QSload_basis_and_row_norms_array (QSprob p, char *cstat, char *rstat,
        double *rownorms),
    QSget_basis_array (QSprob p, char *cstat, char *rstat),
    QSget_basis_and_row_norms_array (QSprob p, char *cstat, char *rstat,
        double *rownorms),
    QSget_binv_row (QSprob p, int indx, double *binvrow),
    QSget_tableau_row (QSprob p, int indx, double *tableaurow),
    QSget_basis_order (QSprob p, int *basorder),
    QSget_status (QSprob p, int *status),
    QSget_solution (QSprob p, double *value, double *x, double *pi,
        double *slack, double *rc),
    QSget_objval (QSprob p, double *value),
    QSget_pi_array (QSprob p, double *pi),
    QSget_rc_array (QSprob p, double *rc),
    QSget_x_array (QSprob p, double *x),
    QSget_slack_array (QSprob p, double *slack),
    QSget_infeas_array (QSprob p, double *pi),
    QSget_colcount (QSprob p),
    QSget_rowcount (QSprob p),
    QSget_nzcount (QSprob p),
    QSget_obj (QSprob p, double *obj),
    QSget_rhs (QSprob p, double *rhs),
    QSget_rows_list (QSprob p, int num, int *rowlist, int **rowcnt,
        int **rowbeg, int **rowind, double **rowval, double **rhs,
        char **sense, char ***names),
    QSget_rows (QSprob p, int **rowcnt, int **rowbeg, int **rowind,
        double **rowval, double **rhs, char **sense, char ***names),
    QSget_columns_list (QSprob p, int num, int *collist, int **colcnt,
        int **colbeg, int **colind, double **colval, double **obj,
        double **lower, double **upper, char ***names),
    QSget_columns (QSprob p, int **colcnt, int **colbeg, int **colind,
        double **colval, double **obj, double **lower, double **upper,
        char ***names),
    QSget_senses (QSprob p, char *senses),
    QSget_rownames (QSprob p, char **rownames),
    QSget_colnames (QSprob p, char **colnames),
    QSget_bound (QSprob p, int colindex, char lu, double *bound),
    QSget_bounds (QSprob p, double *lower, double *upper),
    QSget_intflags (QSprob p, int *intflags),
    QSget_intcount (QSprob p, int *count),
    QSget_column_index (QSprob p, const char *name, int *colindex),
    QSget_row_index (QSprob p, const char *name, int *rowindex),
    QSget_named_x (QSprob p, const char *colname, double *val),
    QSget_named_rc (QSprob p, const char *colname, double *val),
    QSget_named_pi (QSprob p, const char *rowname, double *val),
    QSget_named_slack (QSprob p, const char *rowname, double *val),
    QScompute_row_norms (QSprob p),
    QSwrite_prob (QSprob p, const char *filename, const char *filetype),
    QSwrite_prob_file (QSprob p, FILE *file, const char *filetype),
    QSwrite_basis (QSprob p, QSbas B, const char *filename),
    QStest_row_norms (QSprob p),
    QSset_param (QSprob p, int whichparam, int newvalue),
    QSset_param_double (QSprob p, int whichparam, double newvalue),
    QSget_param (QSprob p, int whichparam, int *value),
    QSget_param_double (QSprob p, int whichparam, double *value);

QSLIB_INTERFACE char* QSget_probname (QSprob p); 
QSLIB_INTERFACE char* QSget_objname (QSprob p); 
QSLIB_INTERFACE char* QSversion (void); 

QSLIB_INTERFACE QSprob
   QScreate_prob (const char *name, int objsense),
   QSread_prob (const char *filename, const char *filetype),
   QSload_prob (const char *probname, int ncols, int nrows,
        int *cmatcnt, int *cmatbeg, int *cmatind, double *cmatval,
        int objsense, double *obj, double *rhs, char *sense,
        double *lower, double *upper, const char **colnames,
        const char **rownames),
   QScopy_prob (QSprob p, const char *newname);

QSLIB_INTERFACE QSbas
   QSget_basis (QSprob p),
   QSread_basis (QSprob p, const char *filename);

#ifdef  __cplusplus 
}
#endif

/****************************************************************************
 *
 * This is the undocumented part of the QSlib interface 
 *
 ****************************************************************************/
/* 
 * functions to facilitate line by line reading from other sources than 
 * files from within MPS/LP parsers  
 * 
 * functions to facilitate the collection of error information instead of 
 * having the parsers print messages to stderr
 *                              by mps/lp format writers
 * 
 * a problem's reporter is used by the solver code to provide important 
 * feedback/progress information
 */

#ifdef WIN32
typedef struct QSLIB_INTERFACE qsline_reader     *QSline_reader;
typedef struct QSLIB_INTERFACE qsformat_error    *QSformat_error;
typedef struct QSLIB_INTERFACE qserror_collector *QSerror_collector;
typedef struct QSLIB_INTERFACE qserror_memory    *QSerror_memory;
#else 
typedef struct qsline_reader      *QSline_reader;
typedef struct qsformat_error     *QSformat_error;
typedef struct qserror_collector  *QSerror_collector;
typedef struct qserror_memory     *QSerror_memory;
#endif 

#ifdef  __cplusplus 
extern "C" {
#endif

/****************************************************************************
 * error collection
 */
#define QS_DATA_ERROR			0
#define QS_DATA_WARN			1
#define QS_MPS_FORMAT_ERROR		2
#define QS_MPS_FORMAT_WARN		3
#define QS_LP_FORMAT_ERROR		4
#define QS_LP_FORMAT_WARN		5
#define QS_INPUT_NERROR        8

QSLIB_INTERFACE const char* QSformat_error_type_string(int tp); 

QSLIB_INTERFACE int 		QSerror_get_type(QSformat_error error); 
QSLIB_INTERFACE const char* QSerror_get_desc(QSformat_error error); 
QSLIB_INTERFACE int 		QSerror_get_line_number(QSformat_error error); 
QSLIB_INTERFACE int 		QSerror_get_pos(QSformat_error error); 
QSLIB_INTERFACE const char* QSerror_get_line(QSformat_error error); 
QSLIB_INTERFACE void        QSerror_print(FILE *f, QSformat_error error); 

QSLIB_INTERFACE 
    QSerror_collector QSerror_collector_new(void *fct, void *dest); 
QSLIB_INTERFACE 
    QSerror_collector QSerror_memory_collector_new(QSerror_memory mem); 
QSLIB_INTERFACE void 
    QSerror_collector_free(QSerror_collector c); 

/****************************************************************************
 * line reader 
 */
QSLIB_INTERFACE QSline_reader 
    QSline_reader_new(void *fct, void* data_src); 
    /* reader->read_line_fct defaults to fgets */

QSLIB_INTERFACE void
    QSline_reader_free(QSline_reader reader); 

QSLIB_INTERFACE void
    QSline_reader_set_error_collector(QSline_reader reader,  
				      QSerror_collector collector); 

QSLIB_INTERFACE char* 
    QSline_reader_get(QSline_reader reader, char *s, int size);

QSLIB_INTERFACE QSprob 
    QSget_prob(QSline_reader reader, 
	       const char* probname, const char *filetype);
    /* the MPS and LP parsers uses the fct from reader 
     * to get to next input line */


/****************************************************************************
 * error memory 
 */
QSLIB_INTERFACE QSerror_memory  QSerror_memory_create(char takeErrorLines); 
QSLIB_INTERFACE void            QSerror_memory_free(QSerror_memory mem);

QSLIB_INTERFACE int            
	QSerror_memory_get_nof(QSerror_memory mem, int error_type);
QSLIB_INTERFACE int            
	QSerror_memory_get_nerrors(QSerror_memory mem);

QSLIB_INTERFACE QSformat_error 
	QSerror_memory_get_last_error(QSerror_memory mem);
QSLIB_INTERFACE QSformat_error            
	QSerror_memory_get_prev_error(QSformat_error e);

/**************************************************************************** 
 * reporter for solver feedback 
 */
QSLIB_INTERFACE void QSset_reporter(QSprob prob, int iterskip, void *fct, void *dest); 

QSLIB_INTERFACE int QSreport_prob (QSprob p, const char *filetype, 
				   QSerror_collector c);

#ifdef  __cplusplus 
}
#endif
#endif /* __QS_QSOPT_H */
