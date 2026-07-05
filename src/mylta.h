#ifndef _LTA_INCLUDED
#define _LTA_INCLUDED
/* 
$Id: lta.h,v 1.7 2007/03/21 12:29:35 das Exp $" 
*/
#include "newcorr.h"
#include <stdio.h>

#define  MAX_ARG        10
#define  MAX_LINE_LEN   2048

#define EndHdr         "END_OF_HEADER"
#define StartAntenna   "*{ Antenna.def"
#define EndAntenna     "*} Antenna"
#define StartSampler   "*{ Sampler.def"
#define EndSampler     "*} Sampler"
#define StartBaseline  "*{ Baseline.def"
#define EndBaseline    "*} Baseline"
#define StartCorrsel   "*{ Corrsel.def"
#define EndCorrsel     "*} Corrsel"
#define StartCorrsys   "*{ Corrsys.def"
#define EndCorrsys     "*} Corrsys"

enum{MAX_KEYS=128};/* Maximum no of header items prased at one time */
enum{VALOFF=10};   /* Keyvalue offset from start of line */
enum{VALLEN=80};   /* Number of bytes reserved for each value */
typedef struct key_struct
{ char key[11];
  char type;
  void *value;
  int  found;
} KeyType;

typedef struct baseline_type
{ int  ant[2],band[2],samp[2];
  char antname[2][VALLEN],bandname[2][VALLEN]; 
} BaselineType;
typedef struct corr_lite
{ char   version[VALLEN],old_version[VALLEN];
  double clock;
  int    sta,statime;
  double t_unit;
  int    ants,f_enable,x_enable,chan_max,pols;
  char   clock_unit[VALLEN],sta_unit[VALLEN],statime_unit[VALLEN];
  char   t_unit_unit[VALLEN];
} CorrLite;
typedef struct samp_info
{ char ant[VALLEN],band[VALLEN];
  int  fft_id;
} SampInfo;
typedef struct vis_info
{ int antennas, samplers, baselines, channels; 
  CorrLite corrlite;
  AntennaParType antenna[MAX_ANTS];
  char  band_names[MAX_BANDS][VALLEN];
  SampInfo samp[MAX_SAMPS];
  BaselineType base[MAX_BASE];
  int   mode,lta;
  short chan_num[MAX_CHANS],samp_num[MAX_SAMPS],samp_rev_num[MAX_SAMPS];
  short s_base_num[MAX_SAMPS],x_base_num[MAX_SAMPS][MAX_SAMPS];
} VisInfo;
typedef struct lta_hdr
{ int  recl,recs;
  char rec_form[VALLEN], obs_mode[VALLEN], version[VALLEN],old_version[VALLEN];
  char d_type[VALLEN],d_size[VALLEN],d_align[VALLEN],int_rep[VALLEN];
  char fl_rep[VALLEN],byte_seq[VALLEN];
  int  flg_off,flg_rec_off,flg_ant_off,flg_smp_off;
  int  flg_bas_off,flg_dat_off;
  int  flg_size,flg_rec_size,flg_ant_size,flg_smp_size;
  int  flg_bas_size,flg_dat_size;
  int  time_off,wt_off,par_off,data_off;
  int  time_size,wt_size,par_size,data_size;
  char data_fmt[VALLEN];
  VisInfo vinfo;
} LtaHdr;

typedef struct scan_hdr
{ int      recs;
  unsigned int antmask,bandmask;
  int    seq;
  char   project[81],code[81],object[81];
  double ra_date,dec_date,mjd_src,dra,ddec;
  double f_step,rf[2],first_lo[2],bb_lo[2];
  int    net_sign[MAX_BANDS];
  float  integ;
  double mjd_ref;
  char   date_obs[81];
  int    year,day,month;
  double ist;
} ScanHdr;
typedef struct scan_hdr_tab
{ ScanHdr      shdr;
  ScanInfoType scan;
  unsigned int start_rec,recs;
} ScanHdrTab;
typedef struct lta_info
{ int recl;
  int lrecs,srecs,lbrecs,sbrecs;
  LtaHdr      lhdr;
  CorrType    corr;
  ScanHdrTab  *stab;
  int         scans;
  int         local_byte_order,data_byte_order;
  int         keep_version;
} LtaInfo;
  
int  cmp_global_hdr(LtaInfo *linfo,LtaInfo *linfo1,char *chkstr);
void convradec(double rad, char *ra, int code);
int  corrgen(CorrType *corr, char *p) ;
int  corr_ref_to_corr2(LtaInfo *linfo,CorrType *corr1); //gsb to old ghb corr
int  corr_ref_to_corr3(LtaInfo *linfo,CorrType *corr1); //gsb to new ghb corr
int  corr_ref_to_corr4(LtaInfo *linfo,CorrType *corr1); //gsb to new ghb corr
int  corr2_to_corr_ref(LtaInfo *linfo); // old ghb corr to gsb
int  corr3_to_corr_ref(LtaInfo *linfo); // new ghb corr to gsb
int  corr4_to_corr_ref(LtaInfo *linfo); // new ghb corr to gsb
int  do_average(char *ob, char *ib,LtaHdr *lhdr);
int  do_avspc(char *dbuf, LtaHdr *lhdr);
int  do_bcal(char *dbuf, char *calbuf, LtaHdr *lhdr);
int  do_median(char *obuf, char *ibuf,int recs, LtaHdr *lhdr);
int  edit_ascii_hdr(char *buf, char *keyword, void *value, char type,char *fmt);
void flipCorr(CorrType *corr) ;
void flipData(char *dbuf, LtaHdr *lh);
int  flipScanInfo(ScanInfoType *);
int  get_lta_hdr(FILE *fp,LtaInfo *linfo);
int  get_self(CorrType *corr, LtaHdr *Lhdr, char *buf, float *self,
	      double fr);
int  get_self_chan(CorrType *corr, LtaHdr *Lhdr, char *buf,float *self,
	     int ichan);
int  get_scan_hdr(FILE *fp, LtaInfo *linfo,int nscan);
long locate(char *needle, int nl, char *haystack, int hl);
int  ltaseek(FILE *fp, int nrec, int recl);
int  make_scantab(FILE *fp, LtaInfo *linfo);
char *mjd2ist_date( double mjd );
int  mk_median(int change_endian,int scan, char *mbuf,LtaInfo *linfo, FILE *fp);
int  parse_ascii_hdr(int nkey,KeyType *keys,char *buf, const char *start, const char *end);
float nr_select(unsigned long k, unsigned long n, float *x);
int  print_ascii_hdr(int nkey, KeyType *keys);
int  print_lta_hdr(LtaHdr *lhdr);
int  print_scan(ScanInfoType *scan);
int  print_scan_hdr(ScanHdr *shdr);
int scan2_to_scan1(ScanInfoType *scan, ScanInfoType *scan1);
int  sort(unsigned long n, float arr[]);
void swap_bytes(unsigned short *p, int n);
void swap_short(unsigned *p, int n);
void swap_long(void *p, int n);
void swap_d(double *p, int n);
void update_hdr(LtaInfo *linfo);
void update_signature(char *dbuf,int scanno, int recno);
void write_hdr(LtaInfo *linfo,FILE *fp);
void write_scan(ScanInfoType *scan1,int lta2, LtaInfo *linfo1,int scan_num,
		char *date_obs,FILE *fp);

#endif


