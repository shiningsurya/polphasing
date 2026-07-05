#include "mylta.h"
#include "newcorr.h"
#include "protocol.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>
#include <unistd.h>

//static char rcsid[]="$Id: ltasub.c,v 1.15 2007/04/18 08:44:28 das Exp $";

#define SETKEY(n,k,t,p) {strcpy(keys[n].key,k);\
                         keys[n].type=t;\
                         keys[n].value=p;\
                         n++;\
                         if(n==MAX_KEYS)\
                         {fprintf(stderr,"%s MAX_KEYS exceeded\n",__FILE__);\
                           return 1;}}

typedef struct xtra_scan_par
{ char rf[VALLEN],first_lo[VALLEN],bb_lo[VALLEN],net_sign[VALLEN];
} XScanPar;
typedef struct xtra_corrlite_par
{ char clock[VALLEN],sta[VALLEN],statime[VALLEN],t_unit[VALLEN];
} XCorrLitePar;

#define RAD2DEG(x)  ((x)*180.0/M_PI)
#define DEG2RAD(x)  ((x)/180.0*M_PI)

char *mjd2ist_date( double mjd )
{
    double J = mjd + 2400000.5 +5.5/24.0; 
    static char date[32] ;
    int month, day;
    long year, a, c, d, x, y, jd;
    double dd;
    static const char *Month[13] ={"Jan", "Feb", "Mar", "Apr",
                  "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" } ;
    if( J < 1721425.5 ) return 0 ; /* January 1.0, 1 A.D.  dont accept BC ! */

    jd = J + 0.5; /* round Julian date up to integer */

    /* Find the number of Gregorian centuries
     * since March 1, 4801 B.C.
     */
    a = (100*jd + 3204500L)/3652425L;

    /* Transform to Julian calendar by adding in Gregorian century years
     * that are not leap years.
     * Subtract 97 days to shift origin of JD to March 1.
     * Add 122 days for magic arithmetic algorithm.
     * Add four years to ensure the first leap year is detected.
     */
    c = jd + 1486;
    if( jd >= 2299160.5 )
        c += a - a/4;
    else
        c += 38;
    /* Offset 122 days, which is where the magic arithmetic
     * month formula sequence starts (March 1 = 4 * 30.6 = 122.4).
     */
    d = (100*c - 12210L)/36525L;
    x = (36525L * d)/100L; /* Days in that many whole Julian years */

    /* Find month and day. */
    y = ((c-x)*100L)/3061L;
    day = c - x - ((306L*y)/10L);
    month = y - 1;
    if( y > 13 ) month -= 12;

    /* Get the year right. */
    year = d - 4715;
    if( month > 2 ) year -= 1;

    a = (jd + 1) % 7; /* Day of the week. */

    dd = day + J - jd + 0.5; /* Fractional part of day. */
    day = dd;
    { int h,m,s ;
      s = (dd-day)*24*3600 ;
      h = s/3600 ;
      m = (s-h*3600)/60 ;
      s -= (h*3600 + m*60) ;
      sprintf(date,"%02d/%s/%04ld %02d:%02d:%02d",day,Month[month-1],year,h,m,s);
    }
    return date ;
}
void convradec(double rad, char *ra, int code)
{
  double ras,temp;
  int    rah,ram;
  char   rasgn,rahr[3],ramin[3],rasec[7];


  
  if(code == 0){/* output hr min sec */
      while(rad < 0.0)
          rad += 360.0;
      temp = (rad - 360.0*(int)floor(rad/360.0))/15.0;                  
  }else{/* output deg min sec */
      while(fabs(rad) > 180.0)
          rad -= (fabs(rad)/rad)*360.0;
      temp = fabs(rad);                      
  }

  rah  = (int)floor(temp);
  if(rah < 10){
    sprintf(rahr,"0%1i",rah);
  }
  else{
    sprintf(rahr,"%2i",rah);
  }
  temp = (temp -rah)*60.0;
  ram  = (int)floor(temp);
  if(ram < 10){
    sprintf(ramin,"0%1i",ram);
  }
  else{
    sprintf(ramin,"%2i",ram);
  }
  ras = (temp - ram)*60.0;
  if(ras < 10.0){
    sprintf(rasec,"0%4.2f",ras);
  }
  else{
    sprintf(rasec,"%4.2f",ras);
  }

  if(code==0){
    sprintf(ra,"%2sh%2sm%5s",rahr,ramin,rasec);
    return;
  }
  rasgn = rad < 0.0 ? '-' : '+';
  sprintf(ra,"%1c%2sd%2s'%5s\"",rasgn,rahr,ramin,rasec);
  return;
}
int print_scan_hdr(ScanHdr *shdr)
{ printf("*** Begin ASCII SCAN HEADER ***\n");
  printf("HDR_RECS= %d\n",shdr->recs);
  printf("ANTMASK = %x\n",shdr->antmask);
  printf("BANDMASK= %x\n",shdr->bandmask);
  printf("SEQ     = %d\n",shdr->seq);
  printf("PROJECT = %s\n",shdr->project);
  printf("CODE    = %s\n",shdr->code);
  printf("OBJECT  = %s\n",shdr->object);
  printf("RA-DATE = %lf\n",shdr->ra_date);
  printf("DEC-DATE= %lf\n",shdr->dec_date);
  printf("MJD_SRC = %lf\n",shdr->mjd_src);
  printf("DRA/DT  = %lf\n",shdr->dra);
  printf("DDEC/DT = %lf\n",shdr->ddec);
  printf("F_STEP  = %lf\n",shdr->f_step);
  printf("RF      = %lf %lf\n",shdr->rf[0],shdr->rf[1]);
  printf("FIRST_LO= %lf %lf\n",shdr->first_lo[0],shdr->first_lo[1]);
  printf("BB_LO   = %lf %lf\n",shdr->bb_lo[0],shdr->bb_lo[1]);
  printf("NET_SIGN= %d %d %d %d\n",shdr->net_sign[0],shdr->net_sign[1],
	            shdr->net_sign[2],shdr->net_sign[3]);
  printf("DATE-OBS= %s\n",shdr->date_obs);
  printf("INTEG   = %lf\n",shdr->integ);
  printf("MJD_REF = %lf\n",shdr->mjd_ref);
  printf("*** End ASCII SCAN HEADER ***\n");
  return 0;
} 
int print_scan(ScanInfoType *scan)
{  SourceParType *src = &scan->source;
   char          ras[32],decs[32];
   printf("*** Begin BINARY SCAN  HEADER ***\n");
   printf("OBJECT   = %s\n",src->object);
   printf("FLUX     = {%f,%f,%f,%f}\n",src->flux.i,src->flux.q,src->flux.u,
	  src->flux.v);
   printf("MJD_SRC     = %12.6f\n",src->mjd0);
   convradec(RAD2DEG(src->ra_mean),ras,0);
   convradec(RAD2DEG(src->dec_mean),decs,1);
   printf("RA_MEAN  = %s\n",ras);
   printf("DEC_MEAN = %s\n",decs);
   convradec(RAD2DEG(src->ra_app),ras,0);
   convradec(RAD2DEG(src->dec_app),decs,1);
   printf("RA_APP   = %s\n",ras);
   printf("DEC_APP  = %s\n",decs);
   convradec(RAD2DEG(src->dra),ras,0);
   convradec(RAD2DEG(src->ddec),decs,1);
   printf("DRA      = %s\n",ras);
   printf("DDEC     = %s\n",decs);
   printf("FREQ     = %-13.3f %-13.3f\n",src->freq[0],src->freq[1]);
   printf("FIRST_LO = %-13.3f %-13.3f\n",src->first_lo[0],src->first_lo[1]);
   printf("BB_LO    = %-13.3f %-13.3f\n",src->bb_lo[0],src->bb_lo[1]);
   printf("REST_FREQ= %-13.3f %-13.3f\n",src->rest_freq[0],src->rest_freq[1]);
   printf("LSR_VEL  = %-13.3f %-13.3f\n",src->lsrvel[0],src->lsrvel[1]);
   printf("CH_WIDTH = %-13.3f\n",src->ch_width);
   printf("NET_SIGN = {%d,%d,%d,%d}\n",src->net_sign[0],src->net_sign[1],
	              src->net_sign[2],src->net_sign[3]);
   printf("ANTMASK  = %x\n",src->antmask);
   printf("BANDMASK = %x\n",src->bandmask);
   printf("CALCODE  = %-4d\n",src->calcode);
   printf("QUAL     = %-4d\n",src->qual);
   printf("*** End BINARY SCAN  HEADER ***\n");
   return 0;
}
int set_lta_keys(LtaHdr *lhdr, KeyType *keys)
{ int nkey;
  nkey=0; 
  /*{ Init.def */
  SETKEY(nkey,"RECL    ",'i',&(lhdr->recl));
  SETKEY(nkey,"HDR_RECS",'i',&(lhdr->recs));
  SETKEY(nkey,"REC_FORM",'s', (lhdr->rec_form));
  SETKEY(nkey,"OBS_MODE",'s', (lhdr->obs_mode));
  SETKEY(nkey,"VERSION ",'S', (lhdr->version));
  SETKEY(nkey,"D_TYPE  ",'S', (lhdr->d_type));
  SETKEY(nkey,"D_SIZE  ",'S', (lhdr->d_size));
  SETKEY(nkey,"D_ALIGN ",'S', (lhdr->d_align));
  SETKEY(nkey,"INT_REP ",'S', (lhdr->int_rep));
  SETKEY(nkey,"FL_REP  ",'s', (lhdr->fl_rep));
  SETKEY(nkey,"BYTE_SEQ",'S', (lhdr->byte_seq));
  /*} Init.def */
  /*{ AddCorr.def */
  SETKEY(nkey,"ANTENNAS",'i',&(lhdr->vinfo.antennas));
  SETKEY(nkey,"SAMPLERS",'i',&(lhdr->vinfo.samplers));
  SETKEY(nkey,"BASELINE",'i',&(lhdr->vinfo.baselines));
  SETKEY(nkey,"CHANNELS",'i',&(lhdr->vinfo.channels));
  SETKEY(nkey,"FLG_OFF ",'i',&(lhdr->flg_off));
  SETKEY(nkey,"FLG_SIZE",'i',&(lhdr->flg_size));
  SETKEY(nkey,"FLGRECOF",'i',&(lhdr->flg_rec_off));
  SETKEY(nkey,"FLGRECSZ",'i',&(lhdr->flg_rec_size));
  SETKEY(nkey,"FLGANTOF",'i',&(lhdr->flg_ant_off));
  SETKEY(nkey,"FLGANTSZ",'i',&(lhdr->flg_ant_size));
  SETKEY(nkey,"FLGSMPOF",'i',&(lhdr->flg_smp_off));
  SETKEY(nkey,"FLGSMPSZ",'i',&(lhdr->flg_smp_size));
  SETKEY(nkey,"FLGBASOF",'i',&(lhdr->flg_bas_off));
  SETKEY(nkey,"FLGBASSZ",'i',&(lhdr->flg_bas_size));
  SETKEY(nkey,"FLGDATOF",'i',&(lhdr->flg_dat_off));
  SETKEY(nkey,"FLGDATSZ",'i',&(lhdr->flg_dat_size));
  SETKEY(nkey,"TIME_OFF",'i',&(lhdr->time_off));
  SETKEY(nkey,"TIMESIZE",'i',&(lhdr->time_size));
  SETKEY(nkey,"WT_OFF  ",'i',&(lhdr->wt_off));
  SETKEY(nkey,"WT_SIZE ",'i',&(lhdr->wt_size));
  SETKEY(nkey,"PAR_OFF ",'i',&(lhdr->par_off));
  SETKEY(nkey,"PAR_SIZE",'i',&(lhdr->par_size));
  SETKEY(nkey,"DATA_OFF",'i',&(lhdr->data_off));
  SETKEY(nkey,"DATASIZE",'i',&(lhdr->data_size));
  SETKEY(nkey,"DATAFMT ",'s', (lhdr->data_fmt));
  /*} AddCorr.def */
  return nkey;
}
int set_scan_keys(ScanHdr *shdr, KeyType *keys, XScanPar *tmp)
{ int nkey;
  nkey=0; 
  /*{ SubArray0 */
  SETKEY(nkey,"HDR_RECS",'i',&(shdr->recs));
  SETKEY(nkey,"ANTMASK ",'x',&(shdr->antmask));
  SETKEY(nkey,"BANDMASK",'x',&(shdr->bandmask));
  SETKEY(nkey,"SEQ     ",'i',&(shdr->seq));
  SETKEY(nkey,"PROJECT ",'s',&(shdr->project));
  SETKEY(nkey,"CODE    ",'s',&(shdr->code));
  SETKEY(nkey,"OBJECT  ",'s',&(shdr->object));
  SETKEY(nkey,"RA-DATE ",'d',&(shdr->ra_date));
  SETKEY(nkey,"DEC-DATE",'d',&(shdr->dec_date));
  SETKEY(nkey,"MJD_SRC ",'d',&(shdr->mjd_src));
  SETKEY(nkey,"DRA/DT  ",'d',&(shdr->dra));
  SETKEY(nkey,"DDEC/DT ",'d',&(shdr->ddec));
  SETKEY(nkey,"F_STEP  ",'d',&(shdr->f_step));
  SETKEY(nkey,"RF      ",'S',&(tmp->rf));
  SETKEY(nkey,"FIRST_LO",'S',&(tmp->first_lo));
  SETKEY(nkey,"BB_LO   ",'S',&(tmp->bb_lo));
  SETKEY(nkey,"NET_SIGN",'S',&(tmp->net_sign));
  /*} SubArray0 */
  /*{ ExtraScan */
  SETKEY(nkey,"DATE-OBS",'S',&(shdr->date_obs));
  SETKEY(nkey,"INTEG   ",'f',&(shdr->integ));
  SETKEY(nkey,"MJD_REF ",'d',&(shdr->mjd_ref));
  /*} ExtraScan */
  return nkey;
} 
int parse_ascii_hdr(int nkey, KeyType *keys, char *buf,const char *start, const char *end)
{ char *p,*q,line[VALLEN];
  int  ngot=0,i,j,found=0;

  for(i=0;i<nkey;i++)keys[i].found=0;
  p=buf;
  if(start!=NULL)
  { while((p=strchr(p,'\n')) != NULL)
    { p++; if(!strncmp(p,EndHdr,strlen(EndHdr))) break;
      if(strncmp(p,start,strlen(start))) continue;
      else {found=1; break;}
    }
    if(!found){fprintf(stderr,"No %s section!\n",start); return -1;}
  }

  while((p=strchr(p,'\n')) != NULL)
  { p++; if(!strncmp(p,EndHdr,strlen(EndHdr))) break;
    if(end!=NULL && !strncmp(p,end,strlen(end))) break;
    for(i=0;i<nkey;i++)
    { if(strchr(p,'=')==NULL)continue;
      if(keys[i].found) continue; 
      if(!strncmp(p,keys[i].key,strlen(keys[i].key)))
      { keys[i].found=1;ngot++;
        strncpy(line,p+VALOFF,VALLEN-VALOFF-1);line[VALLEN-VALOFF-1]='\0';
        switch(keys[i].type)
        { case 'c':sscanf(line,"%s",(char *)keys[i].value);break;
	  case 'f':sscanf(line,"%f",(float *)keys[i].value);break;
          case 'i':sscanf(line,"%i",(int *)keys[i].value);break;
          case 'd':sscanf(line,"%lf",(double *)keys[i].value);break;
	  case 's':sscanf(line,"%s",(char *)keys[i].value);break;
  	  case 'S':for(j=VALLEN-VALOFF-2;j>-1;j--) /* del trailing space*/
	           if(isspace(line[j]))line[j]='\0';
		   else break;
	           q=line;                        /* del leading space*/
		   while(isspace(*q))q++; strcpy((char *)keys[i].value,q);
	           break;
	  case 'x':sscanf(line,"%x",(unsigned int *)keys[i].value);break;
	  default :fprintf(stderr,"Illegal Type %c\n",keys[i].type); return -1;		
	}
	break;
      }
    }
  }
  return ngot;
}
int print_ascii_hdr(int nkey, KeyType *keys)
{ int i;

  printf("** Start of ASCII Header **\n"); 
  for(i=0;i<nkey;i++)
  { switch(keys[i].type)
    { case 'c':printf("%-8.8s %1c %s\n",keys[i].key,keys[i].type,
		      (char *)keys[i].value); break;
      case 'i':printf("%-8.8s %1c %d\n",keys[i].key,keys[i].type,
		      *(int *)keys[i].value); break;
      case 'f':printf("%-8.8s %1c %f\n",keys[i].key,keys[i].type,
		      *(float *)keys[i].value); break;
      case 'd':printf("%-8.8s %1c %f\n",keys[i].key,keys[i].type,
		      *(double *)keys[i].value); break;
      case 's':printf("%-8.8s %1c %s\n",keys[i].key,keys[i].type,
		      (char *)keys[i].value); break;
      case 'S':printf("%-8.8s %1c %s\n",keys[i].key,keys[i].type,
		      (char *)keys[i].value); break;
       default:printf("Illegal Keytype %c\n",keys[i].type);
    }
  }
  printf("** End of ASCII Header **\n"); 
  return 0;
}
void fill(char *line)
{ int i,j;
  j=strlen(line);
  for(i=j;i<VALLEN;i++)line[i]=' ';
  line[VALLEN-1]='\n';
}
int edit_ascii_hdr(char *buf, char *keyword, void *value, char type,char *fmt)
{ char line[VALLEN],*p;
  int  j,found;

  j=VALLEN-VALOFF-1;

  p=buf; found=0;
  while((p=strchr(p,'\n')) != NULL)
  { p++; if(!strncmp(p,"END_OF_HEADER",14)) break;
    { if(strchr(p,'=')==NULL)continue;
      if(!strncmp(p,keyword,strlen(keyword)))
      { found=1;
        switch(type)
        { case 'c': sprintf(line,fmt,*(char *)value);
 	            fill(line);
      	            memcpy(p+VALOFF,line,j);
		    break;
	  case 'f': sprintf(line,fmt,*(float *)value);
                    fill(line);
      	            memcpy(p+VALOFF,line,j);
		    break;
          case 'i': sprintf(line,fmt,*(int *)value);
                    fill(line);
      	            memcpy(p+VALOFF,line,j);
		    break;
          case 'd': sprintf(line,fmt,*(double *)value);
                    fill(line);
      	            memcpy(p+VALOFF,line,j);
		    break;
   	  case 'S':
	  case 's': sprintf(line,fmt,(char *)value);
                    fill(line);
      	            memcpy(p+VALOFF,line,j);
		    break;
  	  default : fprintf(stderr,"Illegal Type %c\n",type);
	            return -1;	
	}
	break;
      }
    }
  }
  if(found) return 0;
  else return -1;

}
int get_corr_lite(char *buf, CorrLite *corrlite)
{ int      nkey; 
  KeyType  keys[MAX_KEYS];
  XCorrLitePar tmp;
  nkey=0; 
  SETKEY(nkey,"VERSION ",'S', (corrlite->version));
  SETKEY(nkey,"CLOCK   ",'S',&(tmp.clock));
  SETKEY(nkey,"STA     ",'S',&(tmp.sta));
  SETKEY(nkey,"STATIME ",'S',&(tmp.statime));
  SETKEY(nkey,"T_UNIT  ",'S',&(tmp.t_unit));
  SETKEY(nkey,"ANTS    ",'i',&(corrlite->ants));
  SETKEY(nkey,"F_ENABLE",'i',&(corrlite->f_enable));
  SETKEY(nkey,"X_ENABLE",'i',&(corrlite->x_enable));
  SETKEY(nkey,"CHAN_MAX",'i',&(corrlite->chan_max));
  SETKEY(nkey,"POLS    ",'i',&(corrlite->pols));
  if(parse_ascii_hdr(nkey,keys,buf,StartCorrsys,EndCorrsys) !=nkey)
  { fprintf(stderr,"Too few keywords found\n"); return -1; }

  sscanf(tmp.clock,  "%lf %s\n",&corrlite->clock,corrlite->clock_unit);
  sscanf(tmp.sta,    "%d %s\n", &corrlite->sta,corrlite->sta_unit);
  sscanf(tmp.statime,"%d %s\n", &corrlite->statime,corrlite->statime_unit);
  sscanf(tmp.t_unit, "%lf %s\n",&corrlite->t_unit,corrlite->t_unit_unit);
  strcpy(corrlite->old_version,corrlite->version);
  return 0;
}
int get_corrsel(char *buf, VisInfo *vinfo)
{ char *p=buf,word[VALLEN];
  int   found=0,i,j;
  int get_range(char *str, short *arr, int max);

  while((p=strchr(p,'\n')) != NULL)
  { p++; if(!strncmp(p,EndHdr,strlen(EndHdr))) break;
    if(strncmp(p,StartCorrsel,strlen(StartCorrsel))) continue;
    else {found=1; break;}
  }
  if(!found){fprintf(stderr,"No %s section!\n",StartCorrsel); return -1;}

  while((p=strchr(p,'\n')) != NULL)
  { p++; if(!strncmp(p,EndHdr,strlen(EndHdr))) break;
    if(!strncmp(p,EndCorrsel,strlen(EndCorrsel))) break;
    if(p[0]=='*') continue;
    if(!strncmp(p,"MODE",4)){sscanf(p+VALOFF,"%d",&vinfo->mode); continue;}
    if(!strncmp(p,"LTA",3)){sscanf(p+VALOFF,"%d",&vinfo->lta); continue;}
    /*  if(!strncmp(p,"SAMP_NUM",8))   SAMP_NUM keyword depreciated!!
        { for(i=0;i<VALLEN;i++){word[i]=*(p+VALOFF+i);if(word[i]=='\n')break;}
        if(i<VALLEN)word[i]='\0';else word[VALLEN-1]='\0';
        if(get_range(word,vinfo->samp_num,MAX_SAMPS)==0)
        { fprintf(stderr,"Bad Sampler Range %s\n",p+VALOFF); return -1;}
        continue;
	}
    */
    if(!strncmp(p,"CHAN_NUM",8))
    { for(i=0;i<VALLEN;i++){word[i]=*(p+VALOFF+i);if(word[i]=='\n')break;}
      if(i<VALLEN)word[i]='\0';else word[VALLEN-1]='\0';
      if(get_range(word,vinfo->chan_num,MAX_CHANS)==0)
      { fprintf(stderr,"Bad Channel Range %s\n",p+VALOFF); return -1;}
      continue;
    }
  }
  for(i=0;i<MAX_SAMPS;i++)
  { vinfo->s_base_num[i]=MAX_SAMPS;
    for(j=0;j<MAX_SAMPS;j++)vinfo->x_base_num[i][j]=MAX_BASE;
  }
  for(i=0;i<MAX_SAMPS;i++)
  { vinfo->samp_rev_num[i]=MAX_SAMPS;
    for(j=0;j<vinfo->samplers;j++)
      if(i==vinfo->samp_num[j])vinfo->samp_rev_num[i]=j; 
  }
  
  for(i=0;i<vinfo->baselines;i++)
  { BaselineType *base=&vinfo->base[i];
    short *rev_num=vinfo->samp_rev_num;
    int s0=rev_num[base->samp[0]],s1=rev_num[base->samp[1]];
    if(s0==s1)vinfo->s_base_num[s0]=i;
    else
    {vinfo->x_base_num[s0][s1]=i;vinfo->x_base_num[s1][s0]=-i;}
  }
  return 0;
}
int get_ant_info(char *buf,AntennaParType *antenna)
{ char *p=buf;
  AntennaParType *ant=antenna; 
  int   i,k,found=0;
  
  while((p=strchr(p,'\n')) != NULL)
  { p++; if(!strncmp(p,EndHdr,strlen(EndHdr))) break;
    if(strncmp(p,StartAntenna,strlen(StartAntenna))) continue;
    else {found=1; break;}
  }
  if(!found){fprintf(stderr,"No %s section!\n",StartAntenna); return -1;}

  k=0;
  while((p=strchr(p,'\n')) != NULL)
  { p++; if(!strncmp(p,EndHdr,strlen(EndHdr))) break;
    if(!strncmp(p,EndAntenna,strlen(EndAntenna))) break;
    if(p[0]=='*') continue;
    if(strncmp(p,"ANT",3)) continue;
    if(k==MAX_ANTS)
    { fprintf(stderr,"MAX_ANTS %d exceeded\n",MAX_ANTS); return -1;}
    sscanf(p+3,"%d %*s %s %lf %lf %lf %lf %lf\n",&i,ant->name,&ant->bx,
	   &ant->by,&ant->bz,&ant->d0[0],&ant->d0[1]);
    if(i!=k){fprintf(stderr,"Antennas out of seq!\n"); return -1;}
    ant++;k++;
  }
  if(k!=MAX_ANTS && k!=(MAX_ANTS-2)) { fprintf(stderr,"WARNING: only %d antenna pars found\n",k); return -1;}
  // GSB MAX_ANTS = 32, GHB MAX_ANTS = 30

  return 0;
}
int get_samp_info(char *buf,VisInfo *vinfo)
{ char *p=buf,word[VALLEN];
  SampInfo *samp=vinfo->samp; 
  int   i,j,k,found=0;
  
  while((p=strchr(p,'\n')) != NULL)
  { p++; if(!strncmp(p,EndHdr,strlen(EndHdr))) break;
    if(strncmp(p,StartSampler,strlen(StartSampler))) continue;
    else {found=1; break;}
  }
  if(!found){fprintf(stderr,"No %s section!\n",StartSampler); return -1;}

  for(k=0;k<MAX_SAMPS;k++)
  { samp[k].ant[0]=samp[k].band[0]='\0';
    vinfo->samp_num[k] = MAX_SAMPS;
  }
  j=k=0;
  while((p=strchr(p,'\n')) != NULL)
  { p++; if(!strncmp(p,EndHdr,strlen(EndHdr))) break;
    if(!strncmp(p,EndSampler,strlen(EndSampler))) break;
    if(p[0]=='*') continue;
    if(strncmp(p,"SMP",3)) continue;
    if(k==MAX_SAMPS) // MAX_SAMPS for GHB =120, GSB=64
    { fprintf(stderr,"MAX_SAMPS %d exceeded\n",MAX_SAMPS); return -1;}
    sscanf(p+VALOFF,"%s",word);
    if(strncmp(word,"None",4))
    { sscanf(p+3,"%d %*s%s %s %d\n",&i,samp->ant,samp->band,&samp->fft_id);
      if(i!=k){fprintf(stderr,"Samplers out of seq!\n"); return -1;}
      vinfo->samp_num[j++] = i; /* read here (SAMP_NUM keyword depreciated) */
    }
    samp++;k++;
  }

  return 0;
}
int get_base_info(char *buf,BaselineType *baseline)
{ char *p=buf;
  BaselineType *base=baseline; 
  int   i,k,found=0;
  
  while((p=strchr(p,'\n')) != NULL)
  { p++; if(!strncmp(p,EndHdr,strlen(EndHdr))) break;
    if(strncmp(p,StartBaseline,strlen(StartBaseline))) continue;
    else {found=1; break;}
  }
  if(!found){fprintf(stderr,"No %s section!\n",StartBaseline); return -1;}

  k=0;
  while((p=strchr(p,'\n')) != NULL)
  { p++; if(!strncmp(p,EndHdr,strlen(EndHdr))) break;
    if(!strncmp(p,EndBaseline,strlen(EndBaseline))) break;
    if(p[0]=='*') continue;
    if(strncmp(p,"BAS",3)) continue;
    if(k==MAX_BASE)
    { fprintf(stderr,"MAX_BASE %d exceeded\n",MAX_BASE); return -1;}
    sscanf(p+3,"%d %*s %d %d %d %d %d %d %s %s %s %s\n",&i,&base->ant[0],
	   &base->band[0],&base->ant[1],&base->band[1],&base->samp[0],
	   &base->samp[1],base->antname[0],base->bandname[0],
	   base->antname[1],base->bandname[1]);
	   
    if(i!=k){fprintf(stderr,"BaseLines out of seq!\n"); return -1;}
    base++;k++;
  }

  return 0;
}

int corr4_to_corr_ref(LtaInfo *linfo)
{ /* GSB to CORR REF*/
  CorrType corr1,*corr=&linfo->corr;
  unsigned long off;
  char *in_p=(char*)corr,*out_p=(char*)&corr1;

  float *mjd_ref;
  double *dmjd_ref;
  int *gen;//,i;
  short *gen1;
  mjd_ref=(float*)malloc(4);
  dmjd_ref=(double*)malloc(8);
  gen=(int*)malloc(4);
//  in_p+=2; // four bytes needs to be skipped after chan num for alignment.
  gen1=(short int*)malloc(2);

#ifdef DEBUG_MODE
  fprintf(stdout, "Calling From CORR4-to-CORR REF\n");
#endif
  memcpy(&corr1,corr,sizeof(CorrType)); 
  off = 8+NAMELEN+8*MAX_BANDS;
  in_p +=off;  out_p +=off;
  memcpy(out_p,in_p,MAX_ANTS*sizeof(AntennaParType)); /* old MAX_ANTS */

  in_p  += MAX_ANTS*sizeof(AntennaParType);
  out_p += MAX_ANTS*sizeof(AntennaParType);

  memcpy(out_p,in_p,64*sizeof(SamplerType)); /* GSB MAX_SAMPS = 64 */
  in_p  += 64*sizeof(SamplerType);
  out_p += MAX_SAMPS*sizeof(SamplerType);  /* MAX SAMPS = 120 defined in master corr */

  memcpy(out_p,in_p,MAX_BASE*sizeof(BaseParType)); /* old_MAX_BASE */
  in_p  += MAX_BASE*sizeof(BaseParType);
  out_p += MAX_BASE*sizeof(BaseParType);

  memcpy(out_p,in_p,sizeof(CorrParType));

#ifdef DEBUG_MODE
  memcpy(gen, in_p, 4);
  fprintf(stdout, "===== >>> MACS       = %d\n", *gen);
  memcpy(gen, in_p+4, 4);
  fprintf(stdout, "===== >>> CHANNELS   = %d\n", *gen);
  memcpy(gen, in_p+8, 4);
  fprintf(stdout, "===== >>> POLS       = %d\n", *gen);
  memcpy(gen, in_p+12, 4);
  fprintf(stdout, "===== >>> STA        = %d\n", *gen);
  memcpy(gen, in_p+16, 4);
  fprintf(stdout, "===== >>> CONTROL    = %d\n", *gen);
  memcpy(gen, in_p+20, 4);
  fprintf(stdout, "===== >>> STATIME    = %d\n", *gen);
  memcpy(gen, in_p+24, 4);
  fprintf(stdout, "===== >>> IABEAM     = %d\n", *gen);
  memcpy(gen, in_p+28, 4);
  fprintf(stdout, "===== >>> PABEAM1    = %d\n", *gen);
  memcpy(gen, in_p+32, 4);
  fprintf(stdout, "===== >>> PABEAM2    = %d\n", *gen);
  memcpy(mjd_ref, in_p+48, 4);
  fprintf(stdout, "===== >>> F_STEP     = %f\n", *mjd_ref);
  memcpy(mjd_ref, in_p+52, 4);
  fprintf(stdout, "===== >>> CLOCK      = %f\n\n", *mjd_ref);
#endif

  in_p  += sizeof(CorrParType);
  out_p += sizeof(CorrParType);

#ifdef ARCH_64BIT
  out_p += 4; // need to skip 4 bytes for 64 bit machine.
#else 
  out_p += 0; // need to skip 4 bytes for 64 bit machine.
#endif

  memcpy(out_p, in_p, 7*sizeof(int)); // antmask, sampler, baselines of daspartype.
#ifdef DEBUG_MODE
  memcpy(gen, in_p, 4);
  fprintf(stdout, "===== >>> ANTMASK    = %x\n", *gen);
  memcpy(gen, in_p+4, 4);
  fprintf(stdout, "===== >>> SAMPLERS   = %d\n", *gen);
  memcpy(gen, in_p+8, 4);
  fprintf(stdout, "===== >>> BASELINES  = %d\n", *gen);
  memcpy(gen, in_p+12, 4);
  fprintf(stdout, "===== >>> CHANNELS  = %d\n", *gen);
  memcpy(gen, in_p+16, 4);
  fprintf(stdout, "===== >>> LTA       = %d\n", *gen);
  memcpy(gen, in_p+20, 4);
  fprintf(stdout, "===== >>> GSBMAXCHA = %d\n", *gen);
  memcpy(gen, in_p+24, 4);
  fprintf(stdout, "===== >>> GSB FSTOP = %d\n", *gen);
#endif
  in_p += (7*sizeof(int));
  out_p += (7*sizeof(int));

  memcpy(out_p, in_p, 3*sizeof(short)); // bandmask, mode, gsb_stokes;
#ifdef DEBUG_MODE
  memcpy(gen1, in_p, 2);
  fprintf(stdout, "===== >>> BANDMASK  = %d\n", *gen1);
  memcpy(gen1, in_p+2, 2);
  fprintf(stdout, "===== >>> MODE      = %d\n", *gen1);
  memcpy(gen1, in_p+4, 2);
  fprintf(stdout, "===== >>> GSB_STOKES= %d\n", *gen1);
#endif
  in_p += (3*sizeof(short));
  out_p += (3*sizeof(short));

  memcpy(out_p, in_p, MAX_CHANS*sizeof(short)); //chan_num[MAX_CHANS]
#ifdef DEBUG_MODE
//  for (i=0;i<MAX_CHANS;i++) {
//  memcpy(gen1, in_p+2*i, 2);
//  fprintf(stdout, "===== >>> CHAN[%4d] = %d\n", i, *gen1);
//  } 
#endif
  in_p += (MAX_CHANS*sizeof(short));
  out_p += (MAX_CHANS*sizeof(short));

#ifdef ARCH_64BIT
  in_p+=2; // four bytes needs to be skipped after chan num for alignment.
  out_p+=6; // four bytes needs to be skipped after chan num for alignment.
#else
  in_p+=2; // four bytes needs to be skipped after chan num for alignment.
  out_p+=2; // four bytes needs to be skipped after chan num for alignment.
#endif

  memcpy(out_p, in_p, 4*sizeof(double)); // mjd_ref, t_unit...
#ifdef DEBUG_MODE
  memcpy(dmjd_ref, in_p, 8);
  fprintf(stdout, "===== >>> MJD_REF   = %lf\n", *dmjd_ref);
  memcpy(dmjd_ref, in_p+8, 8);
  fprintf(stdout, "===== >>> T_UNIT   = %lf\n", *dmjd_ref);
  memcpy(dmjd_ref, in_p+16, 8);
  fprintf(stdout, "===== >>> GSBACQBW = %lf\n", *dmjd_ref);
  memcpy(dmjd_ref, in_p+24, 8);
  fprintf(stdout, "===== >>> GSBFINAL = %lf\n", *dmjd_ref);
#endif
  in_p += (4*sizeof(double));
  out_p += (4*sizeof(double));

//  memcpy(out_p, in_p, sizeof(double)); // dut1
#ifdef DEBUG_MODE
  memcpy(dmjd_ref, in_p, 8);
  fprintf(stdout, "===== >>> DUT1      = %lf\n", *dmjd_ref);
#endif
  in_p += (sizeof(double));
  out_p += (sizeof(double));

//  memcpy(out_p,in_p,sizeof(DasParType)); // daspartype is copied element by element.

  memcpy(corr,&corr1,sizeof(CorrType));
  memcpy(corr->version,"COR30x2 ",8);
  strcpy(linfo->lhdr.version,corr->version);

#ifdef DEBUG_MODE
  fprintf(stdout, "\n========= >>>>>>>>>  ANTENNA   = %s\n", corr->antenna[31].name);
   fprintf(stdout, "========= >>>>>>>>>  ANT BX    = %lf\n", corr->antenna[31].bx);
   fprintf(stdout, "========= >>>>>>>>>  ANT BY    = %lf\n", corr->antenna[31].by);
   fprintf(stdout, "========= >>>>>>>>>  ANT BZ    = %lf\n", corr->antenna[31].bz);
   fprintf(stdout, "========= >>>>>>>>>  ANT DLY   = %lf\n", corr->antenna[31].d0[0]);

  fprintf(stdout, "\n=========== >>>>>>>>>  SAMP ANT  = %d\n", corr->sampler[62].ant_id);
    fprintf(stdout, "=========== >>>>>>>>>  SAMP BAND = %d\n", corr->sampler[62].band);
    fprintf(stdout, "=========== >>>>>>>>>  SAMP FFT  = %d\n", corr->sampler[62].fft_id);
    fprintf(stdout, "=========== >>>>>>>>>  SAMP DPC  = %d\n", corr->sampler[62].dpc);

  fprintf(stdout, "\n========= >>>>>>>>>  MACS      = %d\n", corr->corrpar.macs);
  fprintf(stdout, "========= >>>>>>>>>  CHANNELS  = %d\n", corr->corrpar.channels);
  fprintf(stdout, "========= >>>>>>>>>  POLS      = %d\n", corr->corrpar.pols);
  fprintf(stdout, "========= >>>>>>>>>  STA       = %d\n", corr->corrpar.sta);
  fprintf(stdout, "========= >>>>>>>>>  CNTRL     = %d\n", corr->corrpar.cntrl);
  fprintf(stdout, "========= >>>>>>>>>  STATIME   = %d\n", corr->corrpar.statime);
  fprintf(stdout, "========= >>>>>>>>>  F_STEP    = %lf\n", corr->corrpar.f_step);
  fprintf(stdout, "========= >>>>>>>>>  CLOCK     = %lf\n", corr->corrpar.clock);

  fprintf(stdout, "\n========= >>>>>>>>>  ANTMASK   = %x\n", corr->daspar.antmask);
  fprintf(stdout, "========= >>>>>>>>>  SAMPLERS  = %d\n", corr->daspar.samplers);
  fprintf(stdout, "========= >>>>>>>>>  BASELINES = %d\n", corr->daspar.baselines);
  fprintf(stdout, "========= >>>>>>>>>  CHANNELS  = %d\n", corr->daspar.channels);
  fprintf(stdout, "========= >>>>>>>>>  LTA       = %d\n", corr->daspar.lta);
  fprintf(stdout, "========= >>>>>>>>>  GSBMAXCHAN= %d\n", corr->daspar.gsb_maxchan);
  fprintf(stdout, "========= >>>>>>>>>  GSB FSTOP = %d\n", corr->daspar.gsb_fstop);
  fprintf(stdout, "========= >>>>>>>>>  BANDMASK  = %d\n", corr->daspar.bandmask);
  fprintf(stdout, "========= >>>>>>>>>  MODE      = %d\n", corr->daspar.mode);
  fprintf(stdout, "========= >>>>>>>>>  GSB_STOKES= %d\n", corr->daspar.gsb_stokes);
  fprintf(stdout, "========= >>>>>>>>>  CHAN[127 ]= %d\n", corr->daspar.chan_num[127]);
  fprintf(stdout, "========= >>>>>>>>>  MJD_REF   = %lf\n", corr->daspar.mjd_ref);
  fprintf(stdout, "========= >>>>>>>>>  T_UNIT    = %f\n", corr->daspar.t_unit);
  fprintf(stdout, "========= >>>>>>>>>  GSBACQBW  = %f\n", corr->daspar.gsb_acq_bw);
  fprintf(stdout, "========= >>>>>>>>>  GSBFINALBW= %f\n", corr->daspar.gsb_final_bw);
//  fprintf(stdout, "========= >>>>>>>>>  DUT1      = %lf\n", corr->daspar.dut1);
  fprintf(stdout, "========= CORR4 to CORR REF ============  \n");
#endif

  return 0;
}

int corr2_to_corr_ref(LtaInfo *linfo)
{ /* post polar mode + GHB to post polar mode + GSB*/
  CorrType corr1,*corr=&linfo->corr;
  unsigned long off;
  char *in_p=(char*)corr,*out_p=(char*)&corr1;

  float *mjd_ref;
  double *dmjd_ref;
  int *gen;//,i;
  short *gen1;
  mjd_ref=(float*)malloc(4);
  dmjd_ref=(double*)malloc(8);
  gen=(int*)malloc(4);
  gen1=(short int*)malloc(2);

#ifdef DEBUG_MODE
  fprintf(stdout, "Calling From CORR2-to-CORR REF\n");
#endif
  memcpy(&corr1,corr,sizeof(CorrType)); 
  off = 8+NAMELEN+8*MAX_BANDS;
  in_p +=off;  out_p +=off;

  memcpy(out_p,in_p,30*sizeof(AntennaParType)); /* old MAX_ANTS */
  in_p  += 30*sizeof(AntennaParType);
  out_p += MAX_ANTS*sizeof(AntennaParType);

  memcpy(out_p,in_p,MAX_SAMPS*sizeof(SamplerType)); /* GSB MAX_SAMPS = 64 */
  in_p  += 60*sizeof(SamplerType);
  out_p += MAX_SAMPS*sizeof(SamplerType);  /* MAX SAMPS = 120 defined in master corr */

  memcpy(out_p,in_p,32*MAC_CARDS*sizeof(BaseParType)); /* old_MAX_BASE */
  in_p  += 32*MAC_CARDS*sizeof(BaseParType);
  out_p += MAX_BASE*sizeof(BaseParType);

  memcpy(out_p,in_p,16); /* CorrParType macs, channels, pols, sta */

#ifdef DEBUG_MODE
  memcpy(gen,in_p,4); /* CorrParType macs, channels, pols, sta */
  fprintf(stdout, "===== >>> MACS         = %d\n", *gen);
  memcpy(gen,in_p+4,4); /* CorrParType macs, channels, pols, sta */
  fprintf(stdout, "===== >>> CHANNELS     = %d\n", *gen);
  memcpy(gen,in_p+8,4); /* CorrParType macs, channels, pols, sta */
  fprintf(stdout, "===== >>> POLS         = %d\n", *gen);
  memcpy(gen,in_p+12,4); /* CorrParType macs, channels, pols, sta */
  fprintf(stdout, "===== >>> STA          = %d\n", *gen);
#endif

#ifdef ARCH_64BIT
  in_p += 16;
  out_p += 20;
#else
  in_p += 16;
  out_p += 20;
#endif
  memcpy(out_p,in_p,4); /* CorrParType statime */

#ifdef DEBUG_MODE
  memcpy(gen,in_p,4); /* CorrParType statime */
  fprintf(stdout, "===== >>> STATIME      = %d\n", *gen);
#endif

#ifdef ARCH_64BIT
  in_p += 4;
  out_p += 28;
#else
  in_p += 4;
  out_p += 28;
#endif
  memcpy(out_p,in_p,12); /* CorrParType f_step, clock  dpcmux,clksel,fftmode,macmode */

#ifdef DEBUG_MODE
  memcpy(mjd_ref,in_p,4); /* CorrParType f_step, clock  dpcmux,clksel,fftmode,macmode */
  fprintf(stdout, "===== >>> F_STEP       = %f\n", *mjd_ref);
  memcpy(mjd_ref,in_p+4,4); /* CorrParType f_step, clock  dpcmux,clksel,fftmode,macmode */
  fprintf(stdout, "===== >>> CLOCK        = %f\n", *mjd_ref);
#endif
  
#ifdef ARCH_64BIT
  in_p += 12;
  out_p += 16;
#else
  in_p += 12;
  out_p += 12;
#endif
  memcpy(out_p,in_p,20); /* DasParType antmask, samplers, baselines, channels, lta */

#ifdef DEBUG_MODE
  memcpy(gen,in_p,4); 
  fprintf(stdout, "===== >>> ANTMASK      = %x\n", *gen);
  memcpy(gen,in_p+4,4); 
  fprintf(stdout, "===== >>> SAMPLERS     = %d\n", *gen);
  memcpy(gen,in_p+8,4); 
  fprintf(stdout, "===== >>> BASELINES    = %d\n", *gen);
  memcpy(gen,in_p+12,4); 
  fprintf(stdout, "===== >>> CHANNELS     = %d\n", *gen);
  memcpy(gen,in_p+16,4); 
  fprintf(stdout, "===== >>> LTA          = %d\n", *gen);
#endif

#ifdef ARCH_64BIT
  in_p += 20;
  out_p += 28;
#else
  in_p += 20;
  out_p += 28;
#endif
  memcpy(out_p,in_p,4); /* DasParType bandmask, mode */

#ifdef DEBUG_MODE
  memcpy(gen1,in_p,4); 
  fprintf(stdout, "===== >>> BANDMASK     = %d\n", *gen1);
  memcpy(gen1,in_p+4,4); 
  fprintf(stdout, "===== >>> MODE         = %d\n", *gen1);
#endif

#ifdef ARCH_64BIT
  in_p += 4;
  out_p += 6;
#else
  in_p += 4;
  out_p += 6;
#endif
  memcpy(out_p,in_p,256); // DasParType chan_num, MAX = 256 channel for new ghb corr, 1024 for gsb 
#ifdef DEBUG_MODE
  memcpy(gen1,in_p+20,4);  // channel no. 10
  fprintf(stdout, "===== >>> CH[10]        = %d\n", *gen1);
#endif

#ifdef ARCH_64BIT
  in_p += 256;
  out_p += 2048;
  out_p += 6; // Two bytes extra (?), need to match MJD_REF in right place in out_p.
#else
  in_p += 256;
  out_p += 2048;
  out_p += 2; // Two bytes extra (?), need to match MJD_REF in right place in out_p.
#endif
  memcpy(out_p,in_p,16); /* mjd_ref and t_unit */

#ifdef DEBUG_MODE
  memcpy(dmjd_ref, in_p, 8);
  fprintf(stdout, "===== >>> MJD_REF      = %lf\n", *dmjd_ref);
  memcpy(dmjd_ref, in_p+8, 8);
  fprintf(stdout, "===== >>> T_UNIT       = %lf\n", *dmjd_ref);
#endif

  memcpy(corr,&corr1,sizeof(CorrType));

// for(i=0;i<1024;i++) fprintf(stdout, "========= >>>>>>>>>  CHAN[%04d] %d\n", i, corr->daspar.chan_num[i]);

  memcpy(corr->version,"CORR30  ",8);
  strcpy(linfo->lhdr.version,corr->version);

#ifdef DEBUG_MODE
  fprintf(stdout, "\n========= >>>>>>>>>  ANTENNA   = %s\n", corr->antenna[29].name);
   fprintf(stdout, "========= >>>>>>>>>  ANT BX    = %lf\n", corr->antenna[29].bx);
   fprintf(stdout, "========= >>>>>>>>>  ANT BY    = %lf\n", corr->antenna[29].by);
   fprintf(stdout, "========= >>>>>>>>>  ANT BZ    = %lf\n", corr->antenna[29].bz);
   fprintf(stdout, "========= >>>>>>>>>  ANT DLY   = %lf\n", corr->antenna[29].d0[0]);

  fprintf(stdout, "\n=========== >>>>>>>>>  SAMP ANT  = %d\n", corr->sampler[59].ant_id);
    fprintf(stdout, "=========== >>>>>>>>>  SAMP BAND = %d\n", corr->sampler[59].band);
    fprintf(stdout, "=========== >>>>>>>>>  SAMP FFT  = %d\n", corr->sampler[59].fft_id);
    fprintf(stdout, "=========== >>>>>>>>>  SAMP DPC  = %d\n", corr->sampler[59].dpc);

  fprintf(stdout, "\n========= >>>>>>>>>  ANTMASK   = %x\n", corr->daspar.antmask);
  fprintf(stdout, "========= >>>>>>>>>  SAMPLERS  = %d\n", corr->daspar.samplers);
  fprintf(stdout, "========= >>>>>>>>>  BASELINES = %d\n", corr->daspar.baselines);
  fprintf(stdout, "========= >>>>>>>>>  CHANNELS  = %d\n", corr->daspar.channels);
  fprintf(stdout, "========= >>>>>>>>>  LTA       = %d\n", corr->daspar.lta);
  fprintf(stdout, "========= >>>>>>>>>  BANDMASK  = %d\n", corr->daspar.bandmask);
  fprintf(stdout, "========= >>>>>>>>>  MODE      = %d\n", corr->daspar.mode);
  fprintf(stdout, "========= >>>>>>>>>  CHAN[50  ]= %d\n", corr->daspar.chan_num[50]);
  fprintf(stdout, "========= >>>>>>>>>  MJD_REF   = %lf\n", corr->daspar.mjd_ref);
  fprintf(stdout, "========= >>>>>>>>>  T_UNIT    = %f\n", corr->daspar.t_unit);
//  fprintf(stdout, "========= >>>>>>>>>  DUT1      = %15.12lf \n", (corr->daspar.dut1)*(12.0*3600/M_PI));

  fprintf(stdout, "\n========= >>>>>>>>>  MACS      = %d\n", corr->corrpar.macs);
  fprintf(stdout, "========= >>>>>>>>>  CHANNELS  = %d\n", corr->corrpar.channels);
  fprintf(stdout, "========= >>>>>>>>>  POLS      = %d\n", corr->corrpar.pols);
  fprintf(stdout, "========= >>>>>>>>>  STA       = %d\n", corr->corrpar.sta);
  fprintf(stdout, "========= >>>>>>>>>  STATIME   = %d\n", corr->corrpar.statime);
  fprintf(stdout, "========= >>>>>>>>>  F_STEP    = %lf\n", corr->corrpar.f_step);
  fprintf(stdout, "========= >>>>>>>>>  CLOCK     = %lf\n", corr->corrpar.clock);

#endif

  return 0;
}

int corr3_to_corr_ref(LtaInfo *linfo)
{ /* post polar mode + GHB to post polar mode + GSB*/
  CorrType corr1,*corr=&linfo->corr;
  unsigned long off;
  char *in_p=(char*)corr,*out_p=(char*)&corr1;

//=================================================================
  float *mjd_ref;
  double *dmjd_ref;
  int *gen;//,i;
  short *gen1;

  mjd_ref=(float*)malloc(sizeof(float));
  dmjd_ref=(double*)malloc(sizeof(double));
  gen=(int*)malloc(sizeof(int));
  gen1=(short int*)malloc(sizeof(short));
//=================================================================

#ifdef DEBUG_MODE
  fprintf(stdout, "Calling From CORR3-to-CORR REF\n");
#endif
  memcpy(&corr1,corr,sizeof(CorrType)); 
  off = 8+NAMELEN+8*MAX_BANDS;
  in_p +=off;  out_p +=off;
  memcpy(out_p,in_p,30*sizeof(AntennaParType)); /* old MAX_ANTS */

  in_p  += 30*sizeof(AntennaParType);
  out_p += MAX_ANTS*sizeof(AntennaParType);

  memcpy(out_p,in_p,MAX_SAMPS*sizeof(SamplerType)); /* GSB MAX_SAMPS = 64 */

  in_p  += MAX_SAMPS*sizeof(SamplerType);
  out_p += MAX_SAMPS*sizeof(SamplerType);  /* MAX SAMPS = 120 defined in master corr */

  memcpy(out_p,in_p,MAX_BASE*sizeof(BaseParType)); /* old_MAX_BASE */

  in_p  += MAX_BASE*sizeof(BaseParType);
  out_p += MAX_BASE*sizeof(BaseParType);

  memcpy(out_p,in_p,16); /* CorrParType macs, channels, pols, sta */

#ifdef DEBUG_MODE
  memcpy(gen,in_p,4); /* CorrParType macs, channels, pols, sta */
  fprintf(stdout, "===== >>> MACS         = %d\n", *gen);
  memcpy(gen,in_p+4,4); /* CorrParType macs, channels, pols, sta */
  fprintf(stdout, "===== >>> CHANNELS     = %d\n", *gen);
  memcpy(gen,in_p+8,4); /* CorrParType macs, channels, pols, sta */
  fprintf(stdout, "===== >>> POLS         = %d\n", *gen);
  memcpy(gen,in_p+12,4); /* CorrParType macs, channels, pols, sta */
  fprintf(stdout, "===== >>> STA          = %d\n", *gen);
#endif

  in_p += 16;
  out_p += 20;
  memcpy(out_p,in_p,4); /* CorrParType statime */

#ifdef DEBUG_MODE
  memcpy(gen,in_p,4); /* CorrParType macs, channels, pols, sta */
  fprintf(stdout, "===== >>> STATIME      = %d\n", *gen);
#endif

  in_p += 4;
#ifdef ARCH_64BIT
  out_p += 28;
#else
  out_p += 28;
#endif
  memcpy(out_p,in_p,12); /* CorrParType f_step, clock  dpcmux,clksel,fftmode,macmode */

#ifdef DEBUG_MODE
  memcpy(mjd_ref,in_p,4); /* CorrParType f_step, clock  dpcmux,clksel,fftmode,macmode */
  fprintf(stdout, "===== >>> F_STEP       = %f\n", *mjd_ref);
  memcpy(mjd_ref,in_p+4,4); /* CorrParType f_step, clock  dpcmux,clksel,fftmode,macmode */
  fprintf(stdout, "===== >>> CLOCK        = %f\n", *mjd_ref);
#endif
  
  in_p += 12;
#ifdef ARCH_64BIT
  out_p += 16;
#else
  out_p += 12;
#endif
  memcpy(out_p,in_p,20); /* DasParType antmask, samplers, baselines, channels, lta */

#ifdef DEBUG_MODE
  memcpy(gen,in_p+4,4); 
  fprintf(stdout, "===== >>> SAMPLERS     = %d\n", *gen);
  memcpy(gen,in_p+8,4); 
  fprintf(stdout, "===== >>> BASELINES    = %d\n", *gen);
  memcpy(gen,in_p+12,4); 
  fprintf(stdout, "===== >>> CHANNELS     = %d\n", *gen);
  memcpy(gen,in_p+16,4); 
  fprintf(stdout, "===== >>> LTA          = %d\n", *gen);
#endif

  in_p += 20;
  out_p += 28;
  memcpy(out_p,in_p,4); /* DasParType bandmask, mode */

#ifdef DEBUG_MODE
  memcpy(gen1,in_p,4); 
  fprintf(stdout, "===== >>> BANDMASK     = %d\n", *gen1);
  memcpy(gen1,in_p+4,4); 
  fprintf(stdout, "===== >>> MODE         = %d\n", *gen1);
#endif

  in_p += 4;
  out_p += 6;
  memcpy(out_p,in_p,512); /* DasParType chan_num, MAX = 256 channel for new ghb corr, 1024 for gsb */
#ifdef DEBUG_MODE
//  int i;
//  for (i=0;i<128;i++) {
//  memcpy(gen1, in_p+2*i, 2);
//  fprintf(stdout, "===== >>> CHAN[%4d] = %d\n", i, *gen1);
//  } 
#endif

  in_p += 512;
  out_p += 2048;
#ifdef ARCH_64BIT
  out_p += 6; // Two bytes extra (?), need to match MJD_REF in right place in out_p.
#else
  out_p += 2; // Two bytes extra (?), need to match MJD_REF in right place in out_p.
#endif
  memcpy(out_p,in_p,16); /* mjd_ref and t_unit */

#ifdef DEBUG_MODE
  memcpy(dmjd_ref, in_p, 8);
  fprintf(stdout, "===== >>> MJD_REF      = %lf\n", *dmjd_ref);
  memcpy(dmjd_ref, in_p+8, 8);
  fprintf(stdout, "===== >>> T_UNIT       = %lf\n", *dmjd_ref);
#endif

  in_p += 16;
  out_p += 32;
  memcpy(out_p,in_p,8); /* dut1 */

#ifdef DEBUG_MODE
  memcpy(dmjd_ref, in_p, 8);
  fprintf(stdout, "===== >>> DUT1         = %15.12lf\n", (*dmjd_ref)*(12.0*3600/M_PI));
#endif

  memcpy(corr,&corr1,sizeof(CorrType));

//  for(i=0;i<1024;i++) fprintf(stdout, "========= >>>>>>>>>  CHAN[%04d] %d\n", i, corr->daspar.chan_num[i]);

  memcpy(corr->version,"COR30x2 ",8);
//  strcpy(linfo->lhdr.version,corr->version); // this was causing ltahdr version of ltamerged file to change HST03 to HST00

#ifdef DEBUG_MODE
  fprintf(stdout, "\n========= >>>>>>>>>  ANTENNA   = %s\n", corr->antenna[29].name);
   fprintf(stdout, "========= >>>>>>>>>  ANT BX    = %lf\n", corr->antenna[29].bx);
   fprintf(stdout, "========= >>>>>>>>>  ANT BY    = %lf\n", corr->antenna[29].by);
   fprintf(stdout, "========= >>>>>>>>>  ANT BZ    = %lf\n", corr->antenna[29].bz);
   fprintf(stdout, "========= >>>>>>>>>  ANT DLY   = %lf\n", corr->antenna[29].d0[0]);

  fprintf(stdout, "\n=========== >>>>>>>>>  SAMP ANT  = %d\n", corr->sampler[30].ant_id);
    fprintf(stdout, "=========== >>>>>>>>>  SAMP BAND = %d\n", corr->sampler[30].band);
    fprintf(stdout, "=========== >>>>>>>>>  SAMP FFT  = %d\n", corr->sampler[30].fft_id);
    fprintf(stdout, "=========== >>>>>>>>>  SAMP DPC  = %d\n", corr->sampler[30].dpc);

  fprintf(stdout, "\n========= >>>>>>>>>  ANTMASK   = %x\n", corr->daspar.antmask);
  fprintf(stdout, "========= >>>>>>>>>  SAMPLERS  = %d\n", corr->daspar.samplers);
  fprintf(stdout, "========= >>>>>>>>>  BASELINES = %d\n", corr->daspar.baselines);
  fprintf(stdout, "========= >>>>>>>>>  CHANNELS  = %d\n", corr->daspar.channels);
  fprintf(stdout, "========= >>>>>>>>>  LTA       = %d\n", corr->daspar.lta);
  fprintf(stdout, "========= >>>>>>>>>  BANDMASK  = %d\n", corr->daspar.bandmask);
  fprintf(stdout, "========= >>>>>>>>>  MODE      = %d\n", corr->daspar.mode);
  fprintf(stdout, "========= >>>>>>>>>  CHAN[50  ]= %d\n", corr->daspar.chan_num[50]);
  fprintf(stdout, "========= >>>>>>>>>  MJD_REF   = %lf\n", corr->daspar.mjd_ref);
  fprintf(stdout, "========= >>>>>>>>>  T_UNIT    = %f\n", corr->daspar.t_unit);
//  fprintf(stdout, "========= >>>>>>>>>  DUT1      = %15.12lf \n", (corr->daspar.dut1)*(12.0*3600/M_PI));

  fprintf(stdout, "\n========= >>>>>>>>>  MACS      = %d\n", corr->corrpar.macs);
  fprintf(stdout, "========= >>>>>>>>>  CHANNELS  = %d\n", corr->corrpar.channels);
  fprintf(stdout, "========= >>>>>>>>>  POLS      = %d\n", corr->corrpar.pols);
  fprintf(stdout, "========= >>>>>>>>>  STA       = %d\n", corr->corrpar.sta);
  fprintf(stdout, "========= >>>>>>>>>  STATIME   = %d\n", corr->corrpar.statime);
  fprintf(stdout, "========= >>>>>>>>>  F_STEP    = %lf\n", corr->corrpar.f_step);
  fprintf(stdout, "========= >>>>>>>>>  CLOCK     = %lf\n", corr->corrpar.clock);

#endif

  return 0;
}

int corr_ref_to_corr2(LtaInfo *linfo,CorrType *corr1)
{ /* CORR REF to pre polar mode + GHB*/
  CorrType      *corr=&linfo->corr;
  char          *in_p=(char*)corr,*out_p=(char*)corr1;
  int            i;
  memcpy(corr1,corr,sizeof(CorrType));
  in_p += 8+NAMELEN+8*MAX_BANDS;
  out_p += 8+NAMELEN+8*MAX_BANDS;

  memcpy(out_p, in_p, 30*sizeof(AntennaParType));
  out_p += 30*sizeof(AntennaParType); /* old MAX_ANTS */

  for(i=0;i<corr->daspar.samplers;i++)
  { SamplerType *samp=corr->sampler+i;
    if(samp->dpc!=MAX_SAMPS)
      memcpy(out_p,samp,sizeof(SamplerType));
    out_p+=sizeof(SamplerType);
  }
  for(i=corr->daspar.samplers;i<60;i++) /* old MAX_SAMP */
  { SamplerType *samp=(SamplerType*)out_p;
    samp->dpc=MAX_SAMPS;
    out_p +=sizeof(SamplerType);
  }
  for(i=0;i<corr->daspar.baselines;i++)
  { BaseParType *base=corr->baseline+i;
    if(base->samp[0].dpc==MAX_SAMPS||base->samp[1].dpc==MAX_SAMPS)
      continue;
    memcpy(out_p,base,sizeof(BaseParType));
    out_p+=sizeof(BaseParType);
  }
  for(i=corr->daspar.baselines;i<32*MAC_CARDS;i++) /* old MAX_BASE */
  { BaseParType *base=(BaseParType*)out_p;
    base->samp[0].dpc=base->samp[1].dpc=MAX_SAMPS;
    out_p +=sizeof(BaseParType);
  }

  memcpy(out_p, &corr->corrpar, 4*sizeof(int)); //macs, channels, pols, sta
  out_p+=4*sizeof(int);
  memcpy(out_p, &corr->corrpar.statime, sizeof(int)); //statime
  out_p+=sizeof(int);
  memcpy(out_p, &corr->corrpar.f_step, 2*sizeof(float)); //f_step, clock
  out_p+=2*sizeof(float);
  memcpy(out_p,&corr->corrpar.dpcmux,4*sizeof(char)); //dpcmux,clksel,fftmode,macmode
  out_p+=4*sizeof(char);
  memcpy(out_p,&corr->daspar.antmask, 5*sizeof(int)); //antmask, samplers, baselines, channels, lta
  out_p+=5*sizeof(int);
  memcpy(out_p,&corr->daspar.bandmask, 2*sizeof(short)); //bandmask, mode
  out_p+=2*sizeof(short);
  memcpy(out_p,corr->daspar.chan_num,128*sizeof(short)); // old MAX_CHANS 
  out_p+=128*sizeof(short);
  memcpy(out_p,&corr->daspar.mjd_ref,2*sizeof(double)); // mjd_ref, t_unit
  strcpy(corr1->version,linfo->lhdr.vinfo.corrlite.old_version);
#ifdef DEBUG_MODE
  fprintf(stdout, "\n========= <<<<<<<<<  MACS      = %d\n", corr1->corrpar.macs);
  fprintf(stdout, "========= <<<<<<<<<  CHANNELS  = %d\n", corr1->corrpar.channels);
  fprintf(stdout, "========= <<<<<<<<<  POLS      = %d\n", corr1->corrpar.pols);
  fprintf(stdout, "========= <<<<<<<<<  STA       = %d\n", corr1->corrpar.sta);
  fprintf(stdout, "========= <<<<<<<<<  STATIME   = %d\n", corr1->corrpar.statime);
  fprintf(stdout, "========= <<<<<<<<<  F_STEP    = %lf\n", corr1->corrpar.f_step);
  fprintf(stdout, "========= <<<<<<<<<  CLOCK     = %lf\n", corr1->corrpar.clock);

  fprintf(stdout, "\n========= <<<<<<<<<  ANTMASK   = %x\n", corr1->daspar.antmask);
  fprintf(stdout, "========= <<<<<<<<<  SAMPLERS  = %d\n", corr1->daspar.samplers);
  fprintf(stdout, "========= <<<<<<<<<  LTA       = %d\n", corr1->daspar.lta);
  fprintf(stdout, "========= <<<<<<<<<  BASELINES = %d\n", corr1->daspar.baselines);
  fprintf(stdout, "========= <<<<<<<<<  CHANNELS  = %d\n", corr1->daspar.channels);
  fprintf(stdout, "========= <<<<<<<<<  MJD_REF   = %lf\n", corr1->daspar.mjd_ref);
  fprintf(stdout, "========= <<<<<<<<<  T_UNIT    = %f\n", corr1->daspar.t_unit);
#endif

  return 0;
}

int corr_ref_to_corr3(LtaInfo *linfo,CorrType *corr1)
{ /* CORR REF to pre polar mode + GHB*/
  CorrType      *corr=&linfo->corr;
  char          *in_p=(char*)corr,*out_p=(char*)corr1;
  int            i;
  memcpy(corr1,corr,sizeof(CorrType));
  in_p += 8+NAMELEN+8*MAX_BANDS;
  out_p += 8+NAMELEN+8*MAX_BANDS;

  memcpy(out_p, in_p, 30*sizeof(AntennaParType));
  out_p += 30*sizeof(AntennaParType); /* old MAX_ANTS */

  for(i=0;i<corr->daspar.samplers;i++)
  { SamplerType *samp=corr->sampler+i;
    if(samp->dpc!=MAX_SAMPS)
      memcpy(out_p,samp,sizeof(SamplerType));
    out_p+=sizeof(SamplerType);
  }
  for(i=corr->daspar.samplers;i<120;i++) /* old MAX_SAMP */
  { SamplerType *samp=(SamplerType*)out_p;
    samp->dpc=MAX_SAMPS;
    out_p +=sizeof(SamplerType);
  }
  for(i=0;i<corr->daspar.baselines;i++)
  { BaseParType *base=corr->baseline+i;
    if(base->samp[0].dpc==MAX_SAMPS||base->samp[1].dpc==MAX_SAMPS)
      continue;
    memcpy(out_p,base,sizeof(BaseParType));
    out_p+=sizeof(BaseParType);
  }
  for(i=corr->daspar.baselines;i<DAS_CARDS*32*MAC_CARDS;i++) /* old MAX_BASE */
  { BaseParType *base=(BaseParType*)out_p;
    base->samp[0].dpc=base->samp[1].dpc=MAX_SAMPS;
    out_p +=sizeof(BaseParType);
  }

  memcpy(out_p, &corr->corrpar, 4*sizeof(int)); //macs, channels, pols, sta
  out_p+=4*sizeof(int);
  memcpy(out_p, &corr->corrpar.statime, sizeof(int)); //statime
  out_p+=sizeof(int);
  memcpy(out_p, &corr->corrpar.f_step, 2*sizeof(float)); //f_step, clock
  out_p+=2*sizeof(float);
  memcpy(out_p,&corr->corrpar.dpcmux,4*sizeof(char)); //dpcmux,clksel,fftmode,macmode
  out_p+=4*sizeof(char);
  memcpy(out_p,&corr->daspar.antmask, 5*sizeof(int)); //antmask, samplers, baselines, channels, lta
  out_p+=5*sizeof(int);
  memcpy(out_p,&corr->daspar.bandmask, 2*sizeof(short)); //bandmask, mode
  out_p+=2*sizeof(short);
  memcpy(out_p,corr->daspar.chan_num,256*sizeof(short)); // old MAX_CHANS 
  out_p+=256*sizeof(short);
  memcpy(out_p,&corr->daspar.mjd_ref,2*sizeof(double)); // mjd_ref, t_unit
  out_p+=2*sizeof(double);
//  memcpy(out_p,&corr->daspar.dut1,sizeof(double)); // dut1
//  strcpy(corr1->version,linfo->lhdr.vinfo.corrlite.old_version); // overwriting version HST02 by old version.

#ifdef DEBUG_MODE
  fprintf(stdout, "\n========= <<<<<<<<<  SAMPLERS  = %d\n", corr1->daspar.samplers);
  fprintf(stdout, "========= <<<<<<<<<  LTA       = %d\n", corr1->daspar.lta);
  fprintf(stdout, "========= <<<<<<<<<  BASELINES = %d\n", corr1->daspar.baselines);
  fprintf(stdout, "========= <<<<<<<<<  CHANNELS  = %d\n", corr1->daspar.channels);
  fprintf(stdout, "========= <<<<<<<<<  MJD_REF   = %lf\n", corr1->daspar.mjd_ref);
  fprintf(stdout, "========= <<<<<<<<<  T_UNIT    = %f\n", corr1->daspar.t_unit);
  fprintf(stdout, "========= <<<<<<<<<  F_STEP    = %lf\n", corr1->corrpar.f_step);
  fprintf(stdout, "========= <<<<<<<<<  CLOCK     = %lf\n", corr1->corrpar.clock);
#endif

  return 0;
}

int corr_ref_to_corr4(LtaInfo *linfo,CorrType *corr1)
{ /* CORR REF to pre polar mode + GHB*/
  CorrType      *corr=&linfo->corr;
//  char          *in_p=(char*)corr,*out_p=(char*)corr1;
  char          *out_p=(char*)corr1;
  int            i;
  memcpy(corr1,corr,sizeof(CorrType));
  out_p += 8+NAMELEN+8*MAX_BANDS+MAX_ANTS*sizeof(AntennaParType);
//  fprintf(stdout, "SAMPLER TYPE SIZE = %ld SAMPLERS = %d BASEPARSIZE = %ld BASELINES=%d\n", sizeof(SamplerType), corr->daspar.samplers,sizeof(BaseParType), corr->daspar.baselines);
  for(i=0;i<corr->daspar.samplers;i++)
  { SamplerType *samp=corr->sampler+i;
    if(samp->dpc!=MAX_SAMPS)
      memcpy(out_p,samp,sizeof(SamplerType));
    out_p+=sizeof(SamplerType);
  }
  for(i=corr->daspar.samplers;i<64;i++) /* OLD GSB MAX_SAMP */
  { SamplerType *samp=(SamplerType*)out_p;
    samp->dpc=MAX_SAMPS;
    out_p +=sizeof(SamplerType);
  }

  for(i=0;i<corr->daspar.baselines;i++)
  { BaseParType *base=corr->baseline+i;
    if(base->samp[0].dpc==MAX_SAMPS||base->samp[1].dpc==MAX_SAMPS)
      continue;
    memcpy(out_p,base,sizeof(BaseParType));
    out_p+=sizeof(BaseParType);
  }
  for(i=corr->daspar.baselines;i<DAS_CARDS*32*MAC_CARDS;i++) /* GSB MAX_BASE */
  { BaseParType *base=(BaseParType*)out_p;
    base->samp[0].dpc=base->samp[1].dpc=MAX_SAMPS;
    out_p +=sizeof(BaseParType);
  }
  memcpy(out_p, &corr->corrpar, 9*sizeof(int)); //macs, channels, pols, sta
  out_p+=9*sizeof(int);
  memcpy(out_p, &corr->corrpar.iabeam_res, 5*sizeof(float)); //f_step, clock
  out_p+=5*sizeof(float);
  memcpy(out_p,&corr->corrpar.dpcmux,4*sizeof(char)); //dpcmux,clksel,fftmode,macmode
  out_p+=4*sizeof(char);
  memcpy(out_p,&corr->daspar.antmask, 7*sizeof(int)); //antmask, samplers, baselines, channels, lta
  out_p+=7*sizeof(int);
  memcpy(out_p,&corr->daspar.bandmask, 3*sizeof(short)); //bandmask, mode
  out_p+=3*sizeof(short);
  memcpy(out_p,corr->daspar.chan_num,1024*sizeof(short)); // old MAX_CHANS 
  out_p+=1024*sizeof(short);

#ifdef ARCH_64BIT
  out_p+=2;
#else
#endif
  memcpy(out_p,&corr->daspar.mjd_ref,4*sizeof(double)); // mjd_ref, t_unit
  out_p+=2*sizeof(double);

/*
#ifdef ARCH_64BIT
//    out_p+=224;
#else
#endif
  memcpy(out_p,&corr->corrpar,sizeof(CorrParType));
  out_p+=sizeof(CorrParType);
#ifdef ARCH_64BIT
//    out_p+=4;
#else
#endif
  memcpy(out_p,&corr->daspar,sizeof(DasParType));
*/
  strcpy(corr1->version,linfo->lhdr.vinfo.corrlite.old_version);
#ifdef DEBUG_MODE
  fprintf(stdout, "\n===================================\n");
  fprintf(stdout, "========= >>>>>>>>>  ANTENNA   = %s\n", corr1->antenna[31].name);
  fprintf(stdout, "========= >>>>>>>>>  ANT DLY   = %lf\n", corr1->antenna[31].d0[0]);
  fprintf(stdout, "========= >>>>>>>>>  ANT BX    = %lf\n", corr1->antenna[31].bx);
  fprintf(stdout, "========= >>>>>>>>>  ANT BY    = %lf\n", corr1->antenna[31].by);
  fprintf(stdout, "========= >>>>>>>>>  ANT BZ    = %lf\n", corr1->antenna[31].bz);

  fprintf(stdout, "\n=========== >>>>>>>>>  SAMP ANT  = %d\n", corr1->sampler[40].ant_id);
  fprintf(stdout, "=========== >>>>>>>>>  SAMP BAND = %d\n", corr1->sampler[40].band);
  fprintf(stdout, "=========== >>>>>>>>>  SAMP FFT  = %d\n", corr1->sampler[40].fft_id);
  fprintf(stdout, "=========== >>>>>>>>>  SAMP DPC  = %d\n", corr1->sampler[40].dpc);

  fprintf(stdout, "\n========= >>>>>>>>>  MACS      = %d\n", corr1->corrpar.macs);
  fprintf(stdout, "========= >>>>>>>>>  CHANNELS  = %d\n", corr1->corrpar.channels);
  fprintf(stdout, "========= >>>>>>>>>  POLS      = %d\n", corr1->corrpar.pols);
  fprintf(stdout, "========= >>>>>>>>>  STA       = %d\n", corr1->corrpar.sta);
  fprintf(stdout, "========= >>>>>>>>>  STATIME   = %d\n", corr1->corrpar.statime);
  fprintf(stdout, "========= >>>>>>>>>  F_STEP    = %lf\n", corr1->corrpar.f_step);
  fprintf(stdout, "========= >>>>>>>>>  CLOCK     = %lf\n", corr1->corrpar.clock);

  fprintf(stdout, "\n========= >>>>>>>>>  ANTMASK   = %x\n", corr1->daspar.antmask);
  fprintf(stdout, "========= >>>>>>>>>  SAMPLERS  = %d\n", corr1->daspar.samplers);
  fprintf(stdout, "========= >>>>>>>>>  BASELINES = %d\n", corr1->daspar.baselines);
  fprintf(stdout, "========= >>>>>>>>>  CHANNELS  = %d\n", corr1->daspar.channels);
  fprintf(stdout, "========= >>>>>>>>>  LTA       = %d\n", corr1->daspar.lta);
  fprintf(stdout, "========= >>>>>>>>>  BANDMASK  = %d\n", corr1->daspar.bandmask);
  fprintf(stdout, "========= >>>>>>>>>  MODE      = %d\n", corr1->daspar.mode);
  fprintf(stdout, "========= >>>>>>>>>  CHAN[127  ]= %d\n", corr1->daspar.chan_num[127]);
  fprintf(stdout, "========= >>>>>>>>>  MJD_REF   = %lf\n", corr1->daspar.mjd_ref);
  fprintf(stdout, "========= >>>>>>>>>  T_UNIT    = %f\n", corr1->daspar.t_unit);
//  fprintf(stdout, "========= >>>>>>>>>  DUT1      = %15.12lf \n", (corr1->daspar.dut1)*(12.0*3600/M_PI));
#endif
  return 0;
}

int get_lta_hdr(FILE *fp,LtaInfo *linfo)
{
  char    line[VALLEN],*buf,*p;
  int     recl,hdr_recs, ahdr_recs,bhdr_recs;
  KeyType keys[MAX_KEYS];
  int     nkey;
  rewind(fp);
  /* get the record length, CorrHdr size etc. */
  if(fgets(line,VALLEN,fp) ==NULL)
  {fprintf(stderr,"GlobalHdrFmtErr\n"); return 1;}
  if(strncmp(line,"HDR",3))
  {fprintf(stderr,"No HDR keyword found!\n"); return 1;}
  if(sscanf(line,"%*s %d %d %d\n",&recl,&hdr_recs,&bhdr_recs) !=3)
  {fprintf(stderr,"FmtErr in HDR line\n");}
  rewind(fp);

  if((ahdr_recs=hdr_recs-bhdr_recs)<=0)
  { fprintf(stderr,"Total HdrRecs (%d) <= BinHdrRecs (%d)!\n",
	     hdr_recs,bhdr_recs);
    return 1;
  }

  if((buf=(char*)malloc(recl*(hdr_recs+1)))==NULL)
  { fprintf(stderr,"MallocError\n"); return 1;}

  if(fread(buf,recl,hdr_recs+1,fp)!=hdr_recs+1)
  {fprintf(stderr,"File too short\n"); return 1;}

  nkey=set_lta_keys(&linfo->lhdr,keys);
  if(parse_ascii_hdr(nkey,keys,buf,"*{ Init.def","*} AddCorr") !=nkey)
  { fprintf(stderr,"Too few keywords found\n"); return -1; }
  strcpy(linfo->lhdr.old_version,linfo->lhdr.version);
  if(get_corr_lite(buf,&linfo->lhdr.vinfo.corrlite)) return -1;
  if(get_ant_info(buf,linfo->lhdr.vinfo.antenna)) return -1;
  if(get_samp_info(buf,&linfo->lhdr.vinfo)) return -1;
  if(get_base_info(buf,linfo->lhdr.vinfo.base)) return -1;
  if(get_corrsel(buf,&linfo->lhdr.vinfo)) return -1;

  memcpy(&linfo->corr,buf+recl*ahdr_recs,CorrSize);
//  fprintf(stdout, "========= >>>>>>>>>  %f   %f \n", linfo->corr.daspar.mjd_ref, linfo->corr.daspar.t_unit);
/*
  fprintf(stdout, "CORRSIZE = %d\n", CorrSize);
  fprintf(stdout, "AntennaTypeSize = %d\n", AntennaTypeSize);
  fprintf(stdout, "SamplerTypeSize = %d\n", SamplerTypeSize);
  fprintf(stdout, "BaseParTypeSize = %d\n", BaseParTypeSize);
  fprintf(stdout, "CorrParTypeSize = %d\n", CorrParTypeSize);
  fprintf(stdout, "DasParTypeSize = %d\n\n", DasParTypeSize);
*/
  linfo->keep_version=1;
  memcpy(&linfo->corr,buf+recl*ahdr_recs,CorrSize);
  //
//  fprintf(stdout, "===== >>> VERSION = %s\n", linfo->lhdr.version);
//  if(strstr(linfo->lhdr.version,"COR30x2 ")==NULL && strstr(linfo->lhdr.version,"HST04") ==NULL);
/*
  if(strstr(linfo->lhdr.version,"COR30x2 ")==NULL && strstr(linfo->lhdr.version,"HST03") ==NULL);
  else if(strstr(linfo->lhdr.version,"COR30x2 ")==NULL && strstr(linfo->lhdr.version,"HST03") ==NULL)
    corr2_to_corr_ref(linfo);  // very old correlator corr 30
  else if(strstr(linfo->lhdr.version,"COR30x2 ")!=NULL && strstr(linfo->lhdr.version,"HST03") ==NULL)
    corr3_to_corr_ref(linfo); // new corr post polar corr30x2
  else if(strstr(linfo->lhdr.version,"COR30x2 ")!=NULL && strstr(linfo->lhdr.version,"HST03") !=NULL)
    corr4_to_corr_ref(linfo); // gsb corr cor30x2 HST03
  else if(strstr(linfo->lhdr.version,"COR30x2 ")!=NULL && strstr(linfo->lhdr.version,"HST04") !=NULL)
      ;
*/

//  fprintf(stdout, "===== >>> VERSION = %s\n", linfo->lhdr.version);
//  fprintf(stdout, "VERSION ?? = %s\n", linfo->lhdr.version );

  if(strstr(linfo->lhdr.byte_seq,"Big"))linfo->data_byte_order=BigEndian;
  else linfo->data_byte_order=LittleEndian;

  if(linfo->data_byte_order!=linfo->local_byte_order)flipCorr(&linfo->corr);

  linfo->lrecs=hdr_recs;
  linfo->lbrecs=bhdr_recs;
  linfo->recl=recl;

  p=buf+recl*hdr_recs;
  if(strncmp(p,"SCAN",4))
  {fprintf(stderr,"No SCAN keyword found : GET_LTA_HDR!\n"); return 1;}
  if(sscanf(p,"%*s %d %d\n",&hdr_recs,&bhdr_recs) !=2)
  {fprintf(stderr,"FmtErr in SCAN line\n");}
  linfo->srecs=hdr_recs;
  linfo->sbrecs=bhdr_recs;

  free(buf);
  return 0;
}

int decompose_date(ScanHdr *shdr)
{ const char *mnames[12]={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
		    "Sep","Oct","Nov","Dec"};
  char *p,weekday[VALLEN],month[VALLEN],time[VALLEN]; 
  int i,h,m,day,year;
  float s;

  sscanf(shdr->date_obs,"%s %s %d %s %d\n",weekday,month,&day,
	 time,&year);
  for(i=0;i<strlen(time);i++)if(time[i]==':')time[i]=' ';
  sscanf(time,"%d %d %f\n",&h,&m,&s);
  p=month;for(p=month;*p!='\0';p++)if(!isspace(*p)) break;
  for(i=0;i<12;i++)if(!strncmp(p,mnames[i],3))break;
  if(i==12){fprintf(stderr,"Bad Date %s\n",shdr->date_obs); return -1;}
  shdr->year=year;
  shdr->day=day;
  shdr->month=i;
  shdr->ist=h*3600+m*60+s;
  return 0;
}
int scan1_to_scan2(ScanHdrTab *stab, char *version)
{ long off,off1; 
  SourceParType src;
  SourceParType *src1=&stab->scan.source;
  ScanHdr       *shdr  =&stab->shdr;
  char *p=(char *)&src;
  char *p1=(char *)src1; 
  int t=shdr->year*1000+shdr->month*100+shdr->day; 
  int t0=2001*1000+8*100+15; /* hdr changed on 15/Aug/2001 */

  if(strstr(version,"LTA")==NULL) /* old: have to rely on data_obs*/
  { if(t<t0)
    { memcpy(&src,src1,sizeof(SourceParType));
      off=NAMELEN+4*sizeof(float)+3*sizeof(double);
      off1=NAMELEN+4*sizeof(float)+sizeof(double);
      memcpy(p+off,p1+off1,sizeof(SourceParType)-off);
      memcpy(p+off1,p1+off1,2*sizeof(double)); /* set ra_mean,dec_mean to ra_app,dec_app */
      memcpy(src1,&src,sizeof(SourceParType));
    }
  }

 return 0;
}
int scan2_to_scan1(ScanInfoType *scan, ScanInfoType *scan1)
{ long off,off1; 
  SourceParType *src=&scan->source;
  SourceParType *src1=&scan1->source;
  char *p=(char *)&src;
  char *p1=(char *)src1; 
  
  memcpy(src1,src,sizeof(SourceParType));
  off1=NAMELEN+4*sizeof(float)+3*sizeof(double);
  off=off1+2*sizeof(double);
  memcpy(p1+off1,p+off,sizeof(SourceParType)-off);

 return 0;
}
int get_scan_hdr(FILE *fp, LtaInfo *linfo,int nscan)
{ char  *buf;
  int   hdr_recs=linfo->srecs,bhdr_recs=linfo->sbrecs,ahdr_recs;
  int   recl=linfo->recl;
  int   nkeys,err;
//	char  temp[5];
  KeyType keys[MAX_KEYS];
  XScanPar tmp;
  ScanHdr *shdr=&linfo->stab[nscan].shdr;
  
  tmp.net_sign[0]='\0';  

  if((buf=(char*)malloc(recl*hdr_recs))==NULL)
  { fprintf(stderr,"MallocError\n"); return 1;}

  if(fread(buf,recl,hdr_recs,fp)!=hdr_recs)
  {fprintf(stderr,"File too short\n"); return 1;}


  if(strncmp(buf,"SCAN",4))
  {fprintf(stderr,"No SCAN keyword found :GET_SCAN_HDR\n"); return 1;}
  nkeys=set_scan_keys(shdr,keys,&tmp);
  if((err=parse_ascii_hdr(nkeys,keys,buf,NULL,NULL)) !=nkeys)
  { if(err !=nkeys-1 || strlen(tmp.net_sign))
    { fprintf(stderr,"Too few keywords found\n"); return 1; }
     else strcpy(tmp.net_sign,"1 1 1 1 "); /*EarlyLTAVer:No net_sign in hdr */
  }

  sscanf(tmp.rf,"%lf %lf\n",&shdr->rf[0],&shdr->rf[1]);
  sscanf(tmp.first_lo,"%lf %lf\n",&shdr->first_lo[0],&shdr->first_lo[1]);
  sscanf(tmp.bb_lo,"%lf %lf\n",&shdr->bb_lo[0],&shdr->bb_lo[1]);
  sscanf(tmp.net_sign,"%d %d %d %d\n",&shdr->net_sign[0],&shdr->net_sign[1],
	    &shdr->net_sign[2],&shdr->net_sign[3]);
  if(decompose_date(shdr)) return -1;

  ahdr_recs=hdr_recs-bhdr_recs;
  memcpy(&linfo->stab[nscan].scan,buf+recl*ahdr_recs,sizeof(ScanInfoType));
  if(scan1_to_scan2(&linfo->stab[nscan],linfo->lhdr.version)) return -1;
  if(linfo->data_byte_order!=linfo->local_byte_order)
    flipScanInfo(&linfo->stab[nscan].scan);


  free(buf);
  return 0;
} 
long locate(char *needle, int nl, char *haystack, int hl)
{ int i,j,k;

 for(i=0;i<hl;i++)
 { for(j=0,k=i;j<nl;j++)
   { if(needle[j]!=haystack[k]) break;
     else { k++; if(k==hl) break;}
   }
   if(j==nl) return i;
 }
 
 return -1;
}
int do_avspc(char *dbuf, LtaHdr *lhdr)
{ int   i,j;
  VisInfo *vinfo=&lhdr->vinfo;
  int   doff=lhdr->data_off,channels=lhdr->vinfo.channels;
  float rave,iave,*iw;
  for(i=0;i<vinfo->baselines;i++)
  { iw=(float*)(dbuf+doff+i*channels*2*sizeof(float));
    for(rave=iave=0.0,j=channels/6;j<5*channels/6;j++)
    { rave += iw[2*j]/(2.0/3.0*channels);
      iave += iw[2*j+1]/(2.0/3.0*channels);
    }
    iw[0]=rave;  iw[1]=iave;
  }
  return 0;
}
int do_bcal(char *dbuf, char *calbuf, LtaHdr *lhdr)
{ float *iw,*cw,dr,di,cr,ci,ca;
  int   i;
  int   doff=lhdr->data_off,nd=lhdr->data_size/sizeof(float);

  iw = (float *)(dbuf+doff);
  cw = (float *)(calbuf+doff);
  for(i=0;i<nd;i+=2)
  { cr=cw[i]; ci=cw[i+1]; ca=cr*cr+ci*ci;
    dr=iw[i]; di=iw[i+1];
    iw[i] = (dr*cr+di*ci)/ca;
    iw[i+1] = (di*cr-dr*ci)/ca;
  }
  return 0;
} 
int do_median(char *obuf, char *ibuf,int recs, LtaHdr *lhdr)
{ float *x,*ow,med;
  int   i,j;
  int   doff=lhdr->data_off,nd=lhdr->data_size/sizeof(float);
  int   recl=lhdr->recl;
  if(recs<2){fprintf(stderr,"Need at least 2 records to average\n"); return 1;}
  if((x=(float*)malloc(recs*sizeof(float)))==NULL)
  { fprintf(stderr,"Malloc error\n"); return 1;}
  
  ow=(float*)(obuf+doff);
  for(i=0;i<nd;i++)
  { for(j=0;j<recs;j++) 
      x[j]=*(float*)(ibuf+(j*recl)+doff+i*sizeof(float));
    sort(recs,x-1);
    if(recs%2) med=x[recs/2];
    else med = 0.5*(x[recs/2]+x[recs/2-1]);
    ow[i]=med;
  }
  free(x);
  return 0;
}
int mk_median(int change_endian,int scan, char *mbuf,LtaInfo *linfo, FILE *fp)
{  char *dbuf;
   int i, recs=linfo->stab[scan].recs,recl=linfo->recl;

   if((dbuf=(char*)malloc(recl*recs))==NULL)
   { fprintf(stderr,"Malloc error\n"); return 1;}
   rewind(fp);ltaseek(fp,linfo->stab[scan].start_rec+linfo->srecs,linfo->recl);
   for(i=0;i<recs;i++)
   { fread(dbuf+(i*recl),recl,1,fp);
     if(change_endian)flipData(dbuf+(i*recl),&linfo->lhdr);
   }
   if(do_median(mbuf,dbuf,recs,&linfo->lhdr)) return 1;
   free(dbuf);
   return 0;
}
int do_average(char *ob, char *ib, LtaHdr *lhdr)
{ float *iw,*ow,*iwt,*owt;
  int   i;
  int doff=lhdr->data_off,woff=lhdr->wt_off,nd=lhdr->data_size/sizeof(float);
  iw = (float *)(ib+doff);
  ow = (float *)(ob+doff);
  iwt = (float *)(ib+woff);
  owt = (float *)(ob+woff);
  for(i=0;i<nd;i++)
    ow[i] = (ow[i]*(*owt)+ iw[i]*(*iwt))/(*owt+*iwt);

  *owt += *(iwt);
  return 0;
}
void flipData(char *dbuf, LtaHdr *lh)
{ double t;
  swap_long(dbuf+lh->data_off,lh->data_size/sizeof(float));
  swap_long(dbuf+lh->wt_off,1);
  memcpy(&t,dbuf+lh->time_off,sizeof(double));
  swap_d(&t,1);
  memcpy(dbuf+lh->time_off,&t,sizeof(double));
} 
    
int ltaseek(FILE *fp, int nrec, int recl)
{ int maxrec=(1<<30)/recl;
  int i,k;
  
  k=nrec>0?1:-1;
  if(k*nrec<maxrec)
  { fseek(fp,nrec*recl,SEEK_CUR); return 0;}

  while(k*nrec>0)
  { i=(k*nrec>maxrec)?k*maxrec:nrec;
    fseek(fp,i*recl,SEEK_CUR);
    nrec -= i; 
  }
  return 0;
}
int make_scantab(FILE *fp, LtaInfo *linfo)
{ char *dbuf,*p,seq[16];
  int   i,nrec=0,off,scanno,recl=linfo->recl; 
  ScanHdrTab *stab =linfo->stab;

  if((dbuf=(char*)malloc(recl))==NULL)
  { fprintf(stderr,"Failed to malloc %d bytes\n",recl); return -1;}
  
  fseek(fp,-recl,SEEK_END);
  while(nrec<5)
  { fread(dbuf,1,recl,fp);
    if(!strncmp(dbuf,"DATA",4))break;
    dbuf[recl-1]='\0';
    for(i=0;i<recl;i++)if(dbuf[i]=='\0')dbuf[i]=' ';
    if((p=strstr(dbuf,"DATA"))!=NULL)
    { off=(p-dbuf)-recl; fseek(fp,off-recl,SEEK_CUR); }
    else fseek(fp,-2*recl,SEEK_CUR);
    nrec++;
  }
  if(nrec==5)
  { int zero_scan=0,frac_scan=0;
    fprintf(stderr,"Checking for zero length last scan ");
    rewind(fp);ltaseek(fp,-linfo->srecs,recl);
    fread(dbuf,1,recl,fp);
    if(!strncmp(dbuf,"SCAN",4))
    { fprintf(stderr,"Found...ignoring it\n");
      ltaseek(fp,-2,recl);
      fread(dbuf,1,recl,fp);
      if(!strncmp(dbuf,"DATA",4))zero_scan=1;
    }
    if(!zero_scan)
    { fprintf(stderr,"Failed...Checking for fractional length last scan\n");
      fseek(fp,-recl,SEEK_END);
      fread(dbuf,1,recl,fp);
      dbuf[recl]='\0';for(i=0;i<recl;i++)if(dbuf[i]=='\0')dbuf[i]=' ';
      if((p=strstr(dbuf,"DATA"))!=NULL)
      { strncpy(seq,p+4,10); seq[10]='\0';seq[4]=' ';
        sscanf(seq,"%d %d\n",&scanno,&nrec);
	if(nrec==0)
        { fprintf(stderr,"Found fractional length last scan (ignoring it)\n");
  	  if(scanno==0)
	  { fprintf(stderr,"(File has only one scan with <1 record\n");
	    return -1;
	  }
	  else
          { off=(p-dbuf)-recl; fseek(fp,off-recl,SEEK_CUR);
            ltaseek(fp,-linfo->srecs,recl);
	    fread(dbuf,1,recl,fp);
	    if(!strncmp(dbuf,"DATA",4))frac_scan=1;
	  }
	}
      }
    }
    if(!zero_scan && !frac_scan);
    { fprintf(stderr,"Missing Signature: Can't make ScanTab\n"); return -1;}
  }
  strncpy(seq,dbuf+4,10); seq[10]='\0';seq[4]=' ';
  sscanf(seq,"%d %d\n",&scanno,&nrec);nrec++;
  linfo->scans=scanno+1;
  if((stab=linfo->stab=(ScanHdrTab *)malloc(sizeof(ScanHdrTab)*linfo->scans))==NULL)
  { fprintf(stderr,"Malloc Failure\n"); return -1;}
  for(i=scanno;i>=0;i--)
  { stab[i].recs=nrec;
    ltaseek(fp,-(linfo->srecs+nrec),recl);
//	  fprintf(stdout, "=====>>>>>> REACHED HERE\n");
    if(get_scan_hdr(fp,linfo,i)) return -1;
    ltaseek(fp,-(linfo->srecs+1),recl);
    if(i==0)break;
    fread(dbuf,1,recl,fp);
    if(strncmp(dbuf,"DATA",4))
    { fprintf(stderr,"Checking for zero length scan\n");
      ltaseek(fp,-linfo->srecs,recl);
      fread(dbuf,1,recl,fp);
      if(!strncmp(dbuf,"SCAN",4))
      { nrec=0; ltaseek(fp,(linfo->srecs-1),recl); continue;}
      fprintf(stderr,"Sync Lost: Can't make ScanTab\n"); return -1;
    }
    strncpy(seq,dbuf+4,10); seq[10]='\0';seq[4]=' ';
    sscanf(seq,"%d %d\n",&scanno,&nrec);nrec++;
    if(scanno!=i-1)
    {fprintf(stderr,"Scans Out of Seq: Can't make ScanTab\n"); return -1;}
  }
  
  stab[0].start_rec=linfo->lrecs;
  for(i=1;i<linfo->scans;i++)
    stab[i].start_rec=stab[i-1].start_rec+stab[i-1].recs+linfo->srecs;

  return 0;
}
int print_corr_lite(CorrLite *corrlite)
{ printf("%s\n",StartCorrsys);
  printf("VERSION =%s\n",corrlite->version);
  printf("CLOCK   =%lf %s\n",corrlite->clock,corrlite->clock_unit);
  printf("STA     =%d %s\n",corrlite->sta,corrlite->sta_unit);
  printf("STATIME =%d %s\n",corrlite->statime,corrlite->statime_unit);
  printf("T_UNIT  =%lf %s\n",corrlite->t_unit,corrlite->t_unit_unit);
  printf("ANTS    =%d\n",corrlite->ants);
  printf("F_ENABLE=%d\n",corrlite->f_enable);
  printf("X_ENABLE=%d\n",corrlite->x_enable);
  printf("CHAN_MAX=%d\n",corrlite->chan_max);
  printf("POLS    =%d\n",corrlite->pols);
  printf("%s\n",EndCorrsys);
  return 0;
}
int print_ant_info(VisInfo *vinfo)
{ AntennaParType *ant=&vinfo->antenna[0]; 
  int   i;

  printf("%s\n",StartAntenna);
  for(i=0;i<vinfo->antennas;i++)
  { printf("ANT%02d = %4s %10.2f %10.2f %10.2f %10.2f %10.2f\n",i,ant->name,ant->bx,
	   ant->by,ant->bz,ant->d0[0],ant->d0[1]);
    ant++;
  }
  printf("%s\n",EndAntenna);

  return 0;
}
int print_samp_info(VisInfo *vinfo)
{ SampInfo *samp=&vinfo->samp[0]; 
  int   i;
  printf("%s\n",StartSampler);
  for(i=0;i<vinfo->samplers;i++)
  { if(strlen(samp->ant))
       printf("SMP%02d = %4s %8s %03d\n",i,samp->ant,samp->band,samp->fft_id);
    else 
       printf("SMP%02d = None\n",i);
    samp++;
  }
  printf("%s\n",EndSampler);
  return 0;
}
int print_base_info(VisInfo *vinfo)
{ BaselineType *base=&vinfo->base[0]; 
  int   i;
  printf("%s\n",StartBaseline);
  for(i=0;i<vinfo->baselines;i++)
  { printf("BAS%03d = %2d %2d %2d %2d %2d %2d %4s %8s %4s %8s\n",i,base->ant[0],
	   base->band[0],base->ant[1],base->band[1],base->samp[0],
	   base->samp[1],base->antname[0],base->bandname[0],
	   base->antname[1],base->bandname[1]);
    base++;
  }
  printf("%s\n",EndBaseline);
  return 0;
}
int print_corrsel(VisInfo *vinfo)
{ char word[VALLEN];
  int print_range(char *str, short *val, int n ) ;

  printf("%s\n",StartCorrsel); 
  printf("MODE    = %d\n",vinfo->mode); 
  printf("LTA     = %d\n",vinfo->lta); 
  print_range(word,vinfo->samp_num,vinfo->samplers);
  printf("SAMP_NUM= %s\n",word);
  print_range(word,vinfo->chan_num,vinfo->channels);
//  fprintf(stdout, "word = %s\n", word);
  printf("CHAN_NUM= %s\n",word);
  printf("%s\n",EndCorrsel); 
  return 0;
}
int print_lta_hdr(LtaHdr *lhdr)
{ printf("*** Begin ASCII Global Header ***\n");
  printf("RECL    =%d\n",lhdr->recl);
  printf("HDR_RECS=%d\n",lhdr->recs);
  printf("REC_FORM=%s\n",lhdr->rec_form);
  printf("OBS_MODE=%s\n",lhdr->obs_mode);
  printf("VERSION =%s\n",lhdr->version);
  printf("D_TYPE  =%s\n",lhdr->d_type);
  printf("D_SIZE  =%s\n",lhdr->d_size);
  printf("D_ALIGN =%s\n",lhdr->d_align);
  printf("INT_REP =%s\n",lhdr->int_rep);
  printf("FL_REP  =%s\n",lhdr->fl_rep);
  printf("BYTE_SEQ=%s\n",lhdr->byte_seq);
  printf("ANTENNAS=%d\n",lhdr->vinfo.antennas);
  printf("SAMPLERS=%d\n",lhdr->vinfo.samplers);
  printf("BASELINE=%d\n",lhdr->vinfo.baselines);
  printf("CHANNELS=%d\n",lhdr->vinfo.channels);
  printf("FLG_OFF =%d\n",lhdr->flg_off);
  printf("FLG_SIZE=%d\n",lhdr->flg_size);
  printf("FLGRECOF=%d\n",lhdr->flg_rec_off);
  printf("FLGRECSZ=%d\n",lhdr->flg_rec_size);
  printf("FLGANTOF=%d\n",lhdr->flg_ant_off);
  printf("FLGANTSZ=%d\n",lhdr->flg_ant_size);
  printf("FLGSMPOF=%d\n",lhdr->flg_smp_off);
  printf("FLGSMPSZ=%d\n",lhdr->flg_smp_size);
  printf("FLGBASOF=%d\n",lhdr->flg_bas_off);
  printf("FLGBASSZ=%d\n",lhdr->flg_bas_size);
  printf("FLGDATOF=%d\n",lhdr->flg_dat_off);
  printf("FLGDATSZ=%d\n",lhdr->flg_dat_size);
  printf("TIME_OFF=%d\n",lhdr->time_off);
  printf("TIMESIZE=%d\n",lhdr->time_size);
  printf("WT_OFF  =%d\n",lhdr->wt_off);
  printf("WT_SIZE =%d\n",lhdr->wt_size);
  printf("PAR_OFF =%d\n",lhdr->par_off);
  printf("PAR_SIZE=%d\n",lhdr->par_size);
  printf("DATA_OFF=%d\n",lhdr->data_off);
  printf("DATASIZE=%d\n",lhdr->data_size);
  printf("DATAFMT =%s\n",lhdr->data_fmt);
  if(print_corr_lite(&lhdr->vinfo.corrlite)) return -1;
  if(print_ant_info(&lhdr->vinfo)) return -1;
  if(print_samp_info(&lhdr->vinfo)) return -1;
  if(print_base_info(&lhdr->vinfo)) return -1;
  if(print_corrsel(&lhdr->vinfo)) return -1;


  printf("*** End ASCII Global Header ***\n");
  return 0;
}
void update_signature(char *dbuf,int scanno, int recno)
{
 sprintf(dbuf,"%4s%04d.%05d  ","DATA",scanno,recno);dbuf[16] = ' ' ;
 return;
}
    
int cmp_global_hdr(LtaInfo *linfo,LtaInfo *linfo1,char *chkstr)
{  char *p;
   CorrType *corr=&linfo->corr,*corr1=&linfo1->corr; 
   DasParType *daspar=&corr->daspar,*daspar1=&corr1->daspar;
   int i,j;
   
   for(p=chkstr;*p!='\0';p++)
   { switch(*p)
     { case 'a':if(daspar->antmask!=daspar1->antmask)
                { fprintf(stderr,"antmask mismatch\n"); return -1;}  
		break;
       case 'b':if(daspar->bandmask!=daspar1->bandmask)
   	        { fprintf(stderr,"bandmask mismatch\n"); return -1;}
                break;
       case 'c':if(daspar->channels!=daspar1->channels)
                { fprintf(stderr,"channel number mismatch\n"); return -1;}
                for(i=0;i<daspar->channels;i++)
		  if(daspar->chan_num[i]!=daspar1->chan_num[i])
		  { fprintf(stderr,"channel %d mismatch\n",i); return -1;}
	        break;
       case 'p':if(daspar->baselines!=daspar1->baselines)   /* products */
                { fprintf(stderr,"baseline number mismatch\n"); return -1;}
	        break;
       case 'B':if(daspar->baselines!=daspar1->baselines) 
                { fprintf(stderr,"baseline number mismatch\n"); return -1;}
                for(i=0;i<daspar->baselines;i++)
		{ for(j=0;j<2;j++)
		  { SamplerType *samp=&corr->baseline[i].samp[j];
                    SamplerType *samp1=&corr1->baseline[i].samp[j];
                    if(samp->ant_id != samp1->ant_id ||
		       samp->band != samp1->band)
                    { fprintf(stderr,"baseline %d mismatch\n",i); return -1;}
		  }
		}
                break;
       case 'S':if(daspar->samplers!=daspar1->samplers)
                { fprintf(stderr,"sampler number mismatch\n"); return -1;}
 	        break;
       case 's':if(daspar->samplers!=daspar1->samplers)
                { fprintf(stderr,"sampler number mismatch\n"); return -1;}
                for(i=0;i<daspar->samplers;i++)
		  if(corr->sampler[i].ant_id != corr1->sampler[i].ant_id||
		     corr->sampler[i].band != corr1->sampler[i].band)
		  { fprintf(stderr,"sampler %d mismatch\n",i); return -1;}
                break;
        default:fprintf(stderr,"Illegal Check Type %c\n",*p); return -1;
     }
   }
   
   return 0;
}


