
/*******************************************************************************
*
* File mutils.h
*
* Copyright (C) 2005, 2007, 2008, 2011, 2012, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Utility functions used in main programs.
*
* The externally accessible functions are
*
*   int find_opt(int argc,char *argv[],char *opt)
*     On process 0, this program compares the string opt with the arguments
*     argv[1],..,argv[argc-1] and returns the position of the first argument
*     that matches the string. If there is no matching argument, or if the
*     program is called from another process, the return value is 0.
*
*   int digits(double x,double dx,char *fmt)
*     Assuming x is a value with error dx, this program returns the number n
*     of fractional digits to print so that all significant digits plus two
*     more are shown. The print format fmt has to be "e" or "f" depending on
*     whether the number is to be printed using the "%.ne" or "%.nf" format
*     string. In the second case dx has to be in the range 0<dx<1, and
*     (int)(10^n*dx) is then a two-digit integer that represents the error
*     in the last two digits of the printed value.
*
*   int fdigits(double x)
*     Returns the smallest integer n such that the value of x printed with
*     print format %.nf coincides with x up to a relative error at most a
*     few times the machine precision DBL_EPSILON.
*
*   void check_dir(char* dir)
*     This program checks whether the directory dir is accessible. An
*     error occurs if it is not.
*
*   int name_size(char *format,...)
*     On process 0, this program returns the length of the string that
*     would be printed by calling sprintf(*,format,...). The format
*     string can be any combination of literal text and the conversion
*     specifiers %s, %d and %.nf (where n is a positive integer). When
*     called on other processes, the program does nothing and returns
*     the value of NAME_SIZE.
*
*   long find_section(FILE *stream,char *title)
*     This program scans stream for a line starting with the string "[title]"
*     (after any number of blanks). It terminates with an error message if no
*     such line is found or if there are several of them. The program returns
*     the offset of the line from the beginning of the file and positions the
*     file pointer to the next line.
*
*   long read_line(FILE *stream,char *tag,char *format,...)
*     This program scans stream and reads a line of text in a controlled
*     manner, as described in the notes below. The tag can be the empty
*     string "" and must otherwise be an alpha-numeric word that starts
*     with a letter. If it is not empty, the program searches for the tag
*     in the current section. An error occurs if the tag is not found. The
*     program returns the offset of the line from the beginning of the file
*     and positions the file pointer to the next line.
*
*   int count_tokens(FILE *stream,char *tag)
*     This program finds and reads a line from stream, exactly as read_line()
*     does, and returns the number of tokens found on that line after the tag.
*     Tokens are separated by white space (blanks, tabs or newline characters)
*     and comments (text beginning with #) are ignored. On exit, the file
*     pointer is positioned at the next line.
*
*   void read_iprms(FILE *stream,char *tag,int n,int *iprms)
*     This program finds and reads a line from stream, exactly as read_line()
*     does, reads n integer values from that line after the tag and assigns
*     them to the elements of the array iprms. An error occurs if less than
*     n values are found on the line. The values must be separated by white
*     space (blanks, tabs or newline characters). On exit, the file pointer
*     is positioned at the next line.
*
*   void read_dprms(FILE *stream,char *tag,int n,double *dprms)
*     This program finds and reads a line from stream, exactly as read_line()
*     does, reads n double values from that line after the tag and assigns
*     them to the elements of the array iprms. An error occurs if less than
*     n values are found on the line. The values must be separated by white
*     space (blanks, tabs or newline characters). On exit, the file pointer
*     is positioned at the next line.
*
* Notes:
*
* The programs find_section() and read_line() serve to read structured
* input parameter files (such as the *.in files in the directory main).
*
* Parameter lines that can be read by read_line() must be of the form
*
*   tag v1 v2 ...
*
* where v1,v2,... are data values (strings, integers or floating-point
* numbers) separated by blanks. If the tag is empty, the first data value
* may not be a string. Such lines are read by calling
*
*   read_line(tag,format,&var1,&var2,...)
*
* where var1,var2,... are the variables to which the values v1,v2,... are
* to be assigned. The format string must include the associated sequence
* of conversion specifiers %s, %d, %f or %lf without any modifiers. Other
* tokens are not allowed in the format string, except for additional blanks
* and a newline character at the end of the string (none of these have any
* effect).
*
* The programs find_section() and read_line() ignore blank lines and any text
* appearing after the character #. Lines longer than NAME_SIZE-1 characters are
* not permitted. Each section may occur at most once and, within each section,
* a line tag may not appear more than once. The number of characters written
* to the target string variables is at most NAME_SIZE-1. Buffer overflows are
* thus excluded if the target strings are of size NAME_SIZE or larger.
*
*******************************************************************************/
#define mutils_H


#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <iterator>
#include <sstream>

#include "mutils.hpp"

#define NAME_SIZE 1000

#if ((DBL_MANT_DIG!=53)||(DBL_MIN_EXP!=-1021)||(DBL_MAX_EXP!=1024))
#error : Machine is not compliant with the IEEE-754 standard
#endif

static char line[NAME_SIZE+1];
static char inum[3*sizeof(int)+4];
void gdb_hook(){
    
}


extern void error_root(int test,int no,char *name,char *format,...)
{
   va_list args;
   
   if (test!=0)
   {
      printf("\nError in %s:\n",name);
      va_start(args,format);
      vprintf(format,args);      
      va_end(args);
      printf("\nProgram aborted\n\n");
      exit(no);
   }
}
void error(int test,int no, const char *name,const char *format,...)
{
   va_list args;
   
   if (test!=0)
   {
      printf("\nError in %s:\n",name);
      va_start(args,format);
      vprintf(format,args);      
      va_end(args);
      printf("\nProgram aborted\n\n");
      exit(no);
   }
}

int find_opt(int argc,char *argv[],char *opt)
{
   int k;

   for (k=1;k<argc;k++)
   {
      if (strcmp(argv[k],opt)==0)
         return k;
   }

   return 0;
}


int digits(double x,double dx,char *fmt)
{
   error(dx<0.0,1, "digits [mutils.h]",  "Improper input data (negative error)");

   if (strcmp(fmt,"e")==0)
   {
      if (dx==0.0)
         return DBL_DIG;
      else if (dx>=fabs(x))
         return 1;
      else
         return (int)(floor(1.0+log10(fabs(x)))-floor(log10(dx)));
   }
   else if (strcmp(fmt,"f")==0)
   {
      error((dx==0.0)||(dx>=1.0),1,  "digits [mutils.h]",
             "Improper input data (error out of range for fixed format)");

      return (int)(1.0-floor(log10(dx)));
   }
   else
      error(1,1,"digits [mutils.h]","Unknown data format");

   return 0;
}


int fdigits(double x)
{
   int m,n,ne,k;
   double y,z;

   if (x==0.0)
      return 0;

   y=fabs(x);
   z=DBL_EPSILON*y;
   m=floor(log10(y+z));
   n=0;
   ne=1;

   for (k=0;k<(DBL_DIG-m);k++)
   {
      z=sqrt((double)(ne))*DBL_EPSILON*y;

      if (((y-floor(y))<=z)||((ceil(y)-y)<=z))
         break;

      y*=10.0;
      ne+=1;
      n+=1;
   }

   return n;
}

/* maybe c++ version
void check_dir(char* dir)
{
   int nc,n;
   char *tmp_file;
   FILE *tmp;

   nc=strlen(dir);
   tmp_file=(char*) malloc((nc+6)*sizeof(*tmp_file));
   error_root(tmp_file==NULL,1,"check_dir [mutils.h]",
              "Unable to allocate name string");
   sprintf(tmp_file,"%s/.tmp",dir);

   n=0;
   tmp=fopen(tmp_file,"rb");

   if (tmp==NULL)
   {
      n=1;
      tmp=fopen(tmp_file,"wb");
   }

   error(tmp==NULL,1,"check_dir [mutils.h]",
         "Unable to access directory %s",dir);

   fclose(tmp);
   if (n==1)
      remove(tmp_file);
   free(tmp_file);
}

****C version
void check_dir(char* dir)
{
   int nc,n;
   char *tmp_file;
   FILE *tmp;

   nc=strlen(dir);
   tmp_file=malloc((nc+6)*sizeof(*tmp_file));
   error_root(tmp_file==NULL,1,"check_dir [mutils.h]",
              "Unable to allocate name string");
   sprintf(tmp_file,"%s/.tmp",dir);

   n=0;
   tmp=fopen(tmp_file,"rb");

   if (tmp==NULL)
   {
      n=1;
      tmp=fopen(tmp_file,"wb");
   }

   error(tmp==NULL,1,"check_dir [mutils.h]",
         "Unable to access directory %s",dir);

   fclose(tmp);
   if (n==1)
      remove(tmp_file);
   free(tmp_file);
}
*/

int name_size(char *format,...)
{
   int nlen,ie,n;
   double dmy;
   char *pp,*pc;
   va_list args;

   va_start(args,format);
   pc=format;
   nlen=strlen(format);
   ie=0;
   n=0;

   for (;;)
   {
      pp=strchr(pc,'%');

      if (pp==NULL)
         break;

      pc=pp+1;

      if (pc[0]=='s')
         nlen+=(strlen(va_arg(args,char*))-2);
      else if (pc[0]=='d')
      {
         sprintf(inum,"%d",va_arg(args,int));
         nlen+=(strlen(inum)-2);
      }
      else if (pc[0]=='.')
      {
         if (sscanf(pc,".%d",&n)!=1)
         {
            ie=1;
            break;
         }

         sprintf(inum,".%df",n);
         pp=strstr(pc,inum);

         if (pp!=pc)
         {
            ie=2;
            break;
         }

         nlen+=(n+1-strlen(inum));
         dmy=va_arg(args,double);
         if (dmy<0.0)
            nlen+=1;
      }
      else
      {
         ie=3;
         break;
      }
   }

   va_end(args);
   error(ie!=0,1,"name_size [mutils.h]",
              "Incorrect format string %s (ie=%d)",format,ie);
   return nlen;
}


static int cmp_text(char *text1,char *text2)
{
   size_t n1,n2;
   char *p1,*p2;

   p1=text1;
   p2=text2;

   while (1)
   {
      p1+=strspn(p1," \t\n");
      p2+=strspn(p2," \t\n");
      n1=strcspn(p1," \t\n");
      n2=strcspn(p2," \t\n");

      if (n1!=n2)
         return 0;
      if (n1==0)
         return 1;
      if (strncmp(p1,p2,n1)!=0)
         return 0;

      p1+=n1;
      p2+=n1;
   }
}


static char *get_line(FILE *stream)
{
   char *s,*c;

   s=fgets(line,NAME_SIZE+1,stream);

   if (s!=NULL)
   {
      error(strlen(line)==NAME_SIZE,1,"get_line [mutils.h]",
                 "Input line is longer than NAME_SIZE-1");

      c=strchr(line,'#');
      if (c!=NULL)
         c[0]='\0';
   }

   return s;
}


long find_section(FILE *stream,char *title)
{
   int ie;
   long ofs,sofs;
   char *s,*pl,*pr;


   rewind(stream);
   sofs=-1L;
   ofs=ftell(stream);
   s=get_line(stream);

   while (s!=NULL)
   {
      pl=strchr(line,'[');
      pr=strchr(line,']');

      if ((pl==(line+strspn(line," \t")))&&(pr>pl))
      {
         pl+=1;
         pr[0]='\0';

         if (cmp_text(pl,title)==1)
         {
            error(sofs>=0L,1,"find_section [mutils.h]",
                       "Section [%s] occurs more than once",title);
            sofs=ofs;
         }
      }

      ofs=ftell(stream);
      s=get_line(stream);
   }

   error(sofs==-1L,1,"find_section [mutils.h]",
              "Section [%s] not found",title);
   ie=fseek(stream,sofs,SEEK_SET);
   error(ie!=0,1,"find_section [mutils.h]",
              "Unable to go to section [%s]",title);
   get_line(stream);

   return sofs;
}


static void check_tag(char *tag)
{
   if (tag[0]=='\0')
      return;

   error((strspn(tag," 0123456789.")!=0L)||
              (strcspn(tag," \n")!=strlen(tag)),1,
              "check_tag [mutils.h]","Improper tag %s",tag);
}


static long find_tag(FILE *stream,char *tag)
{
   int ie;
   long tofs,lofs,ofs;
   char *s,*pl,*pr;

   ie=0;
   tofs=-1L;
   lofs=ftell(stream);
   rewind(stream);
   ofs=ftell(stream);
   s=get_line(stream);

   while (s!=NULL)
   {
      pl=strchr(line,'[');
      pr=strchr(line,']');

      if ((pl==(line+strspn(line," \t")))&&(pr>pl))
      {
         if (ofs<lofs)
         {
            ie=0;
            tofs=-1L;
         }
         else
            break;
      }
      else
      {
         pl=line+strspn(line," \t");
         pr=pl+strcspn(pl," \t\n");
         pr[0]='\0';

         if (strcmp(pl,tag)==0)
         {
            if (tofs!=-1L)
               ie=1;
            tofs=ofs;
         }
      }

      ofs=ftell(stream);
      s=get_line(stream);
   }

   error(tofs==-1L,1,"find_tag [mutils.h]","Tag %s not found",tag);
   error(ie!=0,1,"find_tag [mutils.h]",
              "Tag %s occurs more than once in the current section",tag);

   ie=fseek(stream,tofs,SEEK_SET);
   error(ie!=0,1,"find_tag [mutils.h]",
              "Unable to go to line with tag %s",tag);

   return tofs;
}


long read_line(FILE *stream,char *tag,char *format,...)
{
   int is,ic;
   long tofs;
   char *pl,*p;
   va_list args;

   check_tag(tag);

   if (tag[0]!='\0')
   {
      tofs=find_tag(stream,tag);
      get_line(stream);
      pl=line+strspn(line," \t");
      pl+=strcspn(pl," \t\n");
   }
   else
   {
      p=format;
      p+=strspn(p," ");
      error(strstr(p,"%s")==p,1,"read_line [mutils.h]",
                 "String data after empty tag");
      tofs=ftell(stream);
      pl=get_line(stream);
   }

   va_start(args,format);

   for (p=format;;)
   {
      p+=strspn(p," ");
      ic=0;
      is=2;

      if ((p[0]=='\0')||(p[0]=='\n'))
         break;
      else if (p==strstr(p,"%s"))
         ic=sscanf(pl,"%s",va_arg(args,char*));
      else if (p==strstr(p,"%d"))
         ic=sscanf(pl,"%d",va_arg(args,int*));
      else if (p==strstr(p,"%f"))
         ic=sscanf(pl,"%f",va_arg(args,float*));
      else if (p==strstr(p,"%lf"))
      {
         is=3;
         ic=sscanf(pl,"%lf",va_arg(args,double*));
      }
      else
         error(1,1,"read_line [mutils.h]",
                    "Incorrect format string %s on line with tag %s",
                    format,tag);

      error(ic!=1,1,"read_line [mutils.h]",
                 "Missing data item(s) on line with tag %s",tag);

      p+=is;
      pl+=strspn(pl," \t");
      pl+=strcspn(pl," \t\n");
   }

   va_end(args);

   return tofs;
}


long myscanf(int n, char *format,...)
{
   int is,ic=0,i;
   long tofs;
   char *pl,tmp[100];
   char *p;
   va_list args;

   while(ic!=n){
   p=format;
   p+=strspn(p," ");
   
   va_start(args,format);
   for (p=format;;)
   {
      p+=strspn(p," ");
      is=2;

      if ((p[0]=='\0')||(p[0]=='\n'))
         break;
      else if (p==strstr(p,"%s"))
         ic+=scanf("%s",va_arg(args,char*));
      else if (p==strstr(p,"%d"))
         ic+=scanf("%d",va_arg(args,int*));
      else if (p==strstr(p,"%f"))
         ic+=scanf("%f",va_arg(args,float*));
      else if (p==strstr(p,"%lf"))
      {
         is=3; 
         ic+=scanf("%lf",va_arg(args,double*));
      }
      else
         error(1,1,"myscanf [mutils.h]",
                    "Incorrect format string %s on the line to scan",
                    format);
     
        
      p+=is;
   }
   int si=0;
   if (ic!=n){//clean the line
       for (i=0;i<n-ic;i++)
           si+=scanf("%s",tmp);
       ic=0;
   printf("the format required have to be %s try again\n",format);
   }
   
        
   }
   va_end(args);
   tofs=ic;
   return tofs;
}

int count_tokens(FILE *stream,char *tag)
{
   int n;
   char *s;


   check_tag(tag);

   if (tag[0]!='\0')
   {
      find_tag(stream,tag);
      s=get_line(stream);
      s+=strspn(s," \t");
      s+=strcspn(s," \t\n");
   }
   else
      s=get_line(stream);

   s+=strspn(s," \t\n");
   n=0;

   while (s[0]!='\0')
   {
      n+=1;
      s+=strcspn(s," \t\n");
      s+=strspn(s," \t\n");
   }

   return n;
}


void read_iprms(FILE *stream,char *tag,int n,int *iprms)
{
   int nc,ic,i;
   char *s;

   check_tag(tag);

   if (tag[0]!='\0')
   {
      find_tag(stream,tag);
      s=get_line(stream);
      s+=strspn(s," \t");
      s+=strcspn(s," \t\n");
   }
   else
      s=get_line(stream);

   s+=strspn(s," \t\n");
   nc=0;

   while ((s[0]!='\0')&&(nc<n))
   {
      ic=sscanf(s,"%d",&i);

      if (ic==1)
      {
         iprms[nc]=i;
         nc+=1;
         s+=strcspn(s," \t\n");
         s+=strspn(s," \t\n");
      }
      else
         break;
   }

   error(nc!=n,1,"read_iprms [mutils.h]","Incorrect read count");
}


void read_dprms(FILE *stream,char *tag,int n,double *dprms)
{
   int nc,ic;
   double d;
   char *s;

   check_tag(tag);

   if (tag[0]!='\0')
   {
      find_tag(stream,tag);
      s=get_line(stream);
      s+=strspn(s," \t");
      s+=strcspn(s," \t\n");
   }
   else
      s=get_line(stream);

   s+=strspn(s," \t\n");
   nc=0;

   while ((s[0]!='\0')&&(nc<n))
   {
      ic=sscanf(s,"%lf",&d);

      if (ic==1)
      {
         dprms[nc]=d;
         nc+=1;
         s+=strcspn(s," \t\n");
         s+=strspn(s," \t\n");
      }
      else
         break;
   }

   error(nc!=n,1,"read_dprms [mutils.h]","Incorrect read count");
}

void go_to_line(FILE *stream,int line)
{
    char c;
    int i=0,r=1;
    fseek(stream,0,SEEK_SET);

    if (line!=0){
        for (c = getc(stream); c != EOF; c = getc(stream)){
            if (c == '\n') // Increment count if this character is newline   
                    i = i + 1;
            if (i==line){ 
                    r=0;
                    break;
            }
        }
        error(r==1,0,"go_to_line","file does not have enought lines");
    }
    
       
}

void move_line(FILE *stream,int line)
{
    char c;
    int i=0;
    fseek(stream,0,SEEK_CUR);

    if (line!=0){
        for (c = getc(stream); c != EOF; c = getc(stream)){
        if (c == '\n') // Increment count if this character is newline   
                i = i + 1;
        if (i==line) 
                break;
        }
    }
    
}


FILE *open_file(const char * name, const char * option){
//FILE *open_file(const char *option,char *format,...){
       FILE *f =NULL;
       
       f=fopen(name,option);
       error(f==NULL,1,"open_file ",  "Unable to open %s file",name);
       return f;
}

void mysprintf(char *str, size_t size, const char *format, ...){
    int i;
    
    char *str_internal;
    str_internal=(char*) malloc(sizeof(char)*size);
    va_list args;
    va_start(args, format);
    i=vsnprintf (str_internal,size,format, args);
    //perror (str);
    va_end (args);
    error(i>=size,0,"mysprintf","string size=%d larger then the size of char=%d ",i,size);
    snprintf(str,size,"%s",str_internal);
    free(str_internal);
}


template <typename Out>
void split(const std::string &s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    std::vector<std::string>  r;
    std::string empty=" ";
    for (std::string s: elems)
        if (s.empty()==0){
            r.emplace_back(s);
        }
            
    return r;
}

template<typename T>
T convert(std::string& in) {
    T t;
    return t;
}

template<>
int convert<int>(std::string& in) {
    return std::stoi(in);
}
template<>
double convert<double>(std::string& in) {
    return std::stod(in);
}
template<>
std::string convert<std::string>(std::string& in) {
    return in;
}

double timestamp()
{
    struct timeval tm;
    gettimeofday(&tm, NULL);
    return tm.tv_sec + 1.0e-6 * tm.tv_usec;
}
