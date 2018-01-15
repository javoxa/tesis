
/*--- returns the value of ln[G(xx)]-----*/
// log of the Gamma function

float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


/*
float gammp(float a, float x)
{
//    float a,x;

	void gcf(),gser();
//	void nrerror();
	float gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) fprintf(stderr,"Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}
*/


float gammp(float a, float x) //Returns the incomplete gamma function P (a, x).
{
void gcf(float *gammcf, float a, float x, float *gln);
void gser(float *gamser, float a, float x, float *gln);
//void nrerror(char error_text[]);
float gamser,gammcf,gln;
if (x < 0.0 || a <= 0.0) fprintf(stderr,"Invalid arguments in routine gammp");
if (x < (a+1.0)) { //Use the series representation.
gser(&gamser,a,x,&gln);
return gamser;
} 
else { //Use the continued fraction representation
gcf(&gammcf,a,x,&gln);
return 1.0-gammcf; //and take its complement.
}
}






float gammq(float a, float x) //Returns the incomplete gamma function Q(a, x)=1-P (a, x).
{
void gcf(float *gammcf, float a, float x, float *gln);
void gser(float *gamser, float a, float x, float *gln);
//void nrerror(char error_text[]);
float gamser,gammcf,gln;
if (x < 0.0 || a <= 0.0) fprintf(stderr,"Invalid arguments in routine gammq");
if (x < (a+1.0)) {          //Use the series representation
gser(&gamser,a,x,&gln);
return 1.0-gamser;          //and take its complement.
} 
else {                    //Use the continued fraction representation.
gcf(&gammcf,a,x,&gln);
return gammcf;
}
}




#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(float *gammcf,float a, float x,float *gln)
//float *gammcf,*gln,a,x;
{
//	float gammln();
//	void nrerror();
	int i;
	float an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) fprintf(stderr,"a too large, ITMAX too small in gcf\n");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN


#define ITMAX 100
#define EPS 3.0e-7

void gser(float *gamser,float a,float x,float *gln)
//float *gamser,*gln,a,x;
{
//	float gammln();
//	void nrerror();
	int n;
	float sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) fprintf(stderr,"x less than 0 in routine gser\n");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		fprintf(stderr, "a too large, ITMAX too small in routine gser\n");
		return;
	}
}
#undef ITMAX
#undef EPS
