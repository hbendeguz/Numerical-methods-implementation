#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double GravForce(double x,double y, double k){
    double F=-(x*k)/(pow((pow(x,2)+pow(y,2)) ,1.5));
    return F;
}
double KineticEnergy(double vx, double vy, double mass){
   double KE = mass/2*(pow(vx,2)+pow(vy,2));
    return KE;
}
double GravPot(double x, double y, double k, double mass){
    double GP=-k*mass/(pow((pow(x,2)+pow(y,2)),0.5));
    return GP;
}

int main(int argc, const char * argv[]) {
    // Constants...
    double G=(0.667384E-10); //[Si]
    double M=5.972E24;
    double k=G*M; //
   // double k2 = G*
    double Apogeum=1.52097701e11;//[m]
    double VelocityAtApogeum= 29291; // m/s
    double Apogeum2=0.69817079e11;//[m]
    double VelocityAtApogeum2= 38860; // m/s
    double y_2 = 0;
    double y = 0;
    double x = Apogeum;
    double x_2 = Apogeum2;
    double vx=0;
    double vy=VelocityAtApogeum;
    double vy_2=VelocityAtApogeum2;
    
    double Mercury = 0.35E23;
    double Sun = 1.989E30;
    // Timeset.. p periods
    double T=24*28*3600;
    double p=20;
    double dt=T/1700;
    double t=0;
    //x and v temporary variables...
    double x1=0;
    double y1=0;
    double vx1=0;
    double vy1=0;
    
    double x2=0;
    double y2=0;
    double vx2=0;
    double vy2=0;
    double vx_2 = 0;
    
    
    double E;
    //Files..
    FILE *g;
    g=fopen("orbit.txt","w");
   // FILE *f;
   // f=fopen("E(time).txt", "w");
    //Euler-method...
    while(t<T*p){
        x1=x+dt*vx;
        y1=y+dt*vy;
        vx1=vx+dt*GravForce(x, y, k)/M*Sun + dt*GravForce(x_2-x, y_2-y, k)/M*Mercury;
        vy1=vy+dt*GravForce(y, x, k)/M*Sun + dt*GravForce(y_2-y, x_2-x, k)/M*Mercury;
        
        x2=x_2+dt*vx_2;
        y2=y_2+dt*vy_2;
        vx2=vx_2+dt*GravForce(x_2, y_2, k)/M*Sun + dt*GravForce(x-x_2, y-y_2, k);
        vy2=vy_2+dt*GravForce(y_2, x_2, k)/M*Sun + dt*GravForce(y-y_2, x-x_2, k);
        
        
        
       // E=KineticEnergy(vx, vy, Sun)+GravPot(x, y, k, Sun);
        fprintf(g,"%lf %lf %lf %lf \n" ,x1,y1, x2, y2 );
        //fprintf(f,"%lf %lf \n" ,t,E);
        t=t+dt;
        x=x1;
        y=y1;
        vx=vx1;
        vy=vy1;
        x_2=x2;
        y_2=y2;
        vx_2=vx2;
        vy_2=vy2;
        
    }
    fclose(g);
    //fclose(f);
    return 0;
}
