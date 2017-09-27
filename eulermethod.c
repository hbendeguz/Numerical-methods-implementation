//
//  main.c
//  euler
//
//  Created by bendeguz on 2017. 04. 30..
//  Copyright © 2017. bendeguz. All rights reserved.
//


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
    double Apogeum=405500E3;//[m]
    double VelocityAtApogeum= 964; // m/s
    double y = 0;
    double x= Apogeum;
    double vx=0;
    double vy=VelocityAtApogeum;
    double MoonMass = 0.7349E23;
    // Timeset.. p periods
    double T=24*28*3600;
    double p=700;
    double dt=T/1000;
    double t=0;
    //x and v temporary variables...
    double x1=0;
    double y1=0;
    double vx1=0;
    double vy1=0;
    double E;
    //Files..
    FILE *g;
    g=fopen("orbit.txt","w");
    FILE *f;
    f=fopen("E(time).txt", "w");
    //Euler-method...
    while(t<T*p){
        x1=x+dt*vx;
        y1=y+dt*vy;
        vx1=vx+dt*GravForce(x, y, k);
        vy1=vy+dt*GravForce(y, x, k);
        E=KineticEnergy(vx, vy, MoonMass)+GravPot(x, y, k, MoonMass);
        fprintf(g,"%lf %lf \n" ,x1,y1 );
        fprintf(f,"%lf %lf \n" ,t,E);
        t=t+dt;
        x=x1;
        y=y1;
        vx=vx1;
        vy=vy1;
    }
    fclose(g);
    fclose(f);
    return 0;
}
