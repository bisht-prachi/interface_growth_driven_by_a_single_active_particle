/***Prachi Bisht. Program created on 5/11/2018 ***/

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <random>

#define n 256
#define m 256


std::mt19937 gen{static_cast<long unsigned int>(time(0))};
std::uniform_int_distribution<int> dist_n(0,n-1);
std::uniform_int_distribution<int> dist_m(0,m-1);
std::uniform_real_distribution<double> dist_u(0,1);


void gen_lattice_arr(int *lattice)
{
    for (int i=0; i<n; i=i+2) lattice[i] = 1;  
    for (int i=1; i<n; i=i+2) lattice[i] = -1;

//    for (int i=0; i<n; i=i+1) lattice[i] = 1;
//    int r = 0; int y = 0;
//    do
//    {
//        r = dist_n(gen);
//        if(lattice[r]!=-1)
//            lattice[r]=-1;
//        else continue;

//        y++;
//    }while(y<(n/2));


}



void gen_particle_arr(int *particle)
{
    for (int p=0; p<m; p++)
    	particle[p] = dist_n(gen);
}



int search_sites(int *particle, int r)		//tells no. of particles at site r
{
    int number=0;
    for(int i =0; i<m; i++)
    {
        if (particle[i]==r)
        	number+=1;

    }

    return number;
}



void update_surface (double &beta, int *particle, int *lattice, double *height, int (*search_sites)(int*, int))
{
	       
	int r = dist_n(gen);

	int x = search_sites(particle, (r+1));	       
	if (r==n-1) x = search_sites(particle, 0);

	if (x>=1)
	{
	       if(lattice[r]==1)
	       {

		       if(r!=n-1 && lattice[r+1]==-1)
		       {
				   double p = 1/(1 + exp(-2*beta*x));
				   double s = dist_u(gen);
				   if(s <= p)
				   {
				       lattice[r]*=-1;
				       lattice[r+1]*=-1;
				       height[r+1] -= 2;
				       //current -= 1;

				   }
		       }

		       else if(r==n-1 && lattice[0]==-1)
		       {
				   double p = 1/(1 + exp(-2*beta*x));
				   double s = dist_u(gen);
				   if(s <= p)
				   {
				       lattice[r]*=-1;
				       lattice[0]*=-1;
				       height[0] -= 2;
				       //current -= 1;
				   }

		       }

	       }

	       else if(lattice[r]==-1)
	       {

		       if(r!=n-1 && lattice[r+1]==1)
		       {

				   double p = 1*exp(-2*beta*x)/(1 + exp(-2*beta*x));
				   double s = dist_u(gen);
				   if(s < p)
				   {
				       lattice[r]*=-1;
				       lattice[r+1]*=-1;
				       height[r+1] += 2;
				       //current += 1;

				   }


		       }

		       else if(r==n-1 && lattice[0]==1)
		       {
				   double p = 1*exp(-2*beta*x)/(1 + exp(-2*beta*x));
				   double s = dist_u(gen);
				   if(s < p)
				   {
				       lattice[r]*=-1;
				       lattice[0]*=-1;
				       height[0] += 2;
				       //current += 1;
				   }
		       }
	       }

	}
	///end surface update///
}


void update_particle( int *particle, int *lattice)
{
	       int tag = dist_m(gen);

               int dummy = particle[tag];

               if(dummy!=n-1 & dummy!=0)
               {
                    if(lattice[dummy-1]==1 & lattice[dummy]==-1)
                    {

                        double b = dist_u(gen);

                        if( b>0.5)
                        {
                          particle[tag] =  dummy + 1; //shift[tag] += 1;

                        }

                        else if (b<=0.5)
                        {
                          particle[tag] =  dummy - 1; //shift[tag] -= 1;
                        }


                    }

                    else if(lattice[dummy-1]==-1 & lattice[dummy]==-1)
                    {
                         particle[tag] =  dummy + 1; //shift[tag] += 1;
                    }

                    else if(lattice[dummy-1]==1 & lattice[dummy]==1)
                    {
                         particle[tag] =  dummy - 1; //shift[tag] -= 1;
                    }


               }

               else if(dummy==n-1)
               {
                    if(lattice[dummy-1]==1 & lattice[dummy]==-1)
                    {

                        double b = dist_u(gen);

                        if(b>0.5)
                        {
                          particle[tag] =  0; //shift[tag] += 1;
                        }
                        else if (b<=0.5)
                        {
                          particle[tag] =  dummy - 1; //shift[tag] -= 1;
                        }


                    }

                    else if(lattice[dummy-1]==-1 & lattice[dummy]==-1)
                    {
                         particle[tag] =  0;  //shift[tag] += 1;
                    }

                    else if(lattice[dummy-1]==1 & lattice[dummy]==1)
                    {
                         particle[tag] =  dummy - 1;    //shift[tag] -= 1;
                    }


               }

               else if(dummy==0)
               {
                   if(lattice[n-1]==1 & lattice[dummy]==-1)
                    {

                        double b = dist_u(gen);

                        if(b>0.5)
                        {
                          particle[tag] =  dummy + 1; //shift[tag] += 1;
                        }
                        else if (b<=0.5)
                        {
                          particle[tag] =  n - 1; //shift[tag] -= 1;
                        }


                    }

                    else if(lattice[n-1]==-1 & lattice[dummy]==-1)
                    {
                         particle[tag] =  dummy + 1;    //shift[tag] += 1;
                    }

                    else if(lattice[n-1]==1 & lattice[dummy]==1)
                    {
                         particle[tag] =  n - 1;  //shift[tag] -= 1;
                    }


               }
}


void init_height(double *height, int *lattice)
{
	height[0] = 0;
	for (int i = 1; i<n; i++)
	   height[i] = height[i-1] + lattice[i-1];
}


double width (double *height)
{
	double h_avg = 0.0;
	double w = 0.0;
	
	for (int i = 0; i<n; i++)
	{
		h_avg += height[i];
	}h_avg /= n;
	
	for (int i = 0; i<n; i++)
	{
		w += (height[i] - h_avg)*(height[i] - h_avg);
	}w /= n;
	
	return w;
}


int main()
{
    std::ofstream fout, hout, gout, kout, lout;

    int MC = 50000, ENS = 1000, r =0, dummy = 0, k = 0, current = 0;
    int lattice[n]; int particle[m]; double height[n];

    fout.open("width256w0_b-0.5_flat1.dat", std::ios::out);
    
    double width_ss[MC+1] = {0.0};
    double wi = 0.0;
    
    double beta = -0.5;
    
    
 do{
    
    gen_lattice_arr(lattice);
    gen_particle_arr(particle);
    init_height(height, lattice);
    
    wi = width(height);
    width_ss[0] += wi;
    wi = 0.0;
   
    
    for (int j=1; j<=MC; j++)
    {
	    for (int i=0; i<n; i++)
	    {
		update_surface (beta, particle, lattice, height, search_sites);
	       ///end surface update///

		//update_particle( beta, particle, lattice);
		update_particle( particle, lattice);

	    }///microstep loop

	  wi = width(height);
	  width_ss[j] += wi;
	  wi = 0.0;	   

     }///MC loop  

k++;
}while(k<ENS);


	

for(int j = 0; j<=MC; j++)
	fout<<j<<"\t"<<width_ss[j]/(1.0*ENS)<<"\n";

	

gout.close();
fout.close();


return 0;
}
