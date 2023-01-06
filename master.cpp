/***Prachi Bisht. Program created on 28/10/2018 ***/

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <random>

#define n 2048
#define m 1
#define eps 5e-2


std::mt19937 gen{static_cast<long unsigned int>(time(0))};
std::uniform_int_distribution<int> dist_n(0,n-1);
std::uniform_real_distribution<double> dist_u(0,1);


void gen_lattice_arr(int *lattice)
{
	
    for (int i=0; i<n; i=i+2) lattice[i] = 1;  
    for (int i=1; i<n; i=i+2) lattice[i] = -1;

//	int h = n/2;
//    for (int i=0; i<h; i=i+1) lattice[i] = -1;  
//    for (int i=h; i<n; i=i+1) lattice[i] = 1;

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

      particle[0] = n/2;//dist_n(gen);

}



int search_sites(int particle[m], int r)
{
    int number=0;
    for(int i =0; i<m; i++)
    {
        if (particle[i]==r) number+=1;

    }

    return number;
}




int search_particle(int dummy_bond, int particle[m])
{
    if (particle[0]==dummy_bond)
        return 1;
    else
        return 0;

}



void update_surface (int *particle, int *lattice, double *height, int (*search_sites)(int*, int))
{
	       ///begin particle update///
	     int r = dist_n(gen);
	      if(lattice[r]==-1)
		       {

		               if(r!=n-1 && lattice[r+1]==1)
		               {
		               
		               	int x0 = search_sites(particle, (r+1));
				       if (x0>=1)
				       {

		                  
				               lattice[r]*=-1;
				               lattice[r+1]*=-1;  
       			                       height[r+1] += 2;
				       }				   

		               }

		               else if(r==n-1 && lattice[0]==1)
		               {
		               
		               int x0 = search_sites(particle, (0));
				       if (x0>=1)
				       {

		                   
				               lattice[r]*=-1;
				               lattice[0]*=-1;	
       			                       height[0] += 2;
				       }			   

		               }

		       }

		

               ///end surface update///

}


void update_particle( int *particle, int *lattice )
{

               ///begin particle update///
           int chance = dist_n(gen);
           int dummy = particle[0];

            if (chance==dummy)
            {

               if(dummy!=n-1 & dummy!=0)
               {
                    if(lattice[dummy-1]==1 & lattice[dummy]==-1)
                    {

                        double b = dist_u(gen);

                        if( b>0.5)
                        {
                          particle[0] =  dummy + 1;

                        }

                        else if (b<=0.5)
                        {
                          particle[0] =  dummy - 1;
                        }


                    }

                    else if(lattice[dummy-1]==-1 & lattice[dummy]==-1)
                    {
                         particle[0] =  dummy + 1;
                    }

                    else if(lattice[dummy-1]==1 & lattice[dummy]==1)
                    {
                         particle[0] =  dummy - 1;
                    }



               }

               else if(dummy==n-1)
               {
                    if(lattice[dummy-1]==1 & lattice[dummy]==-1)
                    {

                        double b = dist_u(gen);

                        if(b>0.5)
                        {
                          particle[0] =  0;
                        }
                        else if (b<=0.5)
                        {
                          particle[0] =  dummy - 1;
                        }


                    }

                    else if(lattice[dummy-1]==-1 & lattice[dummy]==-1)
                    {
                         particle[0] =  0;
                    }

                    else if(lattice[dummy-1]==1 & lattice[dummy]==1)
                    {
                         particle[0] =  dummy - 1;
                    }



               }

               else if(dummy==0)
               {
                   if(lattice[n-1]==1 & lattice[dummy]==-1)
                    {

                        double b = dist_u(gen);

                        if(b>0.5)
                        {
                          particle[0] =  dummy + 1;
                        }
                        else if (b<=0.5)
                        {
                          particle[0] =  n - 1;
                        }


                    }

                    else if(lattice[n-1]==-1 & lattice[dummy]==-1)
                    {
                         particle[0] =  dummy + 1;
                    }

                    else if(lattice[n-1]==1 & lattice[dummy]==1)
                    {
                         particle[0] =  n - 1;
                    }


               }
            }
}




void init_height(double *height, int *lattice)
{
	height[0] = 0;
	for (int i = 1; i<n; i++)
	   height[i] = height[i-1] + lattice[i-1];
}



int main()
{
    std::ofstream fout, hout, gout, kout, lout, mout, jout, nout, pout;

    int MC = 10000, ENS = 1, r =0, delr = 0, dummy = 0, k = 0;
    double counter = 0.0;
    
//    fout.open("typ_profile2048_t500_1.dat", std::ios::out);
//    gout.open("fwhm256_ew_psu_b-1_.dat", std::ios::out);
//    hout.open("msd2048_ew_psu_b-1_.dat", std::ios::out);

//    kout.open("hmax2048_ew_psu_b-1_.dat", std::ios::out);
//    lout.open("probdx2048_ew_psu_b-1_t500.dat", std::ios::out);
//    mout.open("height2048_t2500.dat", std::ios::out);

//    jout.open("typarea2048_ew_psu_b-1.dat", std::ios::out);
//    jout.open("firstreturn2048.dat", std::ios::out);

//    jout.open("ovisits2048_ew_psu_b-1.dat", std::ios::out);

    nout.open("typflip_stats2048.dat", std::ios::out);

//      pout.open("avg_tilt2048_ew_psu_b-1t8000.dat", std::ios::out);
      
//    int lattice[n] = {0};
//    int particle[m] = {0};
//    double height[n] = {0.0};

    double bias = 0.5;
    int h = n/2;
    
    int timestamp = MC/10 ;

//    double profile_ss[h][timestamp] = {0.0};
    
    
//    double msd_ss[timestamp] = {0.0};


//    double siter_ss[h] = {0.0};


//    double sitea_ss[h] = {0.0};


//    double profile_last[n] = {0.0};

//    double sites_visits_ss[n] = {0.0};

//    double returnstats_ss[MC] = {0.0};

//      double ovisits[timestamp] = {0.0};

//	int avg_tilt[n] = {0};


//	int meanl = 0;
//	int meanr = 0;


//    unsigned int bin_no;
//    int bins = n/2;
//    double probdx[bins] = {0.0};
    
     int flips_ss[MC+1] = {0};
     int flips = 0;

	
//   int lattice[n]; int particle[m];
//   double height[n];
//    
//    gen_lattice_arr(lattice);
//    gen_particle_arr(particle);
//    init_height(height, lattice);
//    
//    int zero = 40000;
//    for (int j = 0; j<zero; j++)
//	{
//	  for (int i=0; i<n; i++)
//          {
//	
//	   update_surface(particle, lattice, height, (*search_sites));
//	   update_particle(particle, lattice);
//	  }	
//	}
    
    
do{

    int lattice[n] = {0};
    int particle[m] = {0};
    double height[n] = {0.0};
//    int sites_visits[n] = {0};

    gen_lattice_arr(lattice);
    gen_particle_arr(particle);
    init_height(height, lattice);  
    
//    double d = 0.0;
//    double y = particle[0];
//    double z = 0.0;

//    delr = 0;
    

//     int siter[n] = {0};
//     int sitea[n] = {0};
//    

//   counter = 0.0;
   flips = 0;
   
//   for (int i=0; i<n; i++)
//   {
//    update_surface (particle, lattice,  height, (*search_sites));
//    update_particle( particle, lattice);
//   }
    
    for (int j=1; j<=MC; j++)
        {     
         

          for (int i=0; i<n; i++)
          {
	       int r = dist_n(gen);	                         
	       
               int x0 = search_sites(particle, (r+1));               
	       if (r==n-1) x0 = search_sites(particle, 0);        
	                                                        
	             	       
               if (x0==1)
               {
               		
//               	sites_visits[r] +=1;
               		
		       if(lattice[r]==-1)
		       {

		               if(r!=n-1 && lattice[r+1]==1)
		               {

				               lattice[r]*=-1;
				               lattice[r+1]*=-1;
       			                       height[r+1] += 2;
       			                       flips += 1;
//				               z = r + 1;
//				   	       d = z - y;
//				   	       y = z;
//				   	       

//				   	       if ( abs(d) > n/2 )
//						{
//						   if ( d > 0 ) d = n - d;
//						   if ( d < 0 ) d = n + d;
//						}
//					    
//					    bin_no = abs(d);
//					    if (bin_no < bins)
//					    probdx[bin_no] +=1.0;

		               }

		               else if(r==n-1 && lattice[0]==1)
		               {
		                
				               lattice[r]*=-1;
				               lattice[0]*=-1;
			                       height[0] += 2;
			                       flips += 1;
//				               z = r + 1;
//				   	       d = z - y;
//				   	       y = z;
//				   	       
//				   	       if ( abs(d) > n/2 )
//						{
//						   if ( d > 0 ) d = n - d;
//						   if ( d < 0 ) d = n + d;
//						}
//					    
//					    bin_no = abs(d);
//					    if (bin_no < bins)
//					    probdx[bin_no] +=1.0;
		               }

		       }
               }


	   int chance = dist_n(gen);
           int dummy = particle[0];

            if (chance==dummy)
            {

               if(dummy!=n-1 & dummy!=0)
               {
                    if(lattice[dummy-1]==1 & lattice[dummy]==-1)
                    {

                        double b = dist_u(gen);

                        if( b>0.5)
                        {
                          particle[0] =  dummy + 1; //delr += 1;

                        }

                        else if (b<=0.5)
                        {
                          particle[0] =  dummy - 1; //delr -= 1;
                        }


                    }

                    else if(lattice[dummy-1]==-1 & lattice[dummy]==-1)
                    {
                         particle[0] =  dummy + 1; //delr += 1;
                    }

                    else if(lattice[dummy-1]==1 & lattice[dummy]==1)
                    {
                         particle[0] =  dummy - 1; //delr -= 1;
                    }


               }

               else if(dummy==n-1)
               {
                    if(lattice[dummy-1]==1 & lattice[dummy]==-1)
                    {

                        double b = dist_u(gen);

                        if(b>0.5)
                        {
                          particle[0] =  0; //delr += 1;
                        }
                        else if (b<=0.5)
                        {
                          particle[0] =  dummy - 1; //delr -= 1;
                        }


                    }

                    else if(lattice[dummy-1]==-1 & lattice[dummy]==-1)
                    {
                         particle[0] =  0;  //delr += 1;
                    }

                    else if(lattice[dummy-1]==1 & lattice[dummy]==1)
                    {
                         particle[0] =  dummy - 1;    //delr -= 1;
                    }


               }

               else if(dummy==0)
               {
                   if(lattice[n-1]==1 & lattice[dummy]==-1)
                    {

                        double b = dist_u(gen);

                        if(b>0.5)
                        {
                          particle[0] =  dummy + 1; //delr += 1;
                        }
                        else if (b<=0.5)
                        {
                          particle[0] =  n - 1; //delr -= 1;
                        }


                    }

                    else if(lattice[n-1]==-1 & lattice[dummy]==-1)
                    {
                         particle[0] =  dummy + 1;    //delr += 1;
                    }

                    else if(lattice[n-1]==1 & lattice[dummy]==1)
                    {
                         particle[0] =  n - 1;  //delr -= 1;
                    }


               }
            }
                                  
//	if (particle[0]==h) counter += 1.0;

     	  }///microstep loop    	  
     	  
     	         
//                 if (j % 10 == 0)
//        	 {
//		       int t = j/10;
//		       ovisits[t] += counter;
////		       counter = 0.0;
//         	 }
            
//            
            
     	  
     	    	  
//     	       	  
//         if (j % 10 == 0)
//         {
//         	int t = j/10;
////         	msd_ss[t-1] += delr*delr;
//         	
//         	for(int i = 0; i<h; i++)
//		   {
//		    	profile_ss[i][t-1] += height[i+h];
//		   }

//         }


	flips_ss[j] += flips;


	
//      	if (particle[0]==h)
//      	{
//      		returnstats_ss[j] += 1.0;
//      		break;
//      	}
      	
      	

      }///MC loop
      
//      for(int i = 0; i<n; i++)
//      {

//	   	sitea_ss[i] += 1.0*sitea[i+h];
//	   	siter_ss[i] += 1.0*siter[i+h];
//	   	profile_last[i] += 1.0*height[i];
//	   	sites_visits_ss[i] += 1.0*sites_visits[i];

//		avg_tilt[i] += lattice[i];
//      }



//	for (int i = h+1; i<n; i=i+2)
//	meanr += lattice[i];
//	for (int i = 0; i<h; i=i+2)
//	meanl += lattice[i];




//	for(int i = 0; i<n; i++)
//        {
//      	fout<<i<<"\t"<<height[i]<<"\n";
//      	}


  k++;
 }while(k<ENS);///history loop
 
// 	std::cout<<MC<<"\t"<<meanl/(1.0*h*ENS)<<"\t"<<meanr/(1.0*h*ENS)<<"\n";
 
 
// 	 for(int i = 0; i<n; i++)
//      		pout<<i<<"\t"<<avg_tilt[i]/(1.0*ENS)<<"\n";
      	
 
 
 
 		for (int j = 0; j<=MC; j++)
			nout<<j<<"\t"<<flips_ss[j]/(1.0*ENS)<<"\n";

 
 
// 		for (int j = 0; j<timestamp; j++)
//			jout<<j*10<<"\t"<<ovisits[j]/(1.0*ENS)<<"\n";


//		for (int j = 1; j<MC; j++)
//			jout<<j<<"\t"<<returnstats_ss[j]/(1.0)<<"\n";

 
		 
//		for (int i = 0; i<bins; i++)
//			lout<<i<<"\t"<<probdx[i]/(1.0*ENS)<<"\n";
 
 
 
 
//	 	for (int j = 0; j<timestamp; j++)
//	 	{
//		 double area = 0.0; 
//		 	for (int i = 0; i<h; i = i+2)
//		 	{
//		 		area +=  (profile_ss[i][j])/(1.0*ENS);
//		 	}
//		 	
//		 	jout<<(1+j)*10<<"\t"<<area<<"\n";
//	 	}
 	
 	
//	 	for (int j = 0; j<timestamp; j++)
//			hout<<(j+1)*10<<"\t"<<msd_ss[j]/(1.0*ENS)<<"\n";
//		


		
		
//		for (int i = 0; i<h; i++)
//		{
//		fout<<i<<"\t";
//	
//			for (int j = 0; j<timestamp; j++)
//			{
//				double s = (profile_ss[i][j])/(1.0*ENS);

//			
//				fout<<s<<"\t";
//			}
//		fout<<"\n";
//		}	

	


	
//		for (int j = 0; j<timestamp; j++)
//		{
//			double s = (profile_ss[0][j])/(1.0*ENS);

//			
//			kout<<10*(j+1)<<"\t"<<s<<"\n";
//		}



		

//		for (int j = 0; j<timestamp; j++)
//		{
//			int site = 0;
//			for (int i = 0; i<h; i = i+2)
//			{
//				double dx = (profile_ss[i][j])/(1.0*ENS);
//			
//				if (abs(dx) <= eps)
//				{
//					site = i;
//					break;
//				}

//			}
//			lout<<10*(j+1)<<"\t"<<site<<"\n";
//			
//		}



		
//		for (int i = 0; i<n; i=i+1)
//		{		

			//mout<<i//<<"\t"<<sitea_ss[i]/(1.0*ENS)
			       //<<"\t"<<siter_ss[i]/(1.0*ENS)
			//fout<<i<<"\t"<<profile_last[i]/(1.0*ENS)
			       //<<"\n";
			       
//			jout<<i<<"\t"<<sites_visits_ss[i]/(1.0*ENS)<<"\n";
//		}
	
	
		
		
		
//	int site = 0;

//	for (int j = 0; j<timestamp; j++)
//	{
//		double hmax = profile_ss[h][j];
//		int i = h;
//		

////		for (int i = h; i>=0; i--)
////		{
////			double delh = profile_ss[i][j] - 0.5*hmax; 
////		
////			if (delh <= 0.0)
////			{
////				site = i;
////				break;
////			}
////			
////			i--;
////		
////		}
////		
//		gout<<(j+1)*100<<"\t"<<(h-site)<<"\n";
//	}
//	
	
	
	

	



gout.close();
fout.close();


return 0;
}
