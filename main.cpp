//
//  main.cpp
//  Thermal SSE
//
//  Created by Aravindh Swaminathan on 13/06/2020.
//  Copyright © 2020 Aravindh Swaminathan. All rights reserved.
//

// SSE with thermal loop updates on XXZ model

#include<iostream>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define min(a,b) (a<b)?a:b
//#include<conio.h>

#define llm 30000 // Array length parameter medium
#define lls 50  // Array length parameter small


using namespace std;


double rnd(int);

double mini(double a, double b)
{
    if(a<b)
       return a;
    else
       return b;
}

struct strtype
{
  int type;
  int bond;

};

struct sites
{
    int x;
    int y;
};


class Thermal
{
   int spin[lls][lls];
   double prob[4];
   int n,L,N,Lx,Ly,Nb,d,NLoops,maxlen,Nbins,MCS;
   /* n is the number of operators
      L is the exapansion Cutoff
      N is the number of atoms
      Nb is the number of bonds
      s is the operator string
      d is the dimension of the lattice
    */
   double beta, h, delta,J,epsilon;
   strtype s[llm];
   sites bsites[2][llm];

 long int X[4*llm], v_first[llm],v_last[llm],vxlist[llm];
   double vxwgt[6];
   double vxprob[4][4][6];
   int vxcode[2][2][2][2];
   int vxnew[4][4][6];
   int vxcodei[4][6]; //tells the kind of spin at leg i of vertex type vx
   int fspn[llm];

   double energy,stamag,sqmag;

   /* Medthod : I want to assign vxcode[1,-1,-1,1]=1, so this will be done as
                vxcode[1][0][1][0] = 1;
       The disallowed ones wll have the value of -1 stored at those memory spots */


public:


   void diagonal_update(); // do the measurements part later, otherwise tested
   //double rnd(int seed);   // done
   void init_thermal();//done
   void gen_linked_list();//done and tested
   int loop_update();//debug
   void calc_probloops();//done and tested, works.
   void makelattice();//done
   int adjust_cutoff(unsigned int loopvar);//done


    Thermal();

    void test();
    void run();
    void test2();



};





Thermal::Thermal()
    {
        
        L=20;
        Lx=4;
        Ly=2;
        N = Lx*Ly;
        n=0;
        d = 2;
        Nb = d * N; // 2 * N
        NLoops=50;
        maxlen=80000;
        beta = 1.0;
        h = 0.0;
        delta = 1.0;
        J = 1.0;
        epsilon = 1.0;
        energy=stamag=sqmag=0;
        Nbins=5;
        MCS=5000;


            /*    for(int i=0;i<llm;i++)
        {
            X[i] = -1;
            v_first[i] = v_last[i] = -1;
            s[i].bond = -1; s[i].type = 0;
            fspn[i] = -1;
            bsites[0][i].x = bsites[0][i].y = bsites[1][i].x = bsites[1][i].y = -1;
        }

        for(int j=0;j<lls;j++)
          for(int k=0;k<lls;k++)
            spin[k][j] = 0; */


        makelattice();
        init_thermal();

        calc_probloops();

}


void Thermal::test2()
{
   /*for(int b=0;b<Nb;b++)
   {
       cout<<"bond number "<<b<<endl;
       cout<<(bsites[0][b].y)*Lx + (bsites[0][b].x)<<" ";
       cout<<(bsites[1][b].y)*Lx + (bsites[1][b].x)<<endl;
   }*/

   int data[100], count[50000];
   int rand_num;
  
   ofstream fout("rantest.txt",ios::out);

   for(int k=0;k<100;k++)
      data[k] = 0;
      
    for(int k=0;k<50000;k++)
       count[k]  = 0;

   for(int j=0;j<50000;j++)
   {
         rand_num = (int)(rnd(0)*100);
              data[rand_num]++;

       }
       
       for(int i=0;i<100;i++)
         count[data[i]]++;
       
   
   
   for(int k=0;k<50000;k++)
      fout<<k<<" "<<count[k]<<endl;
      
    fout.close();

}



void Thermal::run()
    {
        int passed=1,equil=0;
        long int i=0;
        double data[lls];
        double mean=0, stddev=0,meandev=0;

        ofstream fout("nvsiter.txt", ios::out);




            cout<<"The initial state is : "<<endl;

        for(int y=Ly-1;y>=0;y--)
        {
            cout<<endl;
            for(int x=0;x<Lx;x++)
            {
                if(spin[x][y]!=1 && spin[x][y]!=-1)
                  cout<<"Huge error"<<endl;
                  
                spin[x][y]>0?cout<<" + ": cout<<" - ";
            }
        }
          //getch();
         // system("cls");


        for(i=0;i<10000;i++)
        {
            passed = 1;

            fout<<i<<" "<<n<<"\n";


           if(!(i%100))
              cout<<i<<endl;
            diagonal_update(); //done


           if(!n)
             continue;


        //    do
        //    {

            gen_linked_list(); //done and tested

        /*    cout<<"s[p] is :";
        for(int p=0;p<L;p++)
       {
            cout<<"type = "<<s[p].type<<" and bond = "<<s[p].bond<<" with vxlist["<<p<<"] = "<<vxlist[p]<<endl;
       }

       for(int y=Ly-1;y>=0;y--)
       {
          cout<<endl;
        for(int x=0;x<Lx;x++)
            {
                spin[x][y]>0?cout<<" + ": cout<<" - ";
            }
        }
        getch(); system("cls");*/

            passed = loop_update();

            passed = 1;

        //      }while(!passed);


            equil = adjust_cutoff(i); //done




            /*    cout<<"The state after Monte Carlo Step "<<i+1<<" is : "<<endl;

                for(int y=Ly-1;y>=0;y--)
        {
            cout<<endl;
            for(int x=0;x<Lx;x++)
            {
                spin[x][y]>0?cout<<" + ": cout<<" - ";
            }
        }
          getch();
          system("cls");*/


        }

         fout.close();

        cout<<"Here the value of i = "<<i<<endl;
        cout<<"calculating"<<endl;
        //equilibrated, or so I would like to think.......
        for(int j=0;j<Nbins;j++)
        {
            cout<<" i = "<<i<<endl;
            cout<<i/100 + j<<endl;
            energy=0;
          for(i=0;i<MCS;i++)
         {
             passed = 1;
           // cout<<n<<endl;
            diagonal_update();

            energy+=n;
            //do
            //{

            gen_linked_list();

            passed = loop_update();
            passed = 1;

        //    }while(!passed);

           // equil = adjust_cutoff(i);
         }




         energy/=MCS; // <n>
         data[j]= -1.0*energy/beta + (delta*J*1.0/4.0 + h/(2.0*d)+ epsilon)*Nb*1.0;
         cout<<"Bin "<<j<<" :the energy is "<<data[j]<<endl;
        }

        cout<<"\nEnd"<<endl;
        cout<<"L = "<<L<<endl<<"n = "<<n<<endl<<"NLoops = "<<NLoops<<endl;

        for(int y=Ly-1;y>=0;y--)
        {
            cout<<endl;
            for(int x=0;x<Lx;x++)
            {
                
                if(spin[x][y]!=1 && spin[x][y]!=-1)
                  cout<<"Huge error"<<endl;
                  
                spin[x][y]>0?cout<<" + ": cout<<" - ";
            }
        }

        for(i=0;i<Nbins;i++)
        {
            mean+= data[i];
        }
        mean/=Nbins;
       for(i=0;i<Nbins;i++)
       {
             stddev += (mean-data[i])*(mean-data[i]);
             meandev+= fabs(mean - data[i]);
        }
        stddev = (sqrt(stddev))/Nbins;
        meandev/=Nbins;

        cout<<endl<<"The mean energy is "<<mean<<endl;
        cout<<"The standard deviation in the energy is "<<stddev<<endl;
        cout<<"The mean deviation in the energy is "<<meandev<<endl;

    /*    getch();

        cout<<"s[p] is :";
        for(int p=0;p<L;p++)
       {
            cout<<s[p].type<<" "<<s[p].bond<<" "<<vxlist[p]<<endl;
       }*/



    }





 int main()
 {
    double dummy;
    cout<<"Start"<<endl;
    dummy = rnd(5487);
     Thermal ob;
     ob.run();
     return 0;
 }

void Thermal::diagonal_update()
{
   int p=0,b,i,j,check;
   double aprob=1;
   double dprob=1;

   for(p=0;p<L;p++)
   {
        aprob=0;
        if(s[p].type == 0) // unity operator..... consider addition
        {
             b = int((rnd(0))*Nb);
             
        check = spin[bsites[0][b].x][bsites[0][b].y] + spin[bsites[1][b].x][bsites[1][b].y];
        switch(check) // 0: ud, 2:uu, -2:dd
        {
            case -2:
               aprob =  prob[0];
               break;

            case 0:
               aprob =  prob[1];
               break;

            case 2:
                aprob = prob[2];
                break;

            default :
                cout<<"Adding : Something is seriously wrong....... * ~ * "<<endl;
                exit(0);

        }
            if(aprob)
                 aprob *= beta*Nb*1.0/((L-n)*1.0);
               else
              {
                continue;
               }

           /*    if(aprob<1)
               cout<<"aprob is: "<<aprob<<endl;
            else
               cout<<".";*/

        if(rnd(0)<aprob)
        {
            if(!aprob)
              cout<<"Error";
            s[p].bond = b;
            s[p].type = 1;
            n++;              //cout<<rnd(0)<<" "<<mini(aprob,1.0)<<endl;
        }
      }

        else if(s[p].type == 1) // diagonal.... consider removal
         {
           b = s[p].bond;

         check = spin[bsites[0][b].x][bsites[0][b].y] + spin[bsites[1][b].x][bsites[1][b].y];

        switch(check) // 0: ud, 2:uu, -2:dd
        {
            case -2 :
               dprob =  prob[0];
               break;

            case 0 :
               dprob =  prob[1];
               break;

            case 2 :
                dprob = prob[2];
                break;

            default :
                cout<<"Deleting : Something is seriously wrong....... * ~ * "<<endl;
                exit(0);
        }

        if(dprob)
           dprob = ((L-n+1)*1.0)/(beta*Nb*dprob);
        else
           {
               
               cout<<"The state now is "<<endl;
               
                   for(int y=Ly-1;y>=0;y--)
        {
            cout<<endl;
            for(int x=0;x<Lx;x++)
            {
                if(spin[x][y]!=1 && spin[x][y]!=-1)
                  cout<<"Huge error"<<endl;
                  
                spin[x][y]>0?cout<<" + ": cout<<" - ";
            }
        }
               cout<<" p = "<<p<<endl;
               cout<<"n = "<<n<<" and L = "<<L<<" and dprob = "<<dprob<<endl;
               cout<<"The bond is "<<b<<endl;
               cout<<"The spins are "<<spin[bsites[0][b].x][bsites[0][b].y]<< " and "<<spin[bsites[1][b].x][bsites[1][b].y]<<endl;
               cout<<"The bsites are "<<bsites[0][b].x<<" "<<bsites[0][b].y<<" and "<<bsites[1][b].x<<" "<<bsites[1][b].y<<endl;
           cout<<"Yeh yha ho rha hai yahan pe??"<<endl;
           dprob = 1;
           exit(0);
           }

        if(rnd(0)<dprob)
        {
            s[p].bond = -1;
            s[p].type =  0;
            n--;
            //vxlist[p]=-1;
        }

     }

     else // off - diag.... just propagate the state
      {
          b = s[p].bond;
      //    i =  bsites[0][b].y * Lx + bsites[0][b].x ;
      //    j =  bsites[1][b].y * Lx + bsites[1][b].x ;
        spin[bsites[0][b].x][bsites[0][b].y] *= -1;
        spin[bsites[1][b].x][bsites[1][b].y] *= -1;

      }

   }

}


void Thermal::calc_probloops()
{
     int i,j;
     int ss[4],ssprime[4];
//     int *ss1=&ss[0];int *ss2=&ss[1]; int *tt1=&ss[2]; int *tt2=&ss[3] ;
     int ileg, oleg,vx;

     for(ss[0]=0;ss[0]<=1;ss[0]++)
        for(ss[1]=0;ss[1]<=1;ss[1]++)
           for(ss[2]=0;ss[2]<=1;ss[2]++)
              for(ss[3]=0;ss[3]<=1;ss[3]++)
                {
                   vxcode[ss[0]][ss[1]][ss[2]][ss[3]] = -1;
                }
  //cout<<"debug1\n";
     vxcode[1][0][0][1] = 0;
     vxcode[0][1][1][0] = 1;  // 0 and 1 are off diagonal operator vertices
     vxcode[1][0][1][0] = 2;
     vxcode[0][1][0][1] = 3;  // for such a mapping, i -> (i+1)/2
     vxcode[1][1][1][1] = 4;  // so that -1 -> 0 and 1 -> 1
     vxcode[0][0][0][0] = 5;

     vxwgt[0] = vxwgt[1] = prob[3];
     vxwgt[2] = vxwgt[3] = prob[1];
     vxwgt[4] = prob[2];
     vxwgt[5] = prob[0];



     for(vx=0;vx<6;vx++)
     {
        for(ileg=0;ileg<4;ileg++)
        {
            for(oleg=0;oleg<4;oleg++)
            {
                   vxprob[oleg][ileg][vx]=0;
                for(ss[0]=0;ss[0]<=1;ss[0]++)
                    for(ss[1]=0;ss[1]<=1;ss[1]++)
                      for(ss[2]=0;ss[2]<=1;ss[2]++)
                         for(ss[3]=0;ss[3]<=1;ss[3]++)
                             if(vxcode[ss[0]][ss[1]][ss[2]][ss[3]] == vx)
                                {
                                      for(i=0;i<4;i++)
                                        ssprime[i]=ss[i];

                                      ssprime[oleg] = (ssprime[oleg] + 1 ) %2 ;
                                      ssprime[ileg] = (ssprime[ileg] + 1 ) %2 ;
                                      vxnew[oleg][ileg][vx] = vxcode[ssprime[0]][ssprime[1]][ssprime[2]][ssprime[3]];
                                }

           }
       }
      }

  //    cout<<"debug2"<<endl;
     for(vx=0;vx<6;vx++)
      for(ss[0]=0;ss[0]<=1;ss[0]++)
                    for(ss[1]=0;ss[1]<=1;ss[1]++)
                      for(ss[2]=0;ss[2]<=1;ss[2]++)
                         for(ss[3]=0;ss[3]<=1;ss[3]++)
                             {
                               if(vx==vxcode[ss[0]][ss[1]][ss[2]][ss[3]])
                                   for(i=0;i<4;i++)
                                   vxcodei[i][vx]=ss[i];
                             }


     for(vx=0;vx<6;vx++)
     {
         for(ileg=0;ileg<4;ileg++)
         {
             if(vxnew[0][ileg][vx]!=-1)
              vxprob[0][ileg][vx] = vxwgt[vxnew[0][ileg][vx]];
             else
              vxprob[0][ileg][vx] = 0;

             for(oleg=1;oleg<4;oleg++)
             {
                 if(vxnew[oleg][ileg][vx]!=-1)
                    vxprob[oleg][ileg][vx]=vxwgt[vxnew[oleg][ileg][vx]] + vxprob[oleg-1][ileg][vx];
                 else
                    vxprob[oleg][ileg][vx] = vxprob[oleg-1][ileg][vx];
             }
         }
      }


}



int Thermal::loop_update()
{
    long int p,v,v0,b,btype,k,vx,vnext;
    double prob=0;
    int ss1,ss2,tt1,tt2;
    int s1=0,s2=0,t1=0,t2=0,oleg,ileg;
    long int i,j;
    double rnum;
    unsigned int visited=0;
    double avg=0;
    /*
       The legs vertices are numbered by an index v with v = 4p
       By this method. the legs will be numbered 0-3
      */


avg = 0;
for(j=0;j<100;j++) // This is to know how many loops to construct..... modify with Nloops
{
   visited=0;
    do
    {
        v0 = (int)(rnd(0)*4*L);
    }while(X[v0] < 0);
    v=v0;     // cout<<"v0 = "<<v0<<" and X[v0] = "<<X[v]<<endl; getch();
 do
 {
   // cout<<"Now at vertex v = "<<v<<" at iteration "<<j<<endl;  getch();

    visited++;
/*    if(visited>maxlen)
    {
    //    cout<<"return(0);"<<endl;
    //    cout<<"visited = "<<visited<<endl;
        cout<<"loop broken "<<j<<endl;
        //cout<<j<<endl;
        return 0;
    }  */
    //v=4p, so we can use that here to find the level p and hence the bond
    p=v/4;
    ileg=v%4;
//    cout<<"ileg = "<<ileg<<endl; getch();

//    cout<<"visited = "<<visited<<" and ileg = "<<ileg<<" and v = "<<v<<" with X[v] = "<<X[v]<<" and p = "<<p<<endl;getch();

    vx = vxlist[p]; //cout<<vx;

    if(vxlist[p]==-1) cout<<"not again!";

    rnum = rnd(0);
    for(oleg=0;oleg<4;oleg++)
    {
        prob = 1.0*vxprob[oleg][ileg][vx]/vxprob[3][ileg][vx];
        if(rnum<prob)
           {
               vxlist[p]=vxnew[oleg][ileg][vx];
               break;
           }
        // Now we have found the entry and exit legs
    }


//    cout<<"prob = "<<prob<<" vxlist[p](after the exit) = "<<vxlist[p]<<"oleg = "<<oleg<<" rnum = "<<rnum<<endl;getch();

    if(vxlist[p]==-1)
    {
      cout<<"NOOO:: vxlist[p] = -1!!!!"<<endl;
      cout<<"oleg = "<<oleg<<" ,ileg =  "<<ileg<<" ,vx = "<<vx<<" ,p = "<<p<<" ,v= "<<v<<" ,X[v]= "<<X[v]<<" ,v0 = "<<v0<<endl;
      cout<<"j is "<<j<<" ,n is "<<n<<" ,L is "<<L<<endl;
      cout<<"The random number is "<<rnum<<" , and the prob is "<<prob<<endl;
      exit(0);
    }

      vnext = 4*p + oleg;

  //cout<<" after exiting at oleg = "<<oleg<<"Now at vnext = "<<vnext<<endl;   getch();

    /*  flip the spins at ileg and oleg */
    /* LOGIC : if vfirst[i]!=-1, then spin[i] = vxcodei(fspn(i), lvtx[i]) */

  // cout<<"vnext = "<<vnext<<endl;getch();

   if(s[p].type!=1&&s[p].type!=2)
         cout<<"Oh Dear!!"<<endl;

    if((oleg+ileg)%2!=0)// opposite side operation => operator flip
         {
             (s[p].type==1)?(s[p].type=2):(s[p].type=1);  // flipping the operator type using the fact that only hamiltonians act on s[p]
         }


    if(v0==vnext)
      {
     // cout<<"closed, type 1"<<endl;getch();

       break;
       }

     v=X[vnext];
//  cout<<"now at X[vnext] = "<<X[vnext]<<endl; getch();



 }while(v!=v0);

// cout<<"closed of type 2(ish)"<<endl;getch();
      avg += visited;
      avg = avg*NLoops+visited;
    // cout<<"visited @loop "<<j<< " = "<<visited<<endl;  getch();
//    system("cls");

}
// fix this part by a better averaging method....
//here for example a bounce may occur at step 1 , settiing avg=0, hence also encountering divide by zero segmentation faults

avg=(avg*1.0)/NLoops;
for(i=0;i<N;i++)
         {
             if(v_first[i]!=-1)
             {
                 spin[i%Lx][i/Lx] = (vxcodei[v_first[i]%4][vxlist[v_first[i]/4]] * 2) - 1;
             //    cout<<spin[i%Lx][i/Lx]<<endl;
             }
             else
             {
                 if(rnd(0)<0.5)
                 {
                    spin[i%Lx][i/Lx]*=-1;
                }
             }

            //cout<<fspn[i]<<" "<<(vxcodei[fspn[i]][vxlist[v_first[i]/4]] * 2) - 1<<endl;
         }


//cout<<avg<<endl; getch();
NLoops = (N*1.0/avg) *2;
return 1;
}


void Thermal::gen_linked_list()
{
   int v0,b,f;
   int v1,v2;
   int i1,i2;
   int i,j;
   int ss[4];
   int btype;
   int checkarr[llm];
   
   

   for(int p=0;p<L;p++)
     {
         vxlist[p]=-1;
     }

    for(i=0;i<N;i++)
        {
            v_first[i]=v_last[i]=-1;
            checkarr[i]=0;
            fspn[i]=-1;
         }
    for(i=0;i<4*L;i++)
        X[i]=-1;

   for(int p=0;p<L;p++)
   {
      //cout<<"s["<<p<<"].type = "<<s[p].type<<endl; getch();
       if(s[p].type==0)
        continue;


    v0=4*p;
    b=s[p].bond;

    //cout<<"v0 = "<<v0<<" and s["<<p<<"].bond = "<<s[p].bond<<endl;        getch();

    i1 = (bsites[0][b].y)*Lx + bsites[0][b].x;
    i2 = (bsites[1][b].y)*Lx + bsites[1][b].x;

   //cout<<"i1 = "<<i1<<" and i2 = "<<i2<<endl; getch();
/* Atom i has address:
       ( i%Lx , i/Lx)
    Atom (x,y) has index:
        i = y*Lx + x, where i goes from 0 to N-1       */

    v1 = v_last[i1];
    v2 = v_last[i2];

    //cout<<"v1 = "<<v1<<" and v2 = "<<v2<<endl; getch();

    if(v1 != -1)
    {
       X[v1]=v0;
       X[v0]=v1;
    }
    else
        {
           v_first[i1]=v0;
           fspn[i1]=0;
        }

       if(v2 != -1)
    {
       X[v2]=v0+1;
       X[v0+1]=v2;
    }
    else
       {
         v_first[i2]=v0+1;
         fspn[i2]=1;
       }

       //cout<<"X["<<v1<<"] = "<<X[v1]<<" and X["<<v2<<"] = "<<X[v2]<<endl;
      // cout<<"X["<<X[v1]<<"] = "<<X[X[v1]]<<" and X["<<X[v2]<<"] = "<<X[X[v2]]<<endl; getch();

    v_last[i1] = v0 + 2;
    v_last[i2] = v0 + 3;

    //cout<<"v_ last["<<i1<<"] = "<<v_last[i1]<<"and v_last["<<i2<<"] = "<<v_last[i2]<<endl; getch();


// generating vxlist[p], which is needed for the loop update
    btype=s[p].type;
    if(btype!=1 && btype!=2)
       cout<<"All hell has broken loose"<<endl;
       
  


    //if(btype!=0) // not needed
     
    ss[0]=spin[bsites[0][b].x][bsites[0][b].y] * pow(-1,checkarr[i1]);
    ss[1]=spin[bsites[1][b].x][bsites[1][b].y] * pow(-1,checkarr[i2]);

    if(btype==2)
       {
          checkarr[i1]++;
          checkarr[i2]++;
        }

      for(j=0;j<2;j++)
       ss[j] = (ss[j]+1)/2;

    if(btype==1)
    {
        ss[2]=ss[0];
        ss[3]=ss[1];
    }
    else
    {
        ss[2]=ss[1];
        ss[3]=ss[0];
    }

    vxlist[p] = vxcode[ss[0]][ss[1]][ss[2]][ss[3]];

     

   }

   for(i=0;i<N;i++)
   {
       f = v_first[i];
       if(f!=-1)
       {
           X[f] = v_last[i];
           X[v_last[i]] = f;
    }
   }

 /*for(i=0;i<4*L;i++)
   {
       cout<<"i = "<<i<<" and X[i] = "<<X[i]<<endl;
   }
   getch();
   system("cls");*/

}


void Thermal::init_thermal()
{
   int i=0;

   //These are probs(really weights) for adding and deleting diagonal operators
   prob[0] = epsilon;                          //dd
   prob[1] = epsilon + h/(2.0*d) + delta*J/2.0;   //ud
   prob[2] = epsilon + d*h*1.0/(4.0);                //uu

   prob[3] = J/2.0 ; //off diag vertx

// spins are initalized to the Neel State
   for(i=0;i<N;i++)
    (rnd(0)>0.5)?(spin[i%Lx][i/Lx]=-1):(spin[i%Lx][i/Lx]=1);

   for(int p=0;p<L;p++)
   {
        s[p].bond = -1;
        s[p].type = 0;
        vxlist[p]=-1;
   }



}


void Thermal::test()
{
    int ss[4],i,j,k,vx;
     /*   for(ss[0]=0;ss[0]<2;ss[0]++)
          for(ss[1]=0;ss[1]<2;ss[1]++)
            for(ss[2]=0;ss[2]<2;ss[2]++)
            for(ss[3]=0;ss[3]<2;ss[3]++)
        {
            //for(ss[i]=0;ss[i]<2;ss[i]++)

            cout<<"vxcode["<<ss[0]<<"]["<<ss[1]<<"]["<<ss[2]<<"]["<<ss[3]<<"] = "<<vxcode[ss[0]][ss[1]][ss[2]][ss[3]]<<endl;
        }*/

    for(vx=0;vx<6;vx++)
    {
        for(j=0;j<4;j++)
            ss[j] = vxcodei[j][vx];

        cout<<"Vertex type  = "<<vx<<endl;
        cout<<"Check with vxcode = "<<vxcode[ss[0]][ss[1]][ss[2]][ss[3]]<<endl<<endl<<endl;

        cout<<"\t"<<ss[2]<<"        "<<ss[3]<<endl;
        cout<<"\t"<<       "----------"<<endl;
        cout<<"\t"<<       "|        |"<<endl;
        cout<<"\t"<<       "----------"<<endl;
        cout<<"\t"<<ss[0]<<"        "<<ss[1]<<endl;
        cout<<endl;

    /*    for(i=0;i<4;i++)
        {
            for(k=0;k<4;k++)
            {
                cout<<"vxnew["<<k<<"]["<<i<<"]["<<vx<<"] = "<<vxnew[k][i][vx]<<"\t";
                cout<<"vxprob["<<k<<"]["<<i<<"]["<<vx<<"] = "<<vxprob[k][i][vx]<<endl;
            }
            cout<<endl;
        }

    //    getch();
    //    system("cls");*/

    }



}


void Thermal::makelattice() // to generate bsites : 0 <-> 1; 1<-> 2 ;  bonds numbered from 0 to Nb-1
{
  int b, m;
     for(b=0;b<Nb/d;b++) // loop from 0  -> Nb-1  : Nb is even
     {
         if(b%Lx != Lx-1)
         {
        bsites[0][b].x  = b%Lx;
        bsites[1][b].x  = b%Lx +1;
        bsites[0][b].y  = b/Lx;
        bsites[1][b].y  = b/Lx ;
        }
         else                                                    /* Can be combined in a better code */
         {
        bsites[0][b].x  = b%Lx;
        bsites[1][b].x  = 0;
        bsites[0][b].y  = b/Lx;
        bsites[1][b].y  = b/Lx  ;
        }
     }

    // now b = Nb/2 (16)

    while(b<Nb)
    {
      m = b-(Nb/d);
      if(m/Lx != Ly-1)
      {
          bsites[0][b].x  = m%Lx;
        bsites[1][b].x  = m%Lx ;
        bsites[0][b].y  = m/Lx;
        bsites[1][b].y  = m/Lx +1 ;
      }


    else
    {
           bsites[0][b].x  = m%Lx;
        bsites[1][b].x  = m%Lx;
        bsites[0][b].y  = m/Lx;
        bsites[1][b].y  = 0 ;
    }
    b++;
    }

}



int Thermal::adjust_cutoff(unsigned int loopvar)
{
    static int L_old = L;
     static unsigned int lnew=L,pcheck=loopvar,check=0;


    if((loopvar==pcheck+5000) && (lnew == L) )
    {
        cout<<"Equilibrated"<<endl;
        cout<<"n is "<<n<<" and L is "<<L<<endl;
        check = 1;
        return check;
    }



    if(n>(3*L/4))
    {
        L_old = L;
        L = L + (L/4);
        pcheck = loopvar;
        lnew = L;
    }
    else
        return check;


    for (int p = L_old;p<L;p++)
    {
        s[p].bond = -1;
        s[p].type = 0;
        vxlist[p] = -1;
    }

    return check;
}



double rnd(int seed)
   {
/*    linear congruential random number generator with shuffling */
/*    based on ran1 of second edition of "Numerical Recipes" */

const long int a = 16807;
#define m 2147483647  /* 2^31 - 1  */
#define rm 1.0/m
#define q 127773 /*  m = a*q + p  */
const int p = 2836;
const int n = 32;
#define ndiv (1 + (m-1)/n)
#define rmax (1.0 - 1.2e-7)
      static int r[n+1],r0,r1;
      int j,k;
      if (seed != 0)
/*      initialize table of random numbers  */
         {
          r1 = abs(seed);
       for (j = n+8;j>=1;j--)
         {
         k = r1/q;
         r1 = a*(r1-k*q) - p*k;
         if (r1 < 0) r1 = r1 + m;
         if (j < n) r[j] = r1;
           }
         r0 = r[1];
         }
/*     beginning when not initializing    */
/*     compute r1 = mod(a*r1,m) without overflows  */
       k = r1/q;
       r1 = a*(r1 - k*q) - p*k;
       if (r1 < 0) r1 = r1 + m;
       j = 1 + r0/ndiv;
       r0 = r[j];
       r[j] = r1;
       return min(rm*r0,rmax);
   }


