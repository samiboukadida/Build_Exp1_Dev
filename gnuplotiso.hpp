#include <cstring>
// La fonction qui sauvegarde les isolignes pour {\em gnuplot}
// ------------------------------------------
void gnuplot_iso(const char * fileiso, Mesh & Th, KN_<R> & U ,int niso)
{  assert(Th.nv == U.size() && niso != 0);

 ofstream f(fileiso); // fichier des isolignes

 char * fgnuplot=new char[256]; // le nom du script gnuplot
 strcpy(fgnuplot,fileiso);      // on copie le nom de "fileiso"
 strcat(fgnuplot,"-gnuplot");   // on rajoute l'extension "-gnuplot"

 ofstream fl(fgnuplot); // fichier script pour gnuplot
 fl << "plot '"<<fileiso<<"'  w l " <<endl; // première commande du script 

 // sélection des valeurs des isolignes
 R umax=U.max(), umin=U.min();
 if((umax == umin) | (niso == 0)) 
   niso=1;

for(int imu=0;imu<=niso;imu++)
{   R muc=umin+(umax-umin)*imu/niso;
    cout<<"Isoligne No "<<imu<<" = "<<muc<<endl;
    R2 xylast;
    for(int k = 0; k<Th.nt; k++)
    { int m=0;
      for(int i1=0;i1<3;i1++)
      { int i = Th(k,i1);
        int j = Th(k,(i1+1)%3);
        R a=2.;
        if((U[i] - U[j])!=0)
          a = (muc - U[j]) / (U[i] - U[j]);
        else
         {if(muc ==U[j])
	   {f<<(R2)Th.vertices[i]<<endl;//f<<L*(R2)Th.vertices[i]<<endl;
             a=0.;}
         }

        if((a>=0)&&(a<=1))
        { m++;
          xylast= (R2)Th.vertices[i]*a+(R2)Th.vertices[j]*(1.-a); 
          xylast=L*xylast;
          f<<xylast<<endl;}
        }
      if(m)
        f<<endl<<endl;
    }
    fl<<"set label "<<imu+1<<" '"<<muc<<"' at " << xylast.x <<","<<xylast.y <<endl;
}
    f.close();
    fl<<"replot"<<endl;fl.close();
}
