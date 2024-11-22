// exemple de programmation du gradient conjugué preconditionné

// OK pour le complex 

template<class R,class M,class P> 

int GradienConjugue(const M & A,const P & C,const KN_<R> &b,KN_<R> &x,

                    int nbitermax,double eps)

{

   int n=b.N();

   assert(n==x.N());

   KN<R>   g(n), h(n), Ah(n);

   KN<R> & Cg(Ah);         // on utilise Ah pour stocker Cg  

   g = A*x;  

   g -= b;                // g = Ax-b

   Cg = C*g;              // gradient preconditionné 

   h = -Cg; 



   double g2 = real((Cg,conj(g)));

   if (fabs(g2) < eps*eps)      // solution déjà convergée

   {cout << "iter=0" << " ||g||^2 = " << g2 << endl; 

    return 1;}



   double reps2 = eps*eps*fabs(g2); // epsilon relatif 

   

   for (int iter=0;iter<=nbitermax;iter++)

     {      

       Ah = A*h;     

       

       double ro =  - (real((g,conj(h)))/ real((Ah,conj(h))));  // 

       x += R(ro) *h;

       g += R(ro) *Ah; // plus besoin de Ah, on va stocker Cg 

       Cg = C*g;

       double g2p=g2; 

       g2 = real((Cg,conj(g)));

       cout << iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << endl; 

       if (fabs(g2) < reps2) { 

          cout << iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << endl; 

          return 1;// ok, convergence

          }

       R gamma = g2/g2p;       

       h *= gamma;

       h -= Cg;        

     }



   cout << " Non convergence de la methode du gradient conjugue " <<endl;

   return 0; 

}



template <class R> 

class MatriceIdentite: VirtualMatrice<R> { public:

 typedef typename  VirtualMatrice<R>::plusAx plusAx;

 MatriceIdentite() {}; 

 void addMatMul(const  KN_<R>  & x, KN_<R> & Ax) const {   Ax+=x; } 

 plusAx operator*(const KN<R> &  x) const {return plusAx(this,x);} 

}; 

