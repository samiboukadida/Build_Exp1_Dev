%{ //============ Debut Language C ======================================
/* NOTE: Après modification, regénérer avec
1 - "bison -d exp1.y"
2 - g++ -m64 -fpermissive -g -Wall -Wextra -pedantic -Wno-unused-parameter -c exp1.tab.c -o exp1.o

*/
// Include necessary headers
#include <filesystem>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <algorithm>  // For std::replace
#include <cassert>
#include <cmath>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string.h>
//#include <sstream>
  //#include "matl2d.hpp"
  //#include "assertion.hpp"
#include "sfem.hpp"
#include "RNM.hpp"
#include "GC.hpp"

  //#include "RNM_opc.hpp"
  //#include "gmres.hpp"
#include <cctype>
#include <cstdio>
#include <cstring>
#define CHECK_KN

  //For Handling File Path
#include "LexicoEdge.hpp"
#include "HeapSort.hpp"

#include <stdio.h>
#include <stdlib.h>

  using namespace std;
  // Define a global variable to store the mesh file path

  // Declare mesh_load function
  std::string mesh_load(const char* filepath);

  //void InternalError(const char * str) {cerr << str; exit(1);}
  const R L=40.;
#include "gnuplotiso.hpp"

  istream * ccin; // la variable entrer
  typedef double R;
  //Mesh Th("carre.msh"); 
  R xy;
  //istream * fcin (& cin ) ; //fichier input
  // double* fxy = new double [Th.nv]; ;
  char * cmsh=new char[1024];
  char * cf=new char[1024];
  char * cfp=new char[1024];
  char * cg=new char[1024];


  inline void yyerror(const char * msg)
  {  cerr << msg << endl; exit(1);}

  typedef  struct elemrec elemrec;

  struct elemrec
  {
    char *name; 		// name of elembol 
    int type; 			// type of elembol: either VAR or FNCT 
    union {
      double var; // value of a VAR 
      double (*fnctptr)(double); // value of a FNCT 
    } value;
    struct elemrec *next; 	// link field 
  };

  extern elemrec *elem_table;
  elemrec * putelem (const char *elem_name,int elem_type);
  elemrec *getelem (const char *elem_name);
  int yylex(); 
  void init_table();

%}
//============= Fin Language C ===========================================

%union{
  double val; 	 
  elemrec *tptr;
  char *str;  /* Add this for handling strings */
}

%token  EOFFILE
%token  <val> NUM 		
%token  <tptr> VAR FNCT
%token  <str> STRING  /* Associate STRING with the char* type */
%type   <val> exp
%type   <str> file_expr /* Declare the type of file_expr */
%left   '-' '+'
%left   '*' '/'
%left   NEG 			 
%right  '='
%right  '^' 			
 // %token STRING /*For Handling File Path*/
%%
//============ Define the grammar rules here For Handling File Path ================
//==================================================================================
input:
exp EOFFILE { cout << $1 << endl; return 0; }
| line
| file_expr '\n' { /* Handle the file here */ }
;

file_expr:
STRING  {
  /* STRING will contain the path */
  /* Pass this string (file path) to the Mesh loading logic */
  mesh_load($1);  // Ensure this function correctly handles absolute paths
}
;

//input: exp EOFFILE { cout << $1 << endl; return 0;}	// empty
/*		| input line */

line:
'\n'
| exp '\n' { printf ("\t%.10g\n", $1); }
| error '\n' { yyerrok; }
;
exp: 	  NUM { $$ = $1; }
| VAR { $$ = $1->value.var; }
| VAR '=' exp { $$ = $3; $1->value.var = $3;xy=$$; cout << $1->name <<"="<<$1->value.var<<endl;}
| FNCT '(' exp ')' { $$ = (*($1->value.fnctptr))($3); }
| exp '+' exp { $$ = $1 + $3;         }
| exp '-' exp { $$ = $1 - $3; }
| exp '*' exp { $$ = $1 * $3; }
| exp '/' exp { $$ = $1 / $3; }
| '-' exp %prec NEG { $$ = -$2; }
| exp '^' exp { $$ = pow ($1, $3); }
| '(' exp ')' { $$ = $2; }
;
%%
//==================================================================================
//==================== End of grammar ==============================================

//==================== Debut Language C++ ==========================================
/*
*
*/
std::string mesh_load(const char* filepath) {
  std::string meshPath(filepath);
  // Output the original mesh file path for debugging
  std::cout << "Mesh file path: " << meshPath << std::endl;
  // Check if the path starts with "/c/" or similar (Cygwin style) and convert to Windows style
  if (meshPath[0] == '/' && meshPath[2] == '/') {
    // Replace "/c/" with "C:\" and similar patterns
    meshPath[1] = toupper(meshPath[1]);  // Uppercase the drive letter
    meshPath = meshPath.substr(1);  // Remove the leading slash
    std::replace(meshPath.begin(), meshPath.end(), '/', '\\');  // Convert slashes to backslashes
  }
  // Output the converted mesh file path for debugging
  std::cout << "Converted mesh file path: " << meshPath << std::endl;

  return meshPath;
}

elemrec * putelem (const char *elem_name,int elem_type){
  elemrec *ptr;
  ptr = (elemrec *) malloc (sizeof (elemrec));
  ptr->name = (char *) malloc (strlen (elem_name) + 1);
  strcpy (ptr->name,elem_name);
  ptr->type = elem_type;
  ptr->value.var = 0; // set value to 0 even if fctn. 
  ptr->next = (struct elemrec *)elem_table;
  elem_table = ptr;
  return ptr;
}

elemrec *getelem (const char *elem_name){
  elemrec *ptr;
  for (ptr = elem_table; ptr != (elemrec *) 0; ptr = (elemrec *)ptr->next)
    if (strcmp (ptr->name,elem_name) == 0)	return ptr;

  return 0;
}

int  yylex ()
{
  int c;
  // Ignore whitespace, get first nonwhite character. 
  while ((c = ccin->get()) == ' ' || c == '\t');

  if (c == EOF)   return EOFFILE;

  // Char starts a number => parse the number. 
  if (c == '.' || isdigit (c)) {
    ccin->putback(c);
    //scanf ("%lf", &yylval.val);
    *ccin >>yylval.val;
    return NUM;
  }
  // Char starts an identifier => read the name. 
  if (isalpha (c)) {
    elemrec *s;
    static char *elembuf = 0;
    static int length = 0;
    int i;
    // Initially make the buffer long enough for a 40-character elembol name. 
    if (length == 0) length = 40, elembuf = (char *)malloc (length + 1);
    i = 0;
    do {
      // If buffer is full, make it bigger. 
      if (i == length)
        {
          length *= 2;
          elembuf = (char *)realloc (elembuf, length + 1);
        }
      // Add this character to the buffer. 
      elembuf[i++] = c;
      // Get another character. 
      c = ccin->get ();
    } while (c != EOF && isalnum (c));
    ccin->putback(c);
    elembuf[i] = '\0';
    s = getelem (elembuf);
    if (s == 0)
      s = putelem (elembuf, VAR);
    yylval.tptr = s;
    return s->type;
  }
  // Any other character is a token by itself.
  //	xy=yylval.tptr ->value.var;cout<<"je suis la  "<<"xy="<<xy<<endl;
  //cout<<"c="<<c<<endl;
  return c;
}

struct init  {
  const char *fname;
  double (*fnct)(double);
};
															
init arith_fncts[] = {
    { "sin", sin },
    { "cos", cos },
    { "atan", atan },
    { "ln", log },
    { "exp", exp },
    { "sqrt", sqrt },
    { nullptr, nullptr }  // Dernière ligne pour signaler la fin
};


// The elembol table: a chain of 'struct elemrec'. 
					
elemrec *elem_table = (elemrec *)0;

void init_table () 		// puts arithmetic functions in table. 
{ 
  int i;
  elemrec *ptr;
  ptr = putelem("pi", VAR);
  ptr->value.var=4.*atan(1.);  
  
  for (i = 0; arith_fncts[i].fname != 0; i++){
    ptr = putelem (arith_fncts[i].fname, FNCT);
    ptr->value.fnctptr = arith_fncts[i].fnct;
  }
}
/*
* Calculer f(x,y,u)
*/ 
double   interprete (double x,double y,double u, char *f) { 
  xy=0.;
  ccin = new istringstream(f);
  elemrec *ptr;
  ptr = putelem("x", VAR);
  ptr->value.var=x;
  ptr = putelem("y", VAR);
  ptr->value.var=y;
  ptr = putelem("u", VAR);
  ptr->value.var=u;
  yyparse ();
  ptr = putelem("f", VAR);
  return xy;   
  // return ptr->value.var ;
  delete ccin;
  ccin=&cin;
}

//=========================== 
//=======CLASSES=======
//===========================
class MatriceSum2: public VirtualMatrice<R> { public:
  typedef VirtualMatrice<R>::plusAx plusAx;
  typedef VirtualMatrice<R> M;
  const  M & A;
  const  M & B;
  MatriceSum2(const  M & AA,const  M & BB ) : A(AA),B(BB){}
  void addMatMul(const  KN_<R>  & x, KN_<R> & Ax) const {
    Ax += A* x;
    Ax += B* x;
  }
  plusAx operator*(const KN<R> &  x) const {return plusAx(this,x);}
}; 

class MatLap2d: public VirtualMatrice<R> { public:
  //  Matrice $  \beta \nabla u \nabla v $
  const  Mesh & Th;
  const R alpha;
  const R beta;
  const R ccr;
  const R zeta;
  const KN<R> un;
  //const  char * cfp;

  typedef  VirtualMatrice<R>::plusAx plusAx;
  MatLap2d(const  Mesh & T,R aa, R b,R rr,R z, KN<R> unn) : Th(T),alpha(aa),beta(b),ccr(rr), zeta(z),un(unn) {};
  void addMatMul(const  KN_<R>  & x, KN_<R> & Ax) const {
    if ( ccr ) {
      for (int ie=0 ;ie<Th.neb;ie++){ 
        const BoundaryEdge & E(Th.bedges[ie]);
        const Vertex & v0=E[0];
        const Vertex & v1=E[1];
        int lab=E.lab,i0(Th(v0)),i1(Th(v1));
        R x0=x[i0],x1=x[i1];
        R  l( E.length());
        R rrr=l*ccr;
        if (lab==4) Ax[i0] += (x0/3+x1/6)*rrr;  
        if (lab==4) Ax[i1] += (x0/6+x1/3)*rrr;
      }         
    }
    //nouveau term
    if ((zeta) && (std::string(cfp) != "")){  
      for (int k=0;k<Th.nt;k++){   
        const Triangle & K(Th[k]);
        int i0(Th(K[0])),i1(Th(K[1])),i2(Th(K[2]));
        R x0  = x[i0], x1=x[i1], x2=x[i2];
        R un0 = un[i0], un1=un[i1], un2=un[i2];
        R xm0 = (x1+x2)/2, xm1=(x0+x2)/2, xm2=(x0+x1)/2;
        R unm0= (un1+un2)/2, unm1=(un0+un2)/2, unm2=(un0+un1)/2;   
        R cz  = K.area*zeta/6;

        R xvm0=(Th(i1).x+Th(i2).x)/2; R  yvm0=(Th(i1).y+Th(i2).y)/2;
        R xvm1=(Th(i0).x+Th(i2).x)/2; R  yvm1=(Th(i0).y+Th(i2).y)/2;
        R xvm2=(Th(i0).x+Th(i1).x)/2; R  yvm2=(Th(i0).y+Th(i1).y)/2;
  
        R  fpum0=interprete( xvm0, yvm0, unm0,cfp) ;//valeur de la fonction au point m0=milieu de i1 et i2
        R  fpum1=interprete(xvm1, yvm1, unm1,cfp) ;
        R  fpum2=interprete(xvm2, yvm2, unm2,cfp) ;
        if ( K[0].lab ==0 )  Ax[i0] += - ( fpum1*xm1 + fpum2*xm2 )*cz;// ( fs(unm1)*xm1 + fs(unm2)*xm2 )*cz
        if ( K[1].lab ==0 )  Ax[i1] += - ( fpum0*xm0 + fpum2*xm2 )*cz;//( fs(unm0)*xm0 + fs(unm2)*xm2 )*cz;
        if ( K[2].lab ==0 )  Ax[i2] += - ( fpum0*xm0 + fpum1*xm1 )*cz;// ( fs(unm0)*xm0 + fs(unm1)*xm1 )*cz;
      }  
    }

    if ( alpha || beta ) {
      for ( int k=0 ; k<Th.nt ; k++ ) {
          const Triangle & K(Th[k]);
          int i0(Th(K[0])),i1(Th(K[1])),i2(Th(K[2]));//numero globaux des 3 sommets
          R2 H0(K.H(0)),H1(K.H(1)),H2(K.H(2));
          R x0=x[i0],x1=x[i1],x2=x[i2];
          R2 gradx= H0*x0 + H1*x1 + H2*x2;
          R cgg=K.area*alpha;
          R cc=K.area*beta;
          if ( alpha  ) {
            if ( K[0].lab ==0 ) Ax[i0] +=  (gradx,H0)*cgg;
            if ( K[1].lab ==0 ) Ax[i1] +=  (gradx,H1)*cgg;
            if ( K[2].lab ==0 ) Ax[i2] +=  (gradx,H2)*cgg;
          }

          if ( beta ) {
            if ( K[0].lab ==0 )  Ax[i0] += (x0/6+x1/12+x2/12) *cc;
            if ( K[1].lab ==0 )  Ax[i1] += (x0/12+x1/6+x2/12) *cc;
            if ( K[2].lab ==0 )  Ax[i2] += (x0/12+x1/12+x2/6) *cc;
          }
      }
    }
  }
  plusAx operator*(const KN<R> &  x) const {return plusAx(this,x);}
};
//======== FIN CLASSES ===============

//======== Fonction adaptmesh ========
//====================================
void adaptmesh(Mesh & Th,KN<R> un){
  //___calcul d'indice d'erreur nu____
  cout<<"cf="<<cf<<"cfp="<<cfp<<"cg="<<cg<<endl;
 
  KN<R> nu(Th.nt);nu=0.;
  for (int k=0;k<Th.nt;k++){ 
      Triangle & K(Th[k]);
      int i0(Th(K[0])),i1(Th(K[1])),i2(Th(K[2]));
      R Gx=(K[0].x+K[1].x+K[2].x)/3 ,Gy=(K[0].y+K[1].y+K[2].y)/3;
      R un0=un[i0],  un1=un[i1] , un2=un[i2];
      R unG=(un0+un1+un2)/3;
      R Nf=sqrt(K.area*pow(interprete(Gx,Gy,unG,cf),2));//fini norm de f sur K
      //calcul du saut:
      R2 H0(K.H(0)),H1(K.H(1)),H2(K.H(2));
      R2 gradun= H0*un0 + H1*un1 + H2*un2;
      //calcul de nk vecteur unitaire sortant:
      R saut=0,hK=Max(K.lenEdge(0),K.lenEdge(1),K.lenEdge(2));//calcul de hK
      for (int a=0;a<3;a++){
        if(Th.adj[3*k+a] != -1){
          int ttt=Th.adj[3*k+a], tt=ttt/3, aa=ttt%3;
          Triangle & KK(Th[tt]);
          R2 E(KK.H(aa));
          R2 nka=E/sqrt((E,E));
          R gn1=abs((gradun,nka));//gradient normale
          R2 HH0(KK.H(0)),HH1(KK.H(1)),HH2(KK.H(2));//calcul de grad pour KK
          int i0(Th(KK[0])),i1(Th(KK[1])),i2(Th(KK[2]));
          R un0=un[i0];  un1=un[i1] ; un2=un[i2];
          R2 graduna= HH0*un0 + HH1*un1 + HH2*un2;
          R gn2=abs((graduna,-nka));//gradient normale adj
          R ha=K.lenEdge(a);//longeur arete courante
          saut += pow(ha*(gn2-gn1),2);
          if (aa==2) saut=sqrt(saut);
	      }
	    }
    
      nu[k]=hK*Nf+saut;//nuK
      cout<<"nu["<<k<<"]="<<nu[k]<<endl;
  }
  //____calcul de tolerance tol=coef*moyenne des nuK___ici coef=1
  R tol=0.;
  for (int k=0;k<Th.nt;k++) 
    tol += nu[k];

  tol=tol/Th.nt;
  cout<<"tol="<<tol<<endl;
  
  //----- Tableau d'arêtes: le ----
  //-------------------------------
  int nt3 = Th.nt*3;
  LexicoEdge *le = new LexicoEdge[nt3];//tab d'aretes avec  doublons

  for (int k=0,kk=0;k<Th.nt;k++)
    for (int i=0;i<3;i++)
      le[kk++]=LexicoEdge(k,i,Th);

  cout<<"____avant trie le["<<nt3-1<<"].t ="<<le[nt3-1].t <<
    "    le["<<nt3-1<<"].e ="<<le[nt3-1].e <<
    "    le["<<nt3-1<<"].i0="<<le[nt3-1].i0<<
    "    le["<<nt3-1<<"].i1="<<le[nt3-1].i1<<endl;
  

  HeapSort(le,nt3);
  cout<<"____apres trie le["<<0<<"].t ="<<le[0].t <<
    "    le["<<0<<"].e ="<<le[0].e <<
    "    le["<<0<<"].i0="<<le[0].i0<<
    "    le["<<0<<"].i1="<<le[0].i1<<endl;
  cout<<"____apres trie le["<<1<<"].t ="<<le[1].t <<
    "    le["<<1<<"].e ="<<le[1].e <<
    "    le["<<1<<"].i0="<<le[1].i0<<
    "    le["<<1<<"].i1="<<le[1].i1<<endl;

  //----Triangles adjacents----
  //---------------------------
  cout<<"triangle adj"<<endl;
  for (int i=0;i<nt3;i++)
    cout<<"adj["<<i<<"] = "<<Th.adj[i]<<endl;
  
  //_____ tableau  d'aretes sans doublons ______
  int nbarete = nt3;
  for (int i=1;i<nt3;i++)
    if ( le[i-1]==le[i] ) nbarete--;//nb d'aretes 

  LexicoEdge *lee = new LexicoEdge[nbarete];//tab d'aretes sans doublons
  int jj=0;
  for (int i=1;i<nt3;i++) { 
    if (!( le[i-1] ==le[i]) ){
      lee[jj]=le[i-1];jj =jj+1;
	    lee[jj]=le[i];
	  }
  }//____ fin tableau d'aretes sans doublons____
  //____________________________________________
  
  //--------constructions des arêtes adj aux triangles avec nuK>tol => arêtes à decouper --> cutedg
  //--------procedure 1
  int *cutedge = new int[nbarete];//sa taille reelle = nbcutedge ( < nbarete )
  int k=0;
  for (int i=0;i<nbarete;i++) { 
    int t=lee[i].t;
    int e=lee[i].e;
    if ( nu[t] > tol ) {//critère pour découper l'arête
    
      cutedge[k]=3*t+e ; 
      k=k+1;
      cout<<" arêtes à découper"<<endl;
      cout<<"lee["<<i<<"].t ="<<lee[i].t <<"    lee["<<i<<"].e ="<<lee[i].e<<"    lee["<<i<<"].i0 ="<<lee[i].i0<<"    lee["<<i<<"].i1 ="<<lee[i].i1 <<endl;

    }
  }
  int nbcutedge=k;//nb des arêtes à découper
  cout<<"nb des arêtes="<<nbarete<<endl; 
  cout<<"nb des  arêtes à découpe="<<nbcutedge<<endl; 
  for (int i=0;i<nbarete;i++)
    cout<<"lee["<<i<<"].t ="<<lee[i].t  <<"    lee["<<i<<"].e ="<<lee[i].e<<"    lee["<<i<<"].i0 ="<<lee[i].i0<<"    lee["<<i<<"].i1 ="<<lee[i].i1 <<endl;

  Th.decoupe(cutedge,nbcutedge);
}
//========Fin Fonction adaptmesh========

/*
* Helper function to convert Cygwin-style path to Windows-style path
*/
std::string convertCygwinPath(const std::string& cygwinPath) {
  std::string windowsPath = cygwinPath;
  if (windowsPath.rfind("/c/", 0) == 0) {
    windowsPath.replace(0, 3, "C:/");
  }
  std::replace(windowsPath.begin(), windowsPath.end(), '/', '\\');
  return windowsPath;
} 
/*
* Fonction pour normaliser les chemins avec des barres inverses simples
*/ 
std::string normalizePath(const std::string& rawPath) {
    std::string normalized = rawPath;
    // Remplacer les barres inverses simples '\' par des '/'
    std::replace(normalized.begin(), normalized.end(), '\\', '/');
    return normalized;
}
/*
* Fonction pour extraire le nom de base sans extension
*/
std::string getBaseName(const std::string& rawPath) {
    // Normaliser le chemin d'entrée
    std::string normalizedPath = normalizePath(rawPath);
    // Créer un std::filesystem::path avec le chemin normalisé
    std::filesystem::path path(normalizedPath);
    // Retourner le nom du fichier sans extension
    return path.stem().string();
}

//==================================================================
//========= main de exp1.y =========================================
//==================================================================

int main (int argc, char** argv){ 
  ofstream ftab("ftab.txt");
  ccin = & cin; // 
  init_table ();
  ifstream xx(argv[1]);  // Read "data.txt"
  xx >> cmsh >> cf >> cfp >> cg;
  xx.close();

  // Convert the mesh path using mesh_load (returns std::string)
  std::string convertedPath = mesh_load(cmsh);

  std::cout << "Converted mesh file path: " << convertedPath << std::endl;
  std::cout <<"mesh: " << cmsh <<" ________ f: "<< cf <<" ________ fp: "<< cfp <<"_______g:  "<< cg<< std::endl;
  // Now pass the converted path to the Mesh constructor
  // Use `.c_str()` to convert std::string to const char*
  // Extraction du nom de fichier sans extension
  std::string baseName = getBaseName(convertedPath);
  std::cout << "Nom de fichier extrait (baseName) : " << baseName << std::endl;

  Mesh Th(convertedPath.c_str(), 10.0);

  for ( int adapt=0 ; adapt<3 ; adapt++ ) {//RESOLUTION 3 FOIS

    //resolution
    cout<<"resolution"<<endl;
    KN<R> un(Th.nv);un=0.;
    KN<R> vn(Th.nv);vn=0.;
    KN<R> b(Th.nv);b=0.;  //  $ b[i] =  \int_\Omega  f_h w_i $
    KN<R> b1(Th.nv);b1=0.;
    KN<R> b2(Th.nv);b2=0.;

    MatLap2d A1(Th,1.,0.,0.,0.,0);

    KN<R> x(Th.nv);

    //CL un :Dirichlet
    // for (int k=0;k<Th.nv;k++)   if(Th(k).lab==1)    un[k] = T0/Te;
    for (int k=0;k<Th.nv;k++) //g(Th(k))  
      if(Th(k).lab !=0)    un[k] = interprete(Th(k).x,Th(k).y,un[k],cg)  ; 

    R eps=1.e-3, err=2.;
    int iter=0,itermax=50;
    while ((err>eps)&&(iter<itermax)) {
      cout<<"calcul de A2"<<endl;  
      //__________ Matrice du Laplacien
      MatLap2d A2(Th,0.,0.,0.,1.,un); //Matrice M0
      MatriceIdentite<R> Id; // Matrice Identite
      MatriceSum2  A(A1,A2);

      //===== construction de b: second membre =========
      //___________ calcul de  b1 ______________________
      MatLap2d M1(Th,1.,0.,0.,0.,0); b1= M1*un;

      //___________ calcul de b2 _______________________
      cout<<"calcul de b2____"<<endl;
      vn=0.;
      for (int k=0;k<Th.nv;k++)   vn[k]=-interprete(Th(k).x,Th(k).y,un[k],cf)  ;
      //-f(un[k])  ;     //fp(un[k])-xi-beta;
      MatLap2d M2(Th,0.,1.,0.,0.,0); b2= M2*vn;

      //___________ calcul de b3 pour gamma robin ______
      // MatLap2d M3(Th,0,0,rr,0,0); b3= M3*un;
      //====== Valeur du Second membre b ===============
      b=b1+b2;
      
      //==================== resolution ==================
      //==================================================
      x=0.;  //  donne initial  avec les conditions aux limites.
      GradienConjugue(A,Id, b,x,Th.nv,1e-10);
      un = un-x ;  //un+1=un-wn
      //||un+1-un||
      {err=abs(x[0]) ;for (int k=0;k<Th.nv;k++)  if (!( abs(x[k]) < err) )  err=abs(x[k]) ;}

      iter +=1;
      if  (err<=eps) cout<<"sortie avec  ERREUR="<<err<<endl;
      if  (iter>=itermax) cout<<"sortie avec iter ="<< iter <<endl;
    }

    x = un;
    const char * graph=new char[10];
    const char * xsol=new char[10];
    const char * plotiso=new char[10];
    { // a file for gnuplot
  
      //solution approche
      cout<<"je suis a gnuplot"<<endl;
  
      if (adapt==0){ graph=(char*)"uh1";xsol=(char*)"x1.sol";plotiso=(char*)"plotiso1";}
      if (adapt==1){ graph=(char*)"uh2";xsol=(char*)"x2.sol";plotiso=(char*)"plotiso2";}
      if (adapt==2){ graph=(char*)"uh3";xsol=(char*)"x3.sol";plotiso=(char*)"plotiso3";}
      //creation des fichiers:  uh1, uh2, uh3
      ofstream plot(graph);
      for(int it=0;it<Th.nt;it++)
        plot  << (R2) Th[it][0] << " " << x[Th(it,0)]<< endl
              << (R2) Th[it][1] << " " <<x[Th(it,1)] << endl
              << (R2) Th[it][2] << " " << x[Th(it,2)] << endl
              << (R2) Th[it][0] << " " << x[Th(it,0)] << endl << endl << endl;

    }
    // creation des fichiers :  x1.sol, x2.sol, x3.sol
    {
      ofstream f(xsol);
      f << x << endl;
    } 

    gnuplot_iso(plotiso,Th,x,15);

    /* for (int k=0;k<Th.nt;k++) */
    /*   {  Triangle & K(Th[k]); */
    /* 	int i0(Th(K[0])),i1(Th(K[1])),i2(Th(K[2])); */
    /* 	cout<<"Th["<<k<<"]=  "<<i0<<"   "<<i1<<"    "<<i2<<"   "<<K.lab<<endl;} */  
    //=======cout================================================================
    cout <<"fin de résolution pour adapt="<<adapt<<endl;
    cout <<"fin de résolution "<<(adapt+1)<< " fois"<<endl;
    for (int k=0;k<Th.nt;k++) {
        Triangle & K(Th[k]);
        int i0(Th(K[0])),i1(Th(K[1])),i2(Th(K[2]));
        cout<<"Th["<<k<<"]=  "<<i0<<"   "<<i1<<"    "<<i2<<"   "<<K.lab<<endl;
	  }   
    cout<<"Th.nv="<<Th.nv<<"   Th.nt="<<Th.nt<<"   Th.neb="<<Th.neb<<endl;    
    for(int i=0;i<Th.nv;i++)
      cout<<"Th("<<i<<")=  "<<Th(i)<<endl;

    //======= Fin cout ==========================================================
    /* 
     * ==========================================
     * Ajout de l'export des fichiers de maillage 
     * ==========================================
     * Date : 18 novembre 2024
     * Développeur : Sami
     *
     * Description :
     * - Après chaque étape d'adaptation de maillage, le programme exporte un fichier 
     *   représentant le maillage adapté avec les mêmes conventions que le fichier
     *   d'origine.
     * - Le nom des fichiers exportés est généré dynamiquement pour refléter :
     *   1. Le nom d'origine du fichier de maillage.
     *   2. La méthode utilisée (ici `adaptmesh`).
     *   3. La date et l'heure de l'export.
     *   4. Le numéro d'étape d'adaptation (1, 2, etc.).
     *
     * Exemple :
     * - Pour un fichier `carre3.msh` traité le 18 novembre 2024 à 14h30min45s, les fichiers générés seront :
     *   - `carre3_adaptmesh_2024_11_18_14_30_45_1.msh` (après la 1ère adaptation)
     *   - `carre3_adaptmesh_2024_11_18_14_30_45_2.msh` (après la 2ème adaptation)
     *
     * Méthodes impliquées :
     * - `Mesh::exportMesh` : Ajoutée pour écrire un fichier de maillage adapté.
     * - `getBaseName` : Ajoutée pour extraire le nom de base du fichier d'origine.
     *
     * L'objectif de cette modification est de faciliter le suivi des maillages adaptés,
     * tout en conservant une organisation claire et intuitive des fichiers générés.
     */
    //========= Génération du nom du fichier avec date et heure==================
    
    auto now = std::chrono::system_clock::now();
    std::time_t currentTime = std::chrono::system_clock::to_time_t(now);
    std::stringstream dateTimeStream;
    dateTimeStream << std::put_time(std::localtime(&currentTime), "%Y_%m_%d_%H_%M_%S");

    std::string exportFileName = baseName + "_adaptmesh_" + dateTimeStream.str() + "_" + std::to_string(adapt + 1) + ".msh";

    //========= Exporter le maillage ============================================
    Th.exportMesh(exportFileName);

    // Output for verification
    std::cout << "Base Name: " << baseName << std::endl;
    std::cout << "Export File Name: " << exportFileName << std::endl;

    //========= Appel au découpage des arêtes ===== Adaptation du maillage ======    
    if (adapt != 2)   adaptmesh(Th,un);    
  
  }//========fin de résolution 'adapt' fois =====================================

  cout <<"FIN DE RESOLUTION 3 FOIS"<<endl;

  ftab.close();  
  delete cmsh;delete cf;   
  delete cfp;delete cg;
  return 0;
}