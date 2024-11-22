// vérification d'allocation
#include <cmath>
#include <iostream>
#include "assertion.hpp"
#include <cstdlib>
// quelques fonctions utiles
template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}
template<class T> inline T Abs (const T &a){return a <0 ? -a : a;}
template<class T> inline double Norme (const T &a){return sqrt(a*a);}
template<class T> inline void Exchange (T& a,T& b) {T c=a;a=b;b=c;}
template<class T> inline T Max (const T &a,const T & b,const T & c){return Max(Max(a,b),c);}
template<class T> inline T Min (const T &a,const T & b,const T & c){return Min(Min(a,b),c);}

using namespace std;
// Définition de R
typedef double R;

// La classe R2
class R2 {
  friend ostream& operator <<(ostream& f, const R2 & P )
  { f << P.x << ' ' << P.y   ; return f; }
  friend istream& operator >>(istream& f,  R2 & P)
  { f >>  P.x >>  P.y  ; return f; }
  
public:  
  R x,y;
  R2 () :x(0),y(0) {}
  R2 (R a,R b):x(a),y(b)  {}
  R2 (R2 a,R2 b):x(b.x-a.x),y(b.y-a.y)  {}
  R2 operator+(R2 P)const   {return R2(x+P.x,y+P.y);}
  R2 operator+=(R2 P)  {x += P.x;y += P.y;return *this;}
  R2 operator-(R2 P)const   {return R2(x-P.x,y-P.y);}
  R2 operator-=(R2 P) {x -= P.x;y -= P.y;return *this;}
  R2 operator-()const  {return R2(-x,-y);}
  R2 operator+()const  {return *this;}
  R  operator,(R2 P)const  {return  x*P.x+y*P.y;} // produit scalaire
  R  operator^(R2 P)const {return  x*P.y-y*P.x;} // produit mixte
  R2 operator*(R c)const {return R2(x*c,y*c);}
  R2 operator/(R c)const {return R2(x/c,y/c);}
  R2 perp() {return R2(-y,x);} // la perpendiculaire
};
inline R2 operator*(R c,R2 P) {return P*c;}


// La classe Label (réferences de sommets ou triangles) 
class Label {   public: 
  int lab;
  Label(int r=0):lab(r){}
  int onGamma() const { return lab;} 
};

inline ostream& operator <<(ostream& f,const Label & r  )
  { f <<  r.lab ; return f; }
inline istream& operator >>(istream& f, Label & r  )
  { f >>  r.lab ; return f; }



// La classe Vertex (modélisation des sommets) 
class Vertex : public R2,public Label { public:
  Vertex() : R2(),Label(){};
  Vertex(R2 P,int r=0): R2(P),Label(r){}
private:
  Vertex(const Vertex &);
  void operator=(const Vertex &);
};

inline ostream& operator <<(ostream& f, const Vertex & v )
  { f << (R2) v << ' ' << (Label &) v   ; return f; }
inline istream& operator >> (istream& f,  Vertex & v )
  { f >> (R2 &) v >> (Label &) v ; return f; }


// La classe BoundaryEdge (arêtes frontièrs)
class BoundaryEdge: public Label {
public:
  Vertex *vertices[2];
  void set(Vertex * v0,int i0,int i1,int r)
  { vertices[0]=v0+i0; vertices[1]=v0+i1; lab=r; }
  bool in(const Vertex * pv) const {return pv == vertices[0] || pv == vertices[1];}
  BoundaryEdge(){}; // constructor par défaut vide 
  void Draw() const;
  Vertex & operator[](int i) const {ASSERTION(i>=0 && i <2);
  return *vertices[i];}
  R length() const { R2 AB(*vertices[0],*vertices[1]);return sqrt(AB.x * AB.x + AB.y * AB.y);}
private:
  BoundaryEdge(const BoundaryEdge &);   // interdit la construction par copie
  void operator=(const BoundaryEdge &); // interdit l'affectation par copie 
};


// La classe Triangle (modélisation des triangles)
class Triangle: public Label {
  Vertex *vertices[3]; // tableau de trois pointeurs de type Vertex
public:
  R area;
  Triangle(){}; // constructeur par défaut vide
  Vertex & operator[](int i) const {
    ASSERTION(i>=0 && i <3);
    return *vertices[i];} // évaluation des pointeurs 
  void set(Vertex * v0,int i0,int i1,int i2,int r) {
    vertices[0]=v0+i0; vertices[1]=v0+i1; vertices[2]=v0+i2; 
    R2 AB(*vertices[0],*vertices[1]);
    R2 AC(*vertices[0],*vertices[2]);
    area =abs( (AB^AC)*0.5);
    lab=r;
    if (area <0)  cout << " area = " << area <<  "\n" << *vertices[0] << endl  
                        << *vertices[1]  << endl << *vertices[2] << endl;
    ASSERTION(area>=0); }
  
  R2 Edge(int i) const {
    ASSERTION(i>=0 && i <3);
    return R2(*vertices[(i+1)%3],*vertices[(i+2)%3]);}// l'arête opposée au sommet i
  R2 H(int i) const { ASSERTION(i>=0 && i <3);
  R2 E=Edge(i);return E.perp()/(2*area);} // la hauteur
  R lenEdge(int i) const {
    ASSERTION(i>=0 && i <3);
    R2 E=Edge(i);return sqrt(E.x * E.x + E.y * E.y);}
private:
  Triangle(const Triangle &);       // interdit la construction par copie
  void operator=(const Triangle &); // interdit l'affectation  par copie
  
};

inline  int newt(int const  * const old2new,int t,int a)  
 {
 //  routine pour donner le numero nouveau triangle correspondant
 //  a l'arete a du traingle t apres decoupage
 //   Le tableau old2new[3*t+a] contient le numero du nouveau triangle
 //   ou -1 si triangle "nouveau"
 //   ou -2 si arete deja decoupe
 //   ou -3 si nouveau arete
   int tt;
   while ( (tt=old2new[3*t+a])>=0)  t=tt;
   return t;
 }


// La classe Mesh (modï¿½lisation du maillage)
class Mesh { public:
  int nt,nv,neb;
  int ntx,nvx,nebx;  // taille d'alloocation des tableau
  R area;
  Vertex *vertices;
  Triangle *triangles;
  BoundaryEdge  *bedges;  
  int * adj;
  Triangle & operator[](int i) const {return triangles[CheckT(i)];}
  Vertex & operator()(int i) const {return vertices[CheckV(i)];}
  Mesh(const char * filename,double fac=1); // lecture du fichier ".msh"
  Mesh(const std::string& filePath,double fac=1); // lecture du fichier ".msh"
  int operator()(const Triangle & t) const {return CheckT(&t - triangles);}
  int operator()(const Triangle * t) const {return CheckT(t - triangles);}
  int operator()(const Vertex & v) const {return CheckV(&v - vertices);}
  int operator()(const Vertex * v) const{return CheckT(v - vertices);}
  int operator()(int it,int j) const {return (*this)(triangles[it][j]);}
  //  vï¿½rification des depassements de tableau 
  int CheckV(int i) const { ASSERTION(i>=0 && i < nv); return i;} 
  int CheckT(int i) const { ASSERTION(i>=0 && i < nt); return i;}
    
  void BuildAdj();  //  construit le tableau :  adj[3*t+a] = 3*t'+a'  
  // ou   le arete a du triangle t est l'arete a' du triangle t'
  // ou 3*t'+a' = -1 si l'arete est sur le bord
  //  remarque: on pourait mettre -2-numero de l'arete du bord
  // pour envite une boucle stupide  dans decoupe1E
  
  void  decoupe(int * k,int n); // decoupe le tableau k[i] i=0, n-1 des k[i]%3 aretes du triangle k[i]/3

  void exportMesh(const std::string& fileName) const;

private:
 bool  decoupe(int * old2new,int t,int a);  //  decoupe les adjacent a l'arete a du triangle t 
 void  decoupe1T(int * old2new,int t,int a,int m); // decoupe le triangle t / l'arete a en 2 
 Label  decoupe1E(int k,int i0,int i1,int m); // decoupe  l'arete frontiere  en 2 et retourne le lablel de l'arete
   
  Mesh(const Mesh &);              // interdit la construction par copie
  void operator=(const Mesh &);    // interdit l'affectation  par copie
};
