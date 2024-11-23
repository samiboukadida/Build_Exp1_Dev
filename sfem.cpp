//----for handling paths--------
#include <filesystem>
#include <algorithm>
#include <cstring>
//-------------
#include <fstream>
#include <iostream>
#include <map>
#include <functional>
#include "sfem.hpp"
#include "RNM.hpp"

#include "LexicoEdge.hpp"
#include "HeapSort.hpp"
/*struct PairLess : binary_function<pair<int,int>,pair<int,int>, bool>
{ 
  typedef pair<int,int>  Key;
  bool operator()(const Key& x, const Key& y) const { return x.first == y.first ? x.second < y.second : x.first < y.first ;} 
};*/
struct PairLess {
    typedef std::pair<int, int> Key;
    bool operator()(const Key& x, const Key& y) const {
        return x.first == y.first ? x.second < y.second : x.first < y.first;
    }
};


// Le constructeur de la classe Mesh
Mesh::Mesh(const char*  filePath,double fac)
{ // lecture du maillage
  int i,i0,i1,i2,ir;
  //ifstream f(filename);
  //if(!f) {cerr << "Mesh::Mesh Erreur a l'ouverture - fichier " << filename<<endl;exit(1);}
  //cout << "Lecture du fichier  \"" <<filename<<"\""<<  endl;
  // Check if the filename is an absolute path----------------------------------
  // Convert the Windows-style path by replacing backslashes with slashes (if necessary)
  // Check if the file path is an absolute Windows path
  std::ifstream f;

  // Check if the path is an absolute Windows path (e.g., "C:\...")
  if (std::strlen(filePath) > 2 && filePath[1] == ':' && (filePath[2] == '\\' || filePath[2] == '/')) {
      // Absolute Windows path, directly use the provided path
    f.open(filePath);
  } else {
      // Otherwise, treat it as a relative path, prepend the current working directory
      std::filesystem::path fullPath = std::filesystem::current_path() / filePath;
      f.open(fullPath.string());

      // Update filePath to the full path for consistent output
      filePath = fullPath.string().c_str();
  }

  // Error handling if the file can't be opened
  if (!f) {
      std::cerr << "Mesh::Mesh Erreur à l'ouverture - fichier " << filePath << std::endl;
      exit(1);
  }

  // Log that the mesh file is being read
  std::cout << "Lecture du fichier \"" << filePath << "\"" << std::endl;

  // The rest of the Mesh constructor logic remains unchanged...----------------
  f >> nv >> nt >> neb ;
  cout << " Nb de sommets " << nv << " " << " Nb de triangles " 
       << nt << " Nb d'aretes frontiere " << neb << endl;
  nvx = Max((int) (nv*fac),nv);
  ntx = Max((int)( nt*fac),nt);
  nebx = Max((int) (neb*fac),neb);
  
  assert(f.good() && nt && nv) ;
  triangles = new Triangle [ntx];
  vertices = new Vertex[nvx];
  bedges = new BoundaryEdge[nebx];
  adj = new int[ntx*3];
  area=0;
  assert(triangles && vertices);
  for (i=0;i<nv;i++)
    f >> vertices[i] ,  assert(f.good());
  for (i=0;i<nt;i++) { 
    f >> i0 >> i1 >> i2 >> ir;
    assert(f.good() && i0>0 && i0<=nv && i1>0 && i1<=nv && i2>0 && i2<=nv);
    triangles[i].set(vertices,i0-1,i1-1,i2-1,ir); 
    area += triangles[i].area;
  }
  
  for (i=0;i<neb;i++) { 
    f >> i0 >> i1  >> ir;
    assert(f.good() && i0>0 && i0<=nv && i1>0 && i1<=nv );
    bedges[i].set(vertices,i0-1,i1-1,ir);
  }
  BuildAdj();
  cout << " Fin lecture : aire du maillage = " << area <<endl;  
} 
void Mesh::BuildAdj()
{
// construction du tableau des adjacence :
//  remarque:
//  tableau des adjacents  l'arete  a  triangle t ==  l'arete a'  de t'
//   adj[k]=k' ou    k=3*t+a et k'= 3*t'+a'
// si k' < -1 =>   l'arete est l'arete du bord
// remarque on peut mettre dans k" -2-numero de l'arete de bord pour la retourve
// rapidement  et evite une boucle dans decoup1E (bonus)
// ------------------
 int nt3 = nt*3;
  for (int i=0;i<nt3;i++)   adj[i]=-1; // pas d'adjacence par defaut
  const Mesh &Th(*this);   
  // ici mettre la construction du tableau adj (F. hecht)
  LexicoEdge *le = new LexicoEdge[nt3];

  for (int k=0,kk=0;k<Th.nt;k++)
    for (int i=0;i<3;i++)
      le[kk++]=LexicoEdge(k,i,Th);

  HeapSort(le,nt3);
  int nbarete = nt3;
  for (int i=1;i<nt3;i++)
    if ( le[i-1]==le[i] ){
        nbarete--; // arete commune
        int ep=le[i-1].e,kp=le[i-1].t,es=le[i].e,ks=le[i].t;
        int kkp=3*kp+ep;
        int kks=3*ks+es;
        adj[kkp]=kks;adj[kks]=kkp;
    // cout<< "adj["<<kkp<<"]="<<kks<<endl;
    }

  cout << " Nb Arete = " << nbarete << " == " << Th.nt+Th.nv-1 << " si domaine simple " <<  endl;
}

void  Mesh::decoupe(int * k,int n)
{ 
  KN<int> old2new(ntx*3);
  old2new = -1;
   int m=0;
   for (int i=0;i<n;i++)
     if (decoupe(old2new,k[i]/3,k[i]%3)) {m++;cout<<" Nombre d'arete decoupe=  "<<m<<endl;}
   cout << " Nombre d'arete decoupe  " << m << " sur  " << n << endl;
   cout << " Nb de sommets " << nv << " " << " Nb de triangles " 
       << nt << " Nb d'aretes frontiere " << neb << endl;
     
    BuildAdj();  
}
 bool  Mesh::decoupe(int * old2new,int t,int a)
  {
    
    int k=3*t+a;
    int kk= adj[3*t+a]; 
    int tt=kk/3, aa=kk%3; // le triangle adjacent
    if ( old2new[k] <-1) return false; // deja 
    if ( kk>=0 && old2new[kk] <-1) return false; //  deja decoupe    
    assert(nv+1<nvx);
    int m = nv++;
    //__

    int ttt=t;
   
  if (old2new[k]>=0)
   
    {  ttt=old2new[k];cout<<"old2new["<<k<<"]="<< old2new[k]<<endl;}

  else if( old2new[k]==-1) ttt=t;
    const  Triangle &  K(triangles[ttt]);
    const Vertex & A(K[(a+1)%3]), &B(K[(a+2)%3]);

    (R2 &) vertices[m] = ( (R2) A + (R2) B) *.5;     
    vertices[m].lab=0;  // par defaut arete interne  
    decoupe1T(old2new,newt(old2new,t,a),a,m);  // decoupe du nouveau triangle t par l'arete a
    if (kk>=0)  // arete interne 
      decoupe1T(old2new,newt(old2new,tt,aa),aa,m);  // decoupe du nouveau triangle tt par l'arete aa
    else 
      {
      // sur la fontiere on decoup l'arete fontiere 
      //  il est possible de stocker kk= -2- numero de l'arete frontiere pour eviter une boucle (bonus)
      // sinon mettre simple kk=-1;
     // cout << t << " " << a << " " << kk << endl;
      (Label &) vertices[m] =decoupe1E(-kk-2, operator()(A),operator()(B),m); // decoupe de l'arete fontiere
      
      }
    return true;
  }

Label  Mesh::decoupe1E(int k,int i0,int i1,int m)
{ 
  // si k >=0 => numero de l'arete a decoupe sinon on recheche dans toutes les aretes
  // pour trouver l'arete 
  // => algo en neb^2 , peut faire beaucoup mieux 
  // il suffit de parcourir les aretes de bord ayant un meme sommet .
  // ou trie les aretes par numero de sommet et faire une recherche dichotomique
  // FH> etudiant : bonus pour qui le fait 
  Vertex *  v0  = vertices + i0;
  Vertex *  v1  = vertices + i1;
  int istart=0,iend=neb; // borne de la boucle 
  if (k>0) istart=k,iend=k+1; //  pour evite la boucle
  for (int i=istart;i<iend;i++)
    if ( bedges[i].in(v0) &&  bedges[i].in(v1) )
      {
         assert(neb+1 < nebx);
         int lab =  bedges[i].lab;
         bedges[i].set(vertices,i0,m,lab);
         bedges[neb++].set(vertices,m,i1,lab);  
	 cout<<"i   bedges[i]="<<i<<"    "<<bedges[i]<<endl;       
         return  bedges[i];
      }
  cout << " bug arete non trouver! " << istart << " " << iend << " " << i0 << " " << i1 << " " << k << endl;
  assert(0); // on n'a pas trouver l'arete 
}

void  Mesh::decoupe1T(int * old2new,int t,int a,int m)
   {    
     int i0=a,i1=(a+1)%3, i2=(a+2)%3;
     int num[3], num1[3]; // les numero des nouveau triangles
     assert(nt+1 < ntx); 
     int t1 = nt++;
     Triangle & K(triangles[t]);
     int lab=K.lab;
     Triangle & K1(triangles[t1]);     
      num1[0]=num[0]=(*this)(K[0]); // les trois numero de sommets
      num1[1]=num[1]=(*this)(K[1]);
      num1[2]=num[2]=(*this)(K[2]);
      /* 
           i0
          /|\
     t   / | \  t1  cette arete t1, i1 change 
        /  |  \
     i1/___|___\ i2
           m
          a       
      */
    //  num : numero des 3 sommet du nouveau triangle K
    //  num1 : numero des 3 sommet du nouveau triangle K1
      
      //_________    ici  changer num et num1  et old2new  (F. Hecht)

        num[i2]=m  ; num1[i1]=m;
        old2new[3*t+i0]=-2;
        old2new[3*t+i1]=t1;
      K.set(vertices,num[0],num[1],num[2],lab);
   //   cout << " -- "<< endl;
        K1.set(vertices,num1[0],num1[1],num1[2],lab);
      
   }

void Mesh::exportMesh(const std::string& fileName) const {
    std::ofstream outFile(fileName);
    if (!outFile) {
        std::cerr << "Erreur : impossible de créer le fichier " << fileName << std::endl;
        return;
    }

    // Écrire les entêtes (nv, nt, neb)
    outFile << nv << " " << nt << " " << neb << "\n";

    // Écrire les sommets
    for (int i = 0; i < nv; ++i) {
        outFile << vertices[i].x << " " << vertices[i].y << " " << vertices[i].lab << "\n";
    }

    // Écrire les triangles
    for (int i = 0; i < nt; ++i) {
        const Triangle& tri = triangles[i];
        outFile << (*this)(tri[0]) + 1 << " " // Conversion base 0 -> base 1
                << (*this)(tri[1]) + 1 << " "
                << (*this)(tri[2]) + 1 << " "
                << tri.lab << "\n";
    }

    // Écrire les arêtes frontières
    for (int i = 0; i < neb; ++i) {
        const BoundaryEdge& edge = bedges[i];
        outFile << (*this)(edge[0]) + 1 << " " // Conversion base 0 -> base 1
                << (*this)(edge[1]) + 1 << " "
                << edge.lab << "\n";
    }

    std::cout << "Fichier de maillage exporté : " << fileName << std::endl;
}
