
class LexicoEdge { public:
  int t, e; //  numero du triangle, numero de larete dans le triangle
  int i0,i1; // les 2 sommets de l'arete tel que i0 <= i1
  LexicoEdge() :t(0),e(0),i0(0),i1(0) {}
  LexicoEdge(int tt,int ee,const Mesh & Th ) : t(tt),e(ee),i0(Th(tt,(ee+1)%3)),i1(Th(tt,(ee+2)%3)) 
  { if (i0 > i1) {int i=i0;i0=i1;i1=i;}}
};

inline   bool operator<(const LexicoEdge & a,const LexicoEdge & b)
{ return a.i0 == b.i0 ? (a.i1 < b.i1) :  (a.i0 < b.i0);}

inline   bool operator==(const LexicoEdge & a,const LexicoEdge & b)
{ //cout <<  a.i0 << " " << a.i1 << " == " << b.i0 << " " << b.i1 << endl;
   return (a.i1 == b.i1) && (a.i0 == b.i0);}




