template<class R>
inline KN_<R>&  KN_<R>::operator oper (const_R a)  {
    R * l(v);    
    for (long i=0;i<n;i++,l += step) 
      *l oper a;
    return *this;
  }
    
template<class R> 
inline    KNM_<R> & KNM_<R>::operator oper (const_R a)    
{ 
  if(IsVector1() ) 
        KN_<R>::operator oper (a);
  else {  
          KN_<R>  lj(operator()('.',0)); //  (.,.,O)
          for (long  j=0;j<M();++j,++lj) 
             lj oper a;}       
  return *this;
}

template<class R> 
inline    KNMK_<R> & KNMK_<R>::operator oper (const_R a)    
{ 
  if(IsVector1() ) 
        KN_<R>::operator oper (a);
  else {  
          KNM_<R>  lj(operator()('.','.',0)); //  (.,.,O)
          long j=K();
           while(j--)
            {lj oper a;++lj;}
       }
  return *this;
}

template<class R>
inline KN_<R>&  KN_<R>::operator oper (const KN_<const_R> & u)   {
    K_throwassert(u.n == n);
    R * l(v);
    const R *r(u);    
    for (long i=0;i<n;i++,l += step, r += u.step) *l oper *r;
    return *this;
  }
  
template<class R> 
inline    KNM_<R> & KNM_<R>::operator oper (const KNM_<const_R> & u)    
{ 
  if(IsVector1() && u.IsVector1() ) 
        KN_<R>::operator oper(u); // modif FH jan 2004
  else {  
          KN_<R>  lj(operator()('.',0)); //  (.,O)
          KN_<const_R>  uj(u('.',0));
          long  j=M();
          while ( j--)
            { lj oper uj;++lj;++uj;} 
        }      
  return *this;
}


template<class R> 
inline  KNMK_<R> & KNMK_<R>::operator oper (const KNMK_<const_R> & u)    
{ 
  if(IsVector1() && u.IsVector1() ) 
        KN_<R>::operator oper(u); // modif FH 2004
  else {  
          K_throwassert( K() == u.K());
          KNM_<R>  lj(operator()('.','.',0)); //  (.,O)
          KNM_<const_R>  uj(u('.','.',0));
          long j=K();
          while (j--)
           { lj oper uj;++lj;++uj;}
       }
  return *this;
}

//   template<class I> KN_& operator =  (const KN_ITAB<const KN_<R>,I> & u);

template<class T,class I>
KN_ITAB<T,I> & KN_ITAB<T,I>::operator oper (const T & t)
{
  int n= index.N(); 
  K_throwassert( n <= t.n  );
  
  typename T::K  * l(t);
  int step = t.step;
  for (long i=0;i<n;i++,l += step) 
    v[index[i]] oper *l;
  return *this;
}

template<class R>
template<class I,class T>
inline KN_<R>&  KN_<R>::operator oper (const KN_ITAB<T,I> & ui)   {
   const KN<R> & u(ui.v);
    I & itab(ui.index);
    R  * l(v);
    for (long i=0;i<n;i++,l += step)
      *l oper u[itab[i]];
    return *this;
  }
#undef oper
