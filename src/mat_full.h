/*! \addtogroup LA Linear Algebra Package with Sparse Implementations */
/*! @{*/
/*! @file mat_full.h
  @brief Provide functionality for full matrix and vector interactions.
*/

#ifndef _MAT_FULL_H
#define _MAT_FULL_H
#include <vector>
#include <stdexcept>
#include <linear_algebra.decl>
#include <matrix.decl>
#include <refcount.decl>
#include <sparse_vector.decl>
#include <mat_sym_full.decl>
#include <mat_asym_full.decl>
#include <mat_sym_sparse.decl>

namespace linear_algebra
{

//! Class of full column matrices.
/*! 
  'column' specifies that the values are ordered for speedy retrieval per column.
*/

  // Constructors


  //! Constructor of n columns and m rows with uninitialised values.
  template <class TN> mat_full<TN>::mat_full(long n, long m) {
    matrix<TN>::_rows=m; matrix<TN>::_cols=n;
    refvector<refvector<TN> > r(n);
    _size= n*m;
    for(long i=0; i<n; i++)
      {
	refvector<TN> s(m);
	r[i]=s;
      }
    _scols=r;
  }

  //! Constructor of n columns and m rows with initialisation of values. You had better be sure there is no buffer overrun.
  template <class TN> mat_full<TN>::mat_full(long n, long m, const TN *vals[])
  /*!
    @param n Columns of matrix.
    @param m Rows of matrix.
    @param vals array of columns. Each column is a full array.
    
    Each value has an associated row and column index in srows and scols. 
    Therefore they all have the same size.
  */
  {
    _rows=m;
    _cols=n;
    _size=m*n;
    _scols.resize(n);
    for(long i=0;i<_cols;i++) {
      _scols[i].resize(m);
      for(long j=0;j<_rows;j++)
	_scols[i][j]=vals[i][j];
    }
  }

  //! Constructor of n columns and m rows with initialisation of values.
  template <class TN> mat_full<TN>::mat_full(long n, long m, const vector<vector<TN> >& vals)
    /*!
      @param n Dimension of matrix.
      @param scols Vector of column indices.
      @param vals Vector of columns. Each column is a full vector.
      
      Each value has an associated row and column index in srows and scols. 
      Therefore they all have the same size.
    */
    {
      matrix<TN>::_rows=m;
      matrix<TN>::_cols=n;
      _size=m*n;
      if(matrix<TN>::_cols!=vals.size()) 
	{
	  throw domain_error("values and indices don't match up in   mat_full(long n, long m, vector<long> scols, vector<vector<TN> >& vals) ");
	  exit(1);
	}
      _scols.resize(_cols);
      for(long i=0; i<_cols;i++)
	_scols[i]=vals[i];
    }
  
  
  //! Constructor of n columns and m rows with initialisation of values.
  template <class TN> mat_full<TN>::mat_full(long n, long m, const refvector<vector<TN> >& vals)
    /*!
      @param n Dimension of matrix.
      @param scols Vector of column indices.
      @param vals Reference counted vector of columns. Each column is a full vector.
      
      Each value has an associated row and column index in srows and scols. 
      Therefore they all have the same size.
    */
    {
      _rows=m;
      _cols=n;
      _size=m*n;
      if(_cols!=vals->size()) 
	{
	  throw domain_error("values and indices don't match up in   mat_full(long n, long m, vector<long> scols, vector<vector<TN> >& vals) ");
	  exit(1);
	}
      _scols=vals;
    }
  
  //! Constructor of n columns and m rows with initialisation of values.
  template <class TN> mat_full<TN>::mat_full(long n, long m, const refvector<refvector<TN> >& vals)
    /*!
      @param n Dimension of matrix.
      @param scols Reference counted vector of column indices.
      @param vals Reference counted vector of columns. Each column is a full vector.
      
      Each value has an associated row and column index in srows and scols. 
      Therefore they all have the same size.
    */
    {
      _rows=m;
      _cols=n;
      _size=m*n;
      if(_cols!=vals.size()) 
	{
	  throw domain_error("values and indices don't match up in   mat_full(long n, long m, vector<long> scols, vector<vector<TN> >& vals) ");
	  exit(1);
	}
      _scols=vals;
    }
  
  //! Constructor of n columns and m rows with initialisation of values.
  template <class TN> mat_full<TN>::mat_full(long n, long m, const vector<refvector<TN> >& vals)
    /*!
      @param n Dimension of matrix.
      @param scols Reference counted Vector of column indices.
      @param vals Vector of columns. Each column is a full refvector.
      
      Each value has an associated row and column index in srows and scols. 
      Therefore they all have the same size.
    */
    {
      _rows=m;
      _cols=n;
      _size=m*n;
      if(_cols!=vals.size()) 
	{
	  throw domain_error("values and indices don't match up in   mat_full(long n, long m, vector<long> scols, vector<refvector<TN> >& vals) ");
	  exit(1);
	}
      _scols=vals;
    }
  
  template <class TN> mat_full<TN>::mat_full() {_size=0; _rows=0; _cols=0; }
  
  //! Constructor to read a stored matrix from a stream. Aimed at retrieving information from a binary save with mat_full::operator>>().
  template <class TN> mat_full<TN>::mat_full(istream& IN)
  {
    long rs,i,j,k;
#ifdef DEBUG
    const char* string("Full matrix:\n\x0");
    char a;
    i=0;
    IN.read((char*) &a,sizeof(char));
    while(string[i]==a && string[i]!='\n') 
      {
	i++;
	IN.read((char*) &a,sizeof(char));
      }
    if(string[i]!='\n')
      throw domain_error("This is not a full matrix");
#else
    IN.seekg(13,ios::cur);
#endif
    
    
    IN.read((char*) &rs,sizeof(long)/sizeof(char));
    IN.read((char*) &_cols,sizeof(long)/sizeof(char));
    IN.read((char*) &_rows,sizeof(long)/sizeof(char));
#ifdef DEBUG
    if (rs!=_rows*_cols*sizeof(TN)+2*sizeof(long))
      throw domain_error("Dimension of Matrix and saved data do not agree in mat_full(istream& IN)");
#endif
    
    mat_full<TN> M(_cols,_rows);
    
    refvector<char> buffer(_rows*sizeof(TN)/sizeof(char));
    
    TN b;
    char *pa=(char *) &b;
    
    for(j=0;j<_cols;j++)
      {
	IN.read(&buffer[0],_rows*sizeof(TN)/sizeof(char));
	for(i=0;i<_rows;i++)
	  {
	    for(k=0;k<(long) (sizeof(TN)/sizeof(char));k++)
	      pa[k]=buffer[i*sizeof(TN)/sizeof(char)+k];
	    M[j][i]=b;
	  }
      }
    *this=M;
  }
  
  template <class TN> mat_full<TN>::~mat_full() {}
  
  // Operators *********************************************************************************************************
  
  //! Assignment operator for non const objects.
  template <class TN> mat_full<TN>& mat_full<TN>::operator=(mat_full<TN>& a)
  {
    _size=a._size;
    _rows=a._rows;
    _cols=a._cols;
    _scols=a._scols;
    return *this;
  }

  //! Assignment operator for non const objects.
  template <class TN> mat_full<TN>& mat_full<TN>::operator=(const mat_full<TN>& a)
  {
    _size=a._size;
    _rows=a._rows;
    _cols=a._cols;
    _scols.copy(a._scols);
    return *this;
  }


  //! Compares two full matrices.
  template <class TN> bool mat_full<TN>::operator==(const mat_full<TN>& b) const
    {
      if(_size==b._size &&
	 _rows==b._rows &&
	 _cols==b._cols &&
	 _scols==b._scols) return true;
      else return false;
    }
  
  template <class TN> bool mat_full<TN>::same(const mat_full<TN>& b) const
    {
      if ((*this)==b && b._scols.same(_scols)) return true;
      else return false;
    }
  
  //! Return a column through direct access. No integrity checked.
  template <class TN> refvector<TN>& mat_full<TN>::operator[](const long x) {
#ifdef DEBUG
    if(x<0 || x- _scols.size()>=0)
      {
	char v[100];
	snprintf(v,100,"mat_full::operator[](const long x=%li): illegal access.",x);
	throw domain_error(v);
      }
#endif
    return _scols[x]; }

  //! Return a column through direct access. No integrity checked.
  template <class TN> const refvector<TN>& mat_full<TN>::operator[](const long x) const {
#ifdef DEBUG
    if(x<0 || x-_scols.size()>=0) throw domain_error("mat_full::operator[](const long x) const: illegal acces.");
#endif
    return _scols[x];}
  
  //! Return value at (x,y).
  template <class TN> const TN& mat_full<TN>::operator()(const long x,const long y) const
    {
#ifdef DEBUG
      if(x-_scols.size()>=0 || x<0)
	{
	  char v[100];
	  snprintf(v,100,"mat_full::operator()(const long x=%li,const long y=%li) const:Illegal x-access to vector",x,y);
	  throw domain_error(v);
	}
      if(y-_scols[x].size()>=0 || y<0)
	{
	  char v[100];
	  snprintf(v,100,"mat_full::operator()(const long x=%li,const long y=%li) const:Illegal y-access to vector",x,y);
	  throw domain_error(v);
	}
#endif
      return _scols[x][y];
    }
  //! Return value at (x,y).
  template <class TN> TN& mat_full<TN>::operator()(const long x,const long y)
    {
#ifdef DEBUG
      if(x>=_scols.size() || x<0) 
	{
	  char v[100];
	  snprintf(v,100,"mat_full::operator()(const long x=%li,const long y=%li) const:Illegal x-access to vector",x,y);
	  throw domain_error(v);
	}
      if(y>=_scols[x].size() || y<0) 
	{
	  char v[100];
	  snprintf(v,100,"mat_full::operator()(const long x=%li,const long y=%li) const:Illegal y-access to vector",x,y);
	  throw domain_error(v);
	}
#endif
      return _scols[x][y];
    }
  
  //! Multiply a full matrix by a sparse vector.@test
  template <class TN> refvector<TN>& mat_full<TN>::multiply(const sparse_vector<TN>& a, refvector<TN>& r) const
    {
#ifdef DEBUG
      if(_cols != a.dim())
	throw domain_error("Incompatible dimensions in refvector<TN> mat_full::multiply(const sparse_vector<TN>& a) const");
      if(_rows != r.dim())
	throw domain_error("Incompatible dimensions for output in refvector<TN> mat_full::multiply(const sparse_vector<TN>& a) const");
#endif
      long i;
      for(i=0;i<r.dim();i++) r[i]=(TN) 0;
      refvector<double> m(_rows);
      for (i=0; i<a.size(); i++) r+=_scols[a.index(i)].multiply(a[i],m);
      return r;
    }
  
  //! Multiply a full matrix by a sparse vector.@test
  template <class TN> refvector<TN> mat_full<TN>::operator*(const sparse_vector<TN>& a) const
    {
      refvector<TN> r(_rows);
      return multiply(a,r);

    }
  
  //! Multiply a full matrix by a number.@test
  template <class TN> mat_full<TN> mat_full<TN>::operator*(const TN& a) const
    {
      mat_full<TN> r(_cols,_rows);
      return multiply(a,r);
    }
  
  //! Multiply a full matrix by a number.@test
  template <class TN> mat_full<TN>& mat_full<TN>::multiply(const TN& a,mat_full<TN>& r) const
    {
#ifdef DEBUG
      if(r.cols()!=_cols || r.rows() != _rows)
	throw domain_error("mat_full::multiply(const TN& a, ...) called with incorrect output dimensions.");
#endif
      long i;
      for (i=0; i<_cols; i++) _scols[i].multiply(a,r[i]);
      return r;
    }
  
  //! Multiply by an object of TN into left side.
  template <class TN> mat_full<TN>& mat_full<TN>::operator*=(const TN& x)
    {
      for(long i=0; i<_cols; i++)
	{
	  _scols[i]*=x;
	}
      return *this;
    }
  //! Divide by an object of TN into left side. This is very slow and not advisable.
  template <class TN> mat_full<TN>& mat_full<TN>::operator/=(const TN& x)
    {
      long i=0;
      TN x_1=(TN) 1/x;
      while(i<_cols)
	{
	  _scols[i]*=x_1;
	  i++;
	}
      return *this;
    }
  
  //! Multiply an asymmetric, full matrix from the left into this full matrix.
  template <class TN> mat_full<TN>& mat_full<TN>::operator*=(const mat_asym_full<TN>& b)
    {
#ifdef DEBUG
      if(b.cols() != _rows)
	throw domain_error("Incompatible dimensions in mat_full mat_full::operator*=(const mat_asym_full& b)");
      try {
#endif
	long i;
	for(i=0;i<_cols;i++)
	  (*this)[i]=b*(*this)[i];
	_rows=b.rows();
	return (*this);
#ifdef DEBUG
      } catch (domain_error e) {
	cerr << e.what() << endl;
	throw domain_error("called from mat_full::operator*=(const mat_asym_full<TN>& b)");
      }
#endif
    }

  //! Multiply a symmetric, full matrix from the left into this full matrix.
  template <class TN> mat_full<TN>& mat_full<TN>::operator*=(const mat_sym_full<TN>& b)
    {
#ifdef DEBUG
      if(b.cols() != _rows)
	throw domain_error("Incompatible dimensions in mat_full mat_full::operator*=(const mat_sym_full& b)");
      try {
#endif
	long i;
	for(i=0;i<_cols;i++)
	  (*this)[i]=b*(*this)[i];
	_rows=b.rows();
	return (*this);
#ifdef DEBUG
      } catch (domain_error e) {
	cerr << e.what() << endl;
	throw domain_error("called from mat_full::operator*=(const mat_sym_full<TN>& b)");
      }
#endif
    }
  
  //! Multiply two full matrices.
  template <class TN> mat_full<TN>& mat_full<TN>::multiply(const  mat_full<TN>& b, mat_full<TN>& r) const
    {
#ifdef DEBUG
      if( _cols != b._rows)
	{
	  throw domain_error( "Incompatible dimensions in mat_full mat_full::multiply(const mat_full& a, const mat_full& b)");
	  exit(1);
	}
      if( r.cols() != b.cols() || r.rows()!=_rows)
	{
	  throw domain_error( "Incompatible dimensions for output in mat_full mat_full::multiply(const mat_full& a, const mat_full& b)");
	}
      if(b.same(r) || same(r))
	throw domain_error("mat_full::multiply(mat_full, ..) output may not coincide with the arguments.");
      try {
#endif
	long i;
	for(i=0; i<b._cols;i++)
	  multiply(b._scols[i],r[i]);
	return(r);
#ifdef DEBUG
      } catch(domain_error e) { 
	cerr << e.what() << endl;
	throw domain_error(" called from mat_full<TN> mat_full::operator*(const  mat_full<TN>& b) const");
      }
#endif
    }
  
  //! Multiply two full matrices.
  template <class TN> mat_full<TN> mat_full<TN>::operator*(const mat_full<TN>& b) const
    {
      mat_full<TN> r(b.cols(),_rows);
      return multiply(b,r);
    }
  
  //! Multiply the transpose of a full matrix with a full matrix.
  template <class TN> mat_full<TN>& mat_full<TN>::multiply_t_a(const  mat_full<TN>& b, mat_full<TN>& r) const
  {
#ifdef DEBUG
    if( _rows != b._rows)
      {
	throw domain_error( "Incompatible dimensions in mat_full mat_full::multiply_t_a(const mat_full& a, const mat_full& b)");
	exit(1);
      }
    if( r.cols() != b.cols() || r.rows()!=_cols)
      {
	throw domain_error( "Incompatible dimensions for output in mat_full mat_full::multiply_t_a(const mat_full& a, const mat_full& b)");
      }
    if(b.same(r) || same(r))
      throw domain_error("mat_full::multiply_t_a(mat_full, ..) output may not coincide with the arguments.");
    try {
#endif
      long i,j;
      for(i=0; i<r.cols();i++)
	for(j=0;j<r.rows();j++)
	  r[i][j]=(*this)[j]*b[i];
      return(r);
#ifdef DEBUG
    } catch(domain_error e) { 
      cerr << e.what() << endl;
      throw domain_error(" called from mat_full<TN> mat_full::multiply_t_a*(const  mat_full<TN>& b) const");
    }
#endif
  }
  
  //! Multiply the transpose of a full matrix with a full matrix.
  template <class TN> mat_full<TN> mat_full<TN>::multiply_t_a(const mat_full<TN>& b) const
  {
    mat_full<TN> r(b.cols(),_cols);
    return multiply_t_a(b,r);
  }
  
  //! Multiply full matrix with the transpose of a full matrix.
  template <class TN> mat_full<TN>& mat_full<TN>::multiply_a_t(const  mat_full<TN>& b, mat_full<TN>& r) const
  {
#ifdef DEBUG
    if( _cols != b._cols)
      {
	throw domain_error( "Incompatible dimensions in mat_full mat_full::multiply_a_t(const mat_full& a, const mat_full& b)");
	exit(1);
      }
    if( r.cols() != b.rows() || r.rows()!=_rows)
      {
	throw domain_error( "Incompatible dimensions for output in mat_full mat_full::multiply_a_t(const mat_full& a, const mat_full& b)");
      }
    if(b.same(r) || same(r))
      throw domain_error("mat_full::multiply_a_t(mat_full, ..) output may not coincide with the arguments.");
    try {
#endif
      long i,j,k;
      r.zero();
      for(k=0;k<b.cols();k++)
	for(i=0; i<r.cols();i++)
	  for(j=0;j<r.rows();j++)
	    r[i][j]+=(*this)[k][j]*b[k][i];
      return(r);
#ifdef DEBUG
    } catch(domain_error e) { 
      cerr << e.what() << endl;
      throw domain_error(" called from mat_full<TN> mat_full::multiply_a_t(const  mat_full<TN>& b, ) const");
    }
#endif
  }
  
  //! Multiply full matrix with the transpose of a full matrix.
  template <class TN> mat_full<TN> mat_full<TN>::multiply_a_t(const mat_full<TN>& b) const
  {
    mat_full<TN> r(b.rows(),_rows);
    return multiply_a_t(b,r);
  }
  
  //! Multiply a full matrix by a symmetric full matrix.
  template <class TN> mat_full<TN>& mat_full<TN>::multiply(const mat_sym_full<TN>& b,mat_full<TN>& r) const
    {
#ifdef DEBUG
      if(b.rows() != _cols)
	throw domain_error("Incompatible dimensions in mat_full mat_full::multiply(const mat_sym_full<TN>& b) const");
      if(r.rows() != _rows || r.cols()!=b.cols())
	throw domain_error("Incompatible dimensions for output in mat_full mat_full::multiply(const mat_sym_full<TN>& b) const");
      if( same(r))
	throw domain_error("mat_full::multiply(mat_sym_full, ...) output must be != input.");
#endif
      long i, j;
      // Upper triangle
      refvector<TN> a(_rows);
      for (i=0; i< b.cols(); i++)
	{
	  for(j=0;j<a.dim();j++) r[i][j]=(TN) 0;
	  for(j=0; j<=i; j++)
	    r[i]+=_scols[j].multiply(b[i*(i+1)/2+j],a);
	}
      // Lower triangle
      for(i=0; i<b.rows(); i++)
	for(j=0; j<i; j++) r[j]+=_scols[i].multiply(b[i*(i+1)/2+j],a);
      return r;
    }

  //! Multiply a full matrix by a symmetric full matrix.
  template <class TN> mat_full<TN> mat_full<TN>::operator*(const mat_sym_full<TN>& b) const
    {
      mat_full<TN> r(b.cols(),_rows);
      return multiply(b,r);
    }
  
  //! Multiply a full matrix by an antisymmetric full matrix.
  template <class TN> mat_full<TN>& mat_full<TN>::multiply(const mat_asym_full<TN>& b,mat_full<TN>& r) const
    {
#ifdef DEBUG
      if(b.rows() != _cols)
	throw domain_error("Incompatible dimensions in mat_full mat_full::multiply(const mat_sym_full<TN>& b) const");
      if(r.rows() != _rows || r.cols()!=b.cols())
	throw domain_error("Incompatible dimensions for output in mat_full mat_full::multiply(const mat_sym_full<TN>& b) const");
      if( same(r))
	throw domain_error("mat_full::multiply(mat_asym_full, ...) output must be != input.");
#endif
      long i, j;
      
      // Upper triangle
      refvector<TN> a(_rows);
      for(j=0;j<r.cols();j++)
	for(i=0;i<r.rows();i++) r[j][i]=(TN) 0;
      for(i=1;i<b.cols();i++)
	for(j=0;j<i;j++)
	  r[i]+=_scols[j].multiply(b[i*(i+1)/2+j],a);
      
      
      // Lower triangle
      for(j=1;j<b.cols();j++)
	for(i=0;i<j;i++)
	  r[i]-=_scols[j].multiply(b[j*(j+1)/2+i],a);
      
      return r;
    }

  //! Multiply a full matrix by an antisymmetric full matrix.
  template <class TN> mat_full<TN> mat_full<TN>::operator*(const mat_asym_full<TN>& b) const
    {
      mat_full<TN> r(b.cols(),_rows);
      return multiply(b,r);
    }


  //! Multiply a full matrix by a symmetric sparse matrix.
  template <class TN> mat_full<TN>& mat_full<TN>::multiply(const mat_sym_sparse<TN>& b, mat_full<TN>& r) const
    {
#ifdef DEBUG
      if(b.rows() != _cols)
	throw domain_error("Incompatible dimensions in mat_full mat_full::operator*(const mat_sym_sparse<TN>& b) const");
      if(b.cols() != r.cols() || _rows!=r.rows())
	throw domain_error("Incompatible dimensions for output in mat_full mat_full::operator*(const mat_sym_sparse<TN>& b) const");
#endif
      long i,j;
      // Upper triangle
      refvector<TN> a(_rows);
      for (i=0; i< b.size(); i++)
	for(j=0; j<b[i].size(); j++)
	  r[b.index(i)]+=_scols[b[i].index(j)].multiply(b[i][j],a);
      // Lower triangle
      for(i=0; i<b.size(); i++)
	for(j=0; j<i; j++) r[b[i]._index[j]]+=_scols[b._scols[i]].multiply(b[i][j],a);
      
      return r;
    }
  
  //! Multiply a full matrix by a symmetric sparse matrix.
  template <class TN> mat_full<TN> mat_full<TN>::operator*(const mat_sym_sparse<TN>& b) const
    {
      mat_full<TN> r(b.cols(),_rows);
      return multiply(b,r);
    }
  

  //! Multiply a full matrix by a full vector.
  template <class TN> template<class TN2>
    refvector<TN>& mat_full<TN>::multiply(const refvector<TN2>& a, refvector<TN>& r) const
    /*!
      A little care needs to be taken here because any vector for which TN*TN2 is defined this product is also defined.
    */
    {
#ifdef DEBUG
      if(a.size() != _cols)
	throw domain_error("Incompatible dimensions in refvector<TN> mat_full::multiply(const refvector<TN>& a)" );
      if(r.dim() != _rows)
	throw domain_error("Incompatible dimensions for output in refvector<TN> mat_full::multiply(const refvector<TN>& a)" );
      try {
#endif
	refvector<double> b(_rows);
	if(_cols>0) _scols[0].multiply(a[0],r);
	for(long i=1; i<_cols; i++)
	  r+=_scols[i].multiply(a[i],b);
	return r;
#ifdef DEBUG
      } catch(domain_error e) { 
	cerr << e.what() << endl;
	throw domain_error(" called from refvector<TN> mat_full::operator*(const refvector<TN>& a) const");
      }
#endif
    }
  
  //! Multiply a full matrix by a full vector.
  template <class TN> template<class TN2>
    refvector<TN> mat_full<TN>::operator*(const refvector<TN2>& a) const
    {
      refvector<TN> r(_rows);
      return multiply(a,r);
    }
  
  //! Multiplies the \b left hand side from the \underline{ left } by the \b right hand side.
  template <class TN> mat_full<TN>& mat_full<TN>::operator*=(const mat_full<TN>& b)
    {
#ifdef DEBUG
      if( b._cols != _rows)
	{
	  throw domain_error( "Incompatible dimensions in mat_full mat_full::operator*=(const mat_full& b)");
	  exit(1);
	}
      try {
#endif
	long i;
	for(i=0; i<_cols;i++)
	  {
	    _scols[i]=b*_scols[i];
	  }
	_rows=b._rows;
	return *this;
#ifdef DEBUG
      } catch(domain_error e) { 
	cerr << e.what() << endl;
	throw domain_error(" called from mat_full<TN>& mat_full::operator*=(const mat_full<TN>& b)");
      }
#endif
    }

  //! Add two matrices.
  template <class TN> mat_full<TN>& mat_full<TN>::operator+=(const mat_full<TN>& b)
    {
#ifdef DEBUG
      if (b.cols() != _cols || b.rows() != _rows)
	throw domain_error("Incompatible dimensions in mat_full::operator+=(const mat_full<TN>& b)");
#endif
      for(long i=0; i<_cols;i++)
	(*this)[i]+=b[i];
      return *this;
    }


  //! Add two matrices.
  template <class TN> mat_full<TN>& mat_full<TN>::operator+=(const mat_sym_full<TN>& b)
    {
#ifdef DEBUG
      if (b.cols() != _cols || b.rows() != _rows)
	throw domain_error("Incompatible dimensions in mat_full::operator+=(const mat_sym_full<TN>& b)");
#endif
      for(long i=0; i<_cols;i++) {
	for(long j=0;j<=i;j++)
	  (*this)[i][j]+=b[i*(i+1)/2+j];
	for(long j=0;j<i;j++)
	  (*this)[j][i]+=b[i*(i+1)/2+j];
      }
      return *this;
    }

  //! Subtract two matrices.
  template <class TN> mat_full<TN>& mat_full<TN>::operator-=(const mat_full<TN>& b)
    {
#ifdef DEBUG
      if (b.cols() != _cols || b.rows() != _rows)
	throw domain_error("Incompatible dimensions in mat_full::operator-=(const mat_full<TN>& b)");
#endif
      for(long i=0; i<_cols;i++)
	(*this)[i]-=b[i];
      return *this;
    }

  //! Add two matrices.
  template <class TN> mat_full<TN>& mat_full<TN>::operator-=(const mat_sym_full<TN>& b)
    {
#ifdef DEBUG
      if (b.cols() != _cols || b.rows() != _rows)
	throw domain_error("Incompatible dimensions in mat_full::operator-=(const mat_sym_full<TN>& b)");
#endif
      for(long i=0; i<_cols;i++) {
	for(long j=0;j<=i;j++)
	  (*this)[i][j]-=b[i*(i+1)/2+j];
	for(long j=0;j<i;j++)
	  (*this)[j][i]-=b[i*(i+1)/2+j];
      }
      return *this;
    }

  //! Add two matrices.
  template <class TN> mat_full<TN> mat_full<TN>::operator+(const mat_full<TN>& b) const
    {
#ifdef DEBUG
      if (b.cols() != _cols || b.rows() != _rows)
	throw domain_error("Incompatible dimensions in mat_full::operator+=(const mat_full<TN>& b)");
#endif
      mat_full<TN> M(_cols,_rows);
      for(long i=0; i<_cols;i++)
	M[i]=(*this)[i]+b[i];
      return M;
    }
  //! Subtract two matrices.
  template <class TN> mat_full<TN> mat_full<TN>::operator-(const mat_full<TN>& b) const
    {
#ifdef DEBUG
      if (b.cols() != _cols || b.rows() != _rows)
	throw domain_error("Incompatible dimensions in mat_full::operator+=(const mat_full<TN>& b)");
#endif
      mat_full<TN> M(_cols,_rows);
      for(long i=0; i<_cols;i++)
	M[i]=(*this)[i]-b[i];
      return M;
    }
  
  // Normal member functions ****************************************************************************************************

  //! Copy matrix with full depth.
  template <class TN> matrix<TN>& mat_full<TN>::copy(const matrix<TN>& a)
  {
    return copy(a);
  }

  //! Copy matrix with full depth.
  template <class TN> mat_full<TN>& mat_full<TN>::copy(const mat_full<TN>& a)
    {
      _scols.copy(a._scols);
      _rows=a._rows;
      _cols=a._cols;
      return *this;
    }
  
  //! Copy matrix with depth i.
  template <class TN> matrix<TN>& mat_full<TN>::copy(const matrix<TN>& a, long i)
  {
    return copy(a,i);
  }
  //! Copy matrix with depth i.
  template <class TN> mat_full<TN>& mat_full<TN>::copy(const mat_full<TN>& a, const long i)
  {
    _scols.copy(a._scols,i);
    _rows=a._rows;
    _cols=a._cols;
    return *this;
  }

  //! Copy matrix with depth i.
  template <class TN> matrix<TN>& mat_full<TN>::copy(matrix<TN>& a, long i)
  {
    return copy(a,i);
  }
  //! Copy matrix with depth i.
  template <class TN> mat_full<TN>& mat_full<TN>::copy(mat_full<TN>& a, const long i)
  {
    _scols.copy(a._scols,i);
    _rows=a._rows;
    _cols=a._cols;
    return *this;
  }

  //! Adds a value.
  template <class TN> bool mat_full<TN>::add(const long x, const long y, const TN z)
    /*!
      @param z is added to (x,y).
      If the value is zero the associated vectors are expanded.
    */
    {

#ifdef DEBUG
      if(x<_cols && y<_rows)
	{
#endif
	  _scols[x][y]+=z;
	  return true;
#ifdef DEBUG
	}
      throw domain_error("Outside of dimension of mat_full::add " );
      exit(1);
#endif
    }
  
  //! Set a value.
  template <class TN> void mat_full<TN>::set(const long x, const long y, TN z)
    {
      (*this)[x][y]=z;
    }

  //! Simple display of matrix.
  template <class TN> bool mat_full<TN>::display(ostream& os) const
    {
      unsigned long val=0;
      unsigned long v=0;
      while(val< _scols->size())
	{
	  v=0;
	  cout << "( ";
	  while(v<_scols[val]->size())
	    os << _scols[val][v++] << " ";
	  os << " )\n";
	  val++;
	}
      return true;
    }

  //! Simple display of matrix.
  template <class TN> bool mat_full<TN>::display() const
    {
      return display(cout);
    }
  
  
  //! Extend the matrix by n columns and m rows.
  template <class T> mat_full<T>& mat_full<T>::grow(long n, long m)
  {
    for(long i=0;i<_cols;i++)
      _scols[i].resize(_rows+m);
    const refvector<T> N(_rows+m);
    _scols.resize(_cols+n,N);
    _cols+=n;
    _rows+=m;
    return *this;
  }

  //! Returns a mat_sym_full object of the upper triangle of a full matrix.
  template <class TN> mat_sym_full<TN>& mat_full<TN>::sym_upper_triangle(mat_sym_full<TN>& r) const
    {
#ifdef DEBUG
      if(r.cols() != _cols && r.cols()!=_cols)
	throw domain_error("mat_full::sym_upper_triangle(mat_sym_full<TN>& r) const: r has wrong dimension");
      try {
#endif
	long i,j;
	for(i=0; i<r.cols(); i++)
	  for(j=0; j<=i; j++)
	    r[i*(i+1)/2+j]=_scols[i][j];
	
	return r;
#ifdef DEBUG
      } catch(domain_error e) { 
	cerr << e.what() << endl;
	throw domain_error(" called from  mat_sym_full<TN> mat_full::sym_upper_triangle() const");
      }
#endif
    }


  //! Returns a mat_sym_full object of the upper triangle of a full matrix.
  template <class TN> mat_sym_full<TN> mat_full<TN>::sym_upper_triangle() const
  {
    long m;
    if(_rows> _cols) m=_cols;
    else m=_rows;
    mat_sym_full<TN> r(m);
    return sym_upper_triangle(r);
  }

  //! Returns a mat_sym_full object of the lower triangle of a full matrix.
  template <class TN> mat_sym_full<TN>& mat_full<TN>::sym_lower_triangle(mat_sym_full<TN>& r) const
  {
#ifdef DEBUG
      try {
	if(r.cols() != _cols && r.cols()!=_cols)
	  throw domain_error("mat_full::sym_lower_triangle(mat_sym_full<TN>& r) const: r has wrong dimension");
#endif
	long i,j;
	for(i=0; i<r.cols(); i++)
	  for(j=0; j<=i; j++)
	    r[i*(i+1)/2+j]=_scols[j][i];
	
	return r;
#ifdef DEBUG
      } catch(domain_error e) { 
	cerr << e.what() << endl;
	throw domain_error(" called from  mat_sym_full<TN> mat_full::sym_lower_triangle() const");
      }
#endif
  }


  //! Returns a mat_sym_full object of the lower triangle of a full matrix.
  template <class TN> mat_sym_full<TN> mat_full<TN>::sym_lower_triangle() const
  {
    long m;
    if(_rows> _cols) m=_cols;
    else m=_rows;
    mat_sym_full<TN> r(m);
    return sym_lower_triangle(r);
  }
  
  //! Returns a mat_asym_full object of the upper triangle of a full matrix.
  template <class TN> mat_asym_full<TN> mat_full<TN>::asym_upper_triangle() const
    {
      try {
	long i,j;
	long m;
	if(_rows> _cols) m=_cols;
	else m=_rows;
	mat_asym_full<TN> r(m);
	for(i=0; i<m; i++)
	  for(j=0; j<i; j++)
	    r[i*(i+1)/2+j]=_scols[i][j];
	
	return r;
      } catch(domain_error e) { 
	cerr << e.what() << endl;
	throw domain_error(" called from  mat_sym_full<TN> mat_full::sym_upper_triangle() const");
      }
    }

  //! Transposes the full matrix.
  template <class TN> mat_full<TN>& mat_full<TN>::transpose(mat_full<TN>& r) const
  {
#ifdef DEBUG
    if(r.cols()!=_rows || r.rows()!=_cols)
      throw domain_error("");
    try {
#endif
      long i,j;
      for(i=0;i<_cols;i++)
	for(j=0;j<_rows; j++) r[j][i]=_scols[i][j];
      return r;
#ifdef DEBUG
    } catch (domain_error e) {
      cerr << e.what() << endl;
      throw domain_error("called from mat_full::transpose(..)");
    }
#endif
  }

  //! Transposes the full matrix.
  template <class TN> mat_full<TN> mat_full<TN>::transpose() const
  {
    mat_full<TN> r(_rows,_cols);
    return transpose(r);
  }

  //! Release the associated memory.
  template <class TN> void mat_full<TN>::clear()
    {
      refvector<refvector<TN> > a;
      _scols=a;
      _rows=_cols=0;
      _size=0;
    }

  //! Prune the matrix of negligible values.
  template <class TN> mat_full<TN>& mat_full<TN>::prune(TN a)
    {
      long i,j;
      for(i=0;i<_cols;i++)
	for(j=0;j<_rows;j++)
	  if(_scols[i][j]<a && _scols[i][j]>-a) _scols[i][j]=(TN) 0;
      return *this;
    }

  //! Set the matrix to zero.
  template <class TN> mat_full<TN>& mat_full<TN>::zero()
  {
    long i,j;
    linear_algebra::zero(_scols);
    return (*this);
  }

  //! Set the matrix to the identity.
  template <class TN> mat_full<TN>& mat_full<TN>::one()
  {
    long i,j;
    zero();
    for(i=0;i<min(_cols,_rows);i++)
      linear_algebra::one(_scols[i][i]);
    return (*this);
  }
  
  //! Do half a householder and then do gaussian elimination for a left inverse.\todo Throw exceptions for ill-conditioned matrices
  template <class TN> mat_full<double>& mat_full<TN>::gauss_elim(mat_full<TN>& T)
  {
#ifdef DEBUG
    if(_rows != _cols)
      throw domain_error("gauss_elim only works on square matrices");
    if(T.cols()!= _cols || T.rows() != _rows)
      throw domain_error("mat_full::gauss_elim has dimensional quibbles");
    try {
#endif
      long i,j,k;
      T.one();
      TN dummy;
      TN absval;
      TN absval2;

      // half householder
      refvector<TN> v(_rows);
      for(i=0;i<_cols-1;i++)
	{
	  for(j=i;j<v.dim();j++) v[j]=(*this)[i][j];
	  absval=sqrt(v*v);
	  if(v[i]<0.0)
	    {
	      v[i]-=absval;
	      (*this)[i][i]=absval;
	    }
	  else {
	    v[i]+=absval;
	    (*this)[i][i]=-absval;
	  }
	  if(absval>1e-16*fabs((*this)[i+1][i])) {
	    v*=sqrt((TN) 2/(v*v));
	    for(k=i+1;k<_rows;k++)
	      (*this)[i][k]=(TN) 0;
	    for(k=i+1;k<_cols;k++)
	      {
		dummy=(v*(*this)[k]);
		for(j=i;j<_rows;j++)
		  (*this)[k][j]-=v[j]*dummy;
	      }
	    for(k=0;k<T.cols();k++)
	      {
		dummy=T[k]*v;
		for(j=i;j<T.rows();j++)
		  T[k][j]-=v[j]*dummy;
	      }
	  }
	  v[i]=(TN) 0;
	}
      cout << " After householder T:\n";
      T.display();
      
      cout << " After householder this:\n";
      prune(1e-15).display();

      // Gaussian down the line; only one entry due to householder
      // And back up
      for(i=T.cols()-1;i>=0;i--)
	{
	  for(j=0;j<T.cols();j++)
	    {
	      T[j][i]/=(*this)[i][i];
	      for(k=0;k<i;k++)
		T[j][k]-=(*this)[i][k]*T[j][i];
	    }
	}

      return T;
#ifdef DEBUG
    } catch (domain_error e) {
      cerr << e.what() << endl;
      throw domain_error("called from mat_full::gauss_elim");
    }
#endif
  }

} // end linear_algebra::
#include <refcount.h>
#include <matrix.h>
#include <sparse_vector.h>
#include <mat_sym_full.h>
#include <mat_asym_full.h>
#include <mat_sym_sparse.h>

#endif
/*! @}*/
