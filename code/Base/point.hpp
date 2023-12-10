
#ifndef POINT_HPP_
#define POINT_HPP_


#include <iostream>
#include <array>            // std::array

// lost the link from the webpage I borrowed this Point class
// It's distinct from Vector3d because it only uses integer arithmetic and is templated

#define STATIC_ASSERT( e ) static_assert( e, "!(" #e ")" )

template< typename T, int nDimensions = 2 >
class Point
{
private:
    std::array< T, nDimensions > elements_;

public:
    typedef T ValueType;

    T operator[]( int const i )
    {
        return elements_[i];
    }

    T const operator[]( int const i ) const
    {
        return elements_[i];
    }

    bool operator==( Point const other ) const
    {
      int re = 1;
      for( int i = 0; i < nDimensions; ++i )
      {
        re = (re && (elements_[i] == other.elements_[i]));
      }
      return re;
    }

    bool operator!=( Point const other ) const
    {
      return !(*this == other);
    }
    
    void operator+=( Point const other )
    {
        for( int i = 0; i < nDimensions; ++i )
        {
            elements_[i] += other.elements_[i];
        }
    }

    void operator-=( Point const other )
    {
        for( int i = 0; i < nDimensions; ++i )
        {
            elements_[i] -= other.elements_[i];
        }
    }

    // elementwise multiplication
    void operator*=( T const c)
    {
        for( int i = 0; i < nDimensions; ++i )
        {
            elements_[i] *= c;
        }
    }

    // elementwise division
    void operator /=( T const c)
    {
        for( int i = 0; i < nDimensions; ++i )
        {
            elements_[i] /= c;
        }
    }
    
    friend Point operator+( Point const a, Point const b )
    {
        Point ret( a );
        ret += b;
        return ret;
    }

    friend Point operator-( Point const a, Point const b )
    {
        Point ret( a );
        ret -= b;
        return ret;
    }

    Point operator*( T const c) {
      Point ret(*this);
      ret *= c;
      return ret;
    }

    friend Point operator*( T const c, Point const p)
    {
      Point ret(p);
      ret *= c;
      return ret;
    }

    Point operator/ ( T const c) const
    {
      Point ret(*this);
      ret /= c;
      return ret;
    }

    //! Negate a vector
    Point operator-() const
    {
      std::array<T, nDimensions> neg;
      for(int i = 0; i < nDimensions; ++i) {
        neg[i] = -elements_[i];
      }
      return Point(neg);
    }


    friend std::ostream& operator<<(std::ostream& os, const Point<T, nDimensions>& point)
    {
      os << "(";
      int i = 0;
      for (; i < nDimensions-1; ++i)
      {
        os << point.elements_[i] << ",";
      }
      os << point.elements_[i] << ")";
      return os;
    }

    Point(): elements_() {}

    Point( std::array< T, nDimensions> els ) {
      elements_ = els;
    }

    Point( int x, int y )
    {
        STATIC_ASSERT( nDimensions == 2 );
        elements_[0] = x;
        elements_[1] = y;
    }


    Point( int x, int y, int z )
    {
        STATIC_ASSERT( nDimensions == 3 );
        elements_[0] = x;
        elements_[1] = y;
        elements_[2] = z;
    }

};

typedef Point< int, 2 > Point2d;
typedef Point< int, 3 > Point3d;



#endif 
