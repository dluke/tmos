
#include <iostream>


#include <unordered_map>

#define BOOST_TEST_MODULE OnlyStdLib
#include <boost/test/unit_test.hpp>

using std::cout;
using std::endl;

// declaration

struct mapholder {
  static std::unordered_map<int,int> base;
};


// definition

std::unordered_map<int,int> mapholder::base = {{1,2},{3,4}};

BOOST_AUTO_TEST_CASE( point )
{
  cout << mapholder::base[1] << endl;
}



