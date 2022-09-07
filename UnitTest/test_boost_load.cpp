
#include <iostream>

//#include <boost>
#include <boost/optional.hpp>
#include <boost/optional/optional_io.hpp>

#define BOOST_TEST_MODULE OnlyBoostLoad
#include <boost/test/unit_test.hpp>

using std::cout;
using std::endl;

BOOST_AUTO_TEST_CASE( optional )
{
  boost::optional<int> a = boost::none;
  cout << a << endl;
}


