
// -------------------------------------------------------------------
#include <iostream>

// -------------------------------------------------------------------
int main()
{
  typedef Matrix<float, 1, 4> Matrix1f;
  typedef Matrix<float, 4, 4> Matrix4f;

  Matrix1f m1;
  m1(0,0) = 27;
  std::cout << m1;

  Matrix4f m4;
  std::cout << inverse(m4);

  return 0;
}
