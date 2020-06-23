// Wigner little-d functions for half-integer spin particles
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "misc_math.hpp"

// ---------------------------------------------------------------------------
double jpacPhoto::wigner_error(int j, int lam1, int lam2, bool half)
{
  // if (half == true)
  // {
  //   std::cout << "\n";
  //   std::cout << "wigner_d: Argument combination with j = " << j << "/2";
  //   std::cout << ", lam1 = " << lam1 << "/2";
  //   std::cout << ", lam2 = " << lam2 << "/2";
  //   std::cout << " does not exist. Quitting... \n";
  // }
  // else
  // {
  //   std::cout << "\n";
  //   std::cout << "wigner_d: Argument combination with j = " << j;
  //   std::cout << ", lam1 = " << lam1;
  //   std::cout << ", lam2 = " << lam2;
  //   std::cout << " does not exist. Quitting... \n";
  // }
  //
  // exit(0);

  return 0.;
}
// ---------------------------------------------------------------------------
// USING WIKIPEDIA SIGN CONVENTION
// theta is in radians
// lam1 = 2 * lambda and lam2 = 2 * lambda^prime are integers
double jpacPhoto::wigner_d_half(int j, int lam1, int lam2, double theta)
{
  double phase = 1.;

  // If first lam argument is smaller, switch them
  if (abs(lam1) < abs(lam2))
  {
    int temp = lam1;
    lam1 = lam2;
    lam2 = temp;

    phase *= pow(-1., double(lam1 - lam2) / 2.);
  };

  // If first lam is negative, smitch them
  if (lam1 < 0)
  {
    lam1 *= -1;
    lam2 *= -1;

    phase *= pow(-1., double(lam1 - lam2) / 2.);
  }

  std::complex<double> result = 0.;
  switch (j)
  {
    // -------------------------------------------------------------------------
    // j = 1/2
    // -------------------------------------------------------------------------
    case 1:
    {
      if (lam1 == 1)
      {
        if (lam2 == 1)
        {
          result = cos(theta / 2.);
          break;
        }
        else
        {
          result = -sin(theta / 2.);
          break;
        }
      }
      else
      {
        wigner_error(j, lam1, lam2, true);
      }
      break;
    }

    // -------------------------------------------------------------------------
    // j = 3/2
    // -------------------------------------------------------------------------
    case 3:
    {
      // Lambda = 3/2
      if (lam1 == 3)
      {
        switch (lam2)
        {
          case 3:
          {
             result = cos(theta / 2.) / 2.;
             result *= (1. + cos(theta));
             break;
          }
          case 1:
          {
            result = - sqrt(3.) / 2.;
            result *= sin(theta / 2.);
            result *= 1. + cos(theta);
            break;
          }
          case -1:
          {
            result = sqrt(3.) / 2.;
            result *= cos(theta / 2.);
            result *= 1. - cos(theta);
           break;
          }
          case -3:
          {
            result = - sin(theta / 2.) / 2.;
            result *= 1. - cos(theta);
            break;
          }
          default: wigner_error(j, lam1, lam2, true);
        }
      }

      // Lambda = 1/2
      else if (lam1 == 1)
      {
        switch (lam2)
        {
          case 1:
          {
            result = 1. / 2.;
            result *= 3. * cos(theta) - 1.;
            result *= cos(theta / 2.);
            break;
          }
          case -1:
          {
            result = -1. / 2.;
            result *= 3. * cos(theta) + 1.;
            result *= sin(theta / 2.);
            break;
          }
          default: wigner_error(j, lam1, lam2, true);
        }
      }

      // Error
      else
      {
        wigner_error(j, lam1, lam2, true);
      }
      break;
    }

    // -------------------------------------------------------------------------
    // j = 5/2
    // -------------------------------------------------------------------------
    case 5:
    {
      switch (lam1)
      {
        // lambda = 5/2 not yet implemented
        case 5:
        {
          wigner_error(j, lam1, lam2, true);
        }
        // lam1 == 3
        case 3:
        {
          switch (lam2)
          {
            case 3:
            {
              result = -1. / 4.;
              result *= cos(theta / 2.);
              result *= (1. + cos(theta)) * (3. - 5. * cos(theta));
              break;
            }
            case 1:
            {
              result = sqrt(2.) / 4.;
              result *= sin(theta / 2.);
              result *= (1. + cos(theta) * (1. - 5. * cos(theta));
              break;
            }
            case -1:
            {
              result = - sqrt(2.) / 4.;
              result *= cos(theta / 2.);
              result *= (1. - cos(theta) * (1. + 5. * cos(theta));
              break;
            }
            case -3:
            {
              result = 1. / 4.;
              result *= sin(theta / 2.);
              result *= (1. - cos(theta)) * (3. + 5. * cos(theta));
              break;
            }
            default: wigner_error(j, lam1, lam2, true);
          }
          break;
        }
        // lam1 == 1
        case 1:
        {
          switch (lam2)
          {
            case 1:
            {
              result = -1. / 2.;
              result *= cos(theta / 2.);
              result *= (1. + 2. * cos(theta) - 5. * cos(theta)*cos(theta));
              break;
            }
            case -1:
            {
              result = 1. / 2.;
              result *= sin(theta / 2.);
              result *= (1. - 2. * cos(theta) - 5. * cos(theta)*cos(theta));
              break;
            }
            default: wigner_error(j, lam1, lam2, true);
          }
          break;
        }
        default: wigner_error(j, lam1, lam2, true);
      }
      break;
    }
    // Error
    default: wigner_error(j, lam1, lam2, true);
  }

  return phase * result;
};

// ---------------------------------------------------------------------------
std::complex<double> jpacPhoto::wigner_d_int(int j, int lam1, int lam2, std::complex<double> z)
{

  double phase = 1.;
  // If first lam argument is smaller, switch them
  if (abs(lam1) < abs(lam2))
  {
    int temp = lam1;
    lam1 = lam2;
    lam2 = temp;

    phase *= pow(-1., double(lam1 - lam2));
  };

  // If first lam is negative, smitch them
  if (lam1 < 0)
  {
    lam1 *= -1;
    lam2 *= -1;

    phase *= pow(-1., double(lam1 - lam2));
  }

  std::complex<double> result = 0.;
  switch (j)
  {
    // spin - 1
    case 1:
    {
      if (lam1 == 1)
      {
        switch (lam2)
        {
          case 1:
          {
            result = (1. + cos(theta)) / 2.;
            break;
          }
          case 0:
          {
            result = - sin(theta) / sqrt(2.);
            break;
          }
          case -1:
          {
            result = (1. - cos(theta)) / 2.;
            break;
          }
          default: wigner_error(j, lam1, lam2, false);
        }
      }
      else if (lam1 == 0)
      {
        result = cos(theta);
      }
      else
      {
        wigner_error(j, lam1, lam2, false);
      }
      break;
    }
  // Error
  default: wigner_error(j, lam1, lam2, false);
}

return phase * result;
};
