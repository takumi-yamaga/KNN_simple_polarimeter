// MyConstants.h

namespace mass{
  // mass in MeV/c2
  
  const double pi_plus = 139.57018;
  const double pi_zero = 134.9766;
  
  const double proton  = 938.272081;
  const double neutron = 939.565413;

  const double lambda      = 1115.683;
  const double sigma_plus  = 1189.37;
  const double sigma_zero  = 1192.642;
  const double sigma_minus = 1197.449;

}

void PrintVector(const TVector3& vec){
    std::cout << vec.X() 
      << ", " << vec.Y()
      << ", " << vec.Z()
      << std::endl;
}
