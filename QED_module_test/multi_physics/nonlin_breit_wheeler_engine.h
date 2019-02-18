#ifndef __PMP_BREITWHEELER__
#define __PMP_BREITWHEELER__

#include <cstdint>
#include <random>

namespace picsar{
  namespace multi_physics{

    class nonlin_breit_wheeler_engine{
    public:
      nonlin_breit_wheeler_engine(int64_t seed);

    private:
        std::mt19937_64 rng;


    };
    
  }
}

#endif
