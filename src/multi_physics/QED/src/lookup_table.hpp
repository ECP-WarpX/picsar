#ifndef __PMP_LOOKUP__
#define __PMP_LOOKUP__

#include <vector>
#include <array>

namespace picsar{
    namespace multi_physics{

     template< size_t dims, class TC, class TD>
      class lookup_table{
      public:
          lookup_table();

      private:
          std::array<std::vector<TC>, dims> coords;
          std::vector<TD> raw_data;
      };

    }

    void small_test(){
        lookup_table<3,double,double> lookup;
    }


}



#endif __PMP_LOOKUP__
