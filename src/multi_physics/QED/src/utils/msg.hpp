#ifndef PICSAR_MULTIPHYSICS_MSG
#define PICSAR_MULTIPHYSICS_MSG

//Should be included by all the src files of the library
#include "../qed_commons.h"

#include<iostream>

namespace picsar{
namespace multi_physics{
namespace utils{

    //A message function which writes to a stream (provided that the
    //pointer is not NULL)
    void msg(const std::string& msg, std::ostream* stream)
    {
        if(stream != nullptr)
            *stream << msg;
    }

    //An error message function which writes to stderr
    void err(const std::string& msg)
    {
         std::cerr << msg;
    }
}
}
}

#endif //PICSAR_MULTIPHYSICS_MSG
