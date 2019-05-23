#ifndef __PICSAR_MULTIPHYSICS_MSG__
#define __PICSAR_MULTIPHYSICS_MSG__

//This .hpp file contains small helper functions to write a message to an
//ostream (which might be std::cout) or to std::err for error messages

#include<string>
#include<iostream>

//Should be included by all the src files of the library
#include "qed_commons.h"

//############################################### Declaration

namespace picsar{
    namespace multi_physics{
        //A message function which writes to a stream (provided that the
        //pointer is not NULL)
        PXRMP_FORCE_INLINE
        void msg(const std::string& msg, std::ostream* stream);

        //An error message function which writes to stderr
        PXRMP_FORCE_INLINE
        void err(const std::string& msg);
    }
}

//############################################### Implementation

PXRMP_FORCE_INLINE
void picsar::multi_physics::msg(const std::string& msg,
    std::ostream* stream)
    {
        if(stream != nullptr)
            *stream << msg;
    }

PXRMP_FORCE_INLINE
void picsar::multi_physics::err(const std::string& msg)
    {
        std::cerr << msg;
    }

#endif //__PICSAR_MULTIPHYSICS_MSG__
