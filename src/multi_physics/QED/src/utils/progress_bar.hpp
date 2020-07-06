#ifndef PICSAR_MULTIPHYSICS_PROGRESS_BAR
#define PICSAR_MULTIPHYSICS_PROGRESS_BAR

//Should be included by all the src files of the library
#include "../qed_commons.h"

#include <string>
#include <iostream>

namespace picsar{
namespace multi_physics{
namespace utils{

template<typename OutStream = std::ostream>
void draw_progress(
    const int i, const int how_many,
    const std::string text,
    const int up_freq = 1,
    bool last=false,
    OutStream& out = &std::cout)
{
    if (i % up_freq != 0 && i != how_many)
        return;

    const auto bar_width = 50;
    const auto progress = (i*1.0/how_many);
    const auto pos = static_cast<int>(bar_width*progress);
    out << text << " [";
    for (int j = 0; j < bar_width; ++j) {
        if (j < pos) out << "=";
        else if (j == pos) out << ">";
        else out << " ";
    }
    out << "] " << static_cast<int>(progress * 100.0);
    if(last) out <<" %\n";
    else out <<" %\r";
    out.flush();
}

}
}
}

#endif //PICSAR_MULTIPHYSICS_PROGRESS_BAR