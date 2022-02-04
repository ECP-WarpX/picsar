#ifndef PICSAR_MULTIPHYSICS_PROGRESS_BAR
#define PICSAR_MULTIPHYSICS_PROGRESS_BAR

//Should be included by all the src files of the library
#include "picsar_qed/qed_commons.h"

#include <string>
#include <iostream>

namespace picsar{
namespace multi_physics{
namespace utils{

    /**
    * A simple progress bar
    *
    * @param[in] i current progress index
    * @param[in] how_many maximum value of the progress index
    * @param[in] text an optional text to append to the progress bar
    * @param[in] up_freq frequency at which the progress bar is updated
    * @param[in] last if true the last character is a new line instead of a carriage return
    * @param[in] out the std::ostream where the progress bar is drawn
    */
    inline
    void draw_progress(
        const int i, const int how_many,
        const std::string text = "",
        const int up_freq = 1,
        bool last=false,
        std::ostream& out = std::cout)
        {
            if (i % up_freq != 0 && i != how_many)
            return;

            const auto bar_length = 50;
            const auto progress = (i*1.0/how_many);
            const auto pos = static_cast<int>(bar_length*progress);
            out << " [";
            for (int j = 0; j < bar_length; ++j) {
                if (j < pos) out << "=";
                else if (j == pos) out << ">";
                else out << " ";
            }
            const int progress_percentage = static_cast<int>(progress * 100.0);
            out << "] " << progress_percentage << "%  " << text ;
            if(last) out <<"\n";
            else out <<"\r";
            out.flush();
        }

}
}
}

#endif //PICSAR_MULTIPHYSICS_PROGRESS_BAR
