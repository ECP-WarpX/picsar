#ifndef __PMP_LOOKUP__
#define __PMP_LOOKUP__

#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

#include <iostream>

//This source file contains a templated lookup table (work in progress)
//It allows to do this:
//vector<double> c1{1,2,3,4,5,6,7,8,9,10};
//vector<double> c2{-1.0,0,1.0};
//vector<double> c3{-100,0.0,100};
//picsar::multi_physics::lookup_table<3, double, double> table{c1, c2, c3};

namespace picsar{
    namespace multi_physics{

        template<typename COORD_TYPE, typename DATA_TYPE>
        class lookup_table{
        public:
            template <typename... COORDS>
            lookup_table(COORDS... _coords);

            size_t get_dims();

            size_t get_data_size();

        private:
            std::vector<std::vector<COORD_TYPE>> coords;
            size_t total_size = 0;
            std::vector<DATA_TYPE> raw_data;

            template <typename COORDS>
            void init_lookup_table(COORDS _coords);

            template <typename COORDS, typename... REST>
            void init_lookup_table(COORDS _coords, REST... _other_coords);
        };

        //##############################################################Implementation

        template<typename COORD_TYPE, typename DATA_TYPE>
        template <typename... COORDS>
        lookup_table<COORD_TYPE, DATA_TYPE>::lookup_table(COORDS... _coords){
            init_lookup_table(_coords...);

            for(auto& clist: coords){
                if (!std::is_sorted(clist.begin(), clist.end())){
                    throw std::logic_error(std::string("Coordinates must be sorted!"));
                }
            }

            total_size = 1;
            for(auto& clist: coords)
                total_size *= clist.size();

            raw_data.resize(total_size);
        }


        template<typename COORD_TYPE, typename DATA_TYPE>
        template<typename COORDS>
        void lookup_table<COORD_TYPE, DATA_TYPE>::init_lookup_table(COORDS _coords){
            coords.push_back(_coords);
        }
        template<typename COORD_TYPE, typename DATA_TYPE>
        template<typename COORDS, typename... REST>
        void lookup_table<COORD_TYPE, DATA_TYPE>::init_lookup_table(COORDS _coords, REST... _other_coords){
            coords.push_back(_coords);
            init_lookup_table(_other_coords...);
        }

        template<typename COORD_TYPE, typename DATA_TYPE>
        size_t lookup_table<COORD_TYPE, DATA_TYPE>::get_dims(){
            return coords.size();
        }

        template<typename COORD_TYPE, typename DATA_TYPE>
        size_t lookup_table<COORD_TYPE, DATA_TYPE>::get_data_size(){
            return raw_data.size();
        }


    }
}



#endif
