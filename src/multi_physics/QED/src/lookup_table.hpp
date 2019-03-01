#ifndef __PMP_LOOKUP__
#define __PMP_LOOKUP__

#include <stdexcept>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <array>

#include <iostream>

//This source file contains a templated lookup table (work in progress)
//It allows to do this:
//vector<double> c1{1,2,3,4,5,6,7,8,9,10};
//vector<double> c2{-1.0,0,1.0};
//vector<double> c3{-100,0.0,100};
//picsar::multi_physics::lookup_table<3, double, double> table{c1, c2, c3};

namespace picsar{
    namespace multi_physics{

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        class lookup_table{
        public:
            lookup_table();

            template <typename... COORDS>
            lookup_table(COORDS... _coords);

            template <typename... COORDS>
            lookup_table(COORDS... _coords, std::vector<DATA_TYPE>&& _raw_data);

            bool fill_at(std::array<COORD_TYPE, DIM> where, DATA_TYPE what);
            bool get_at(std::array<COORD_TYPE, DIM> where, DATA_TYPE& what);

            size_t get_dims();
            size_t get_data_size();
            std::vector<COORD_TYPE> get_coords(size_t dim);

            void print_on_disk(std::string file_name);

        private:
            std::array<std::vector<COORD_TYPE>, DIM> coords;
            size_t total_size = 0;
            std::vector<DATA_TYPE> raw_data;

            template <typename COORDS>
            void init_lookup_table(size_t counter, COORDS _coords);

            template <typename COORDS, typename... REST>
            void init_lookup_table(size_t counter, COORDS _coords, REST... _other_coords);

            void do_coords_checks();
            bool where_to_pos(const std::array<COORD_TYPE, DIM>& where, std::array<size_t, DIM>& pos);
            void compute_total_size();

            PXR_FORCE_INLINE
            size_t compute_index(std::array<size_t, DIM> pos);
        };

        //##############################################################Implementation

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        lookup_table<DIM, COORD_TYPE, DATA_TYPE>::lookup_table(){}

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        template <typename... COORDS>
        lookup_table<DIM, COORD_TYPE, DATA_TYPE>::lookup_table(COORDS... _coords){
            static_assert(sizeof...(_coords) == DIM, "Wrong number of arguments.");
            init_lookup_table(0, _coords...);
            do_coords_checks();
            compute_total_size();

            total_size = 1;
            for(auto const& clist: coords)
                total_size *= clist.size();

            raw_data.resize(total_size);
        }

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        template <typename... COORDS>
        lookup_table<DIM, COORD_TYPE, DATA_TYPE>::lookup_table(COORDS... _coords, std::vector<DATA_TYPE>&& _raw_data){
            static_assert(sizeof...(_coords) == DIM, "Wrong number of arguments.");
            init_lookup_table(0, _coords...);
            do_coords_checks();
            compute_total_size();

            if(total_size == _raw_data.size())
                raw_data = _raw_data;
            else
                throw std::logic_error(std::string("Raw data has the wrong size!"));

            raw_data.resize(total_size);
        }

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        bool lookup_table<DIM, COORD_TYPE, DATA_TYPE>::fill_at(std::array<COORD_TYPE, DIM> where, DATA_TYPE what){
            std::array<size_t, DIM> pos;
            bool success = where_to_pos(where, pos);
            if(success)
                raw_data[compute_index(pos)] = what;

            return success;
        }

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        bool lookup_table<DIM, COORD_TYPE, DATA_TYPE>::get_at(std::array<COORD_TYPE, DIM> where, DATA_TYPE& what){
            std::array<size_t, DIM> pos;
            bool success = where_to_pos(where, pos);
            if(success)
                what = raw_data[compute_index(pos)];

            return success;
        }

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        template<typename COORDS>
        void lookup_table<DIM, COORD_TYPE, DATA_TYPE>::init_lookup_table(size_t counter, COORDS _coords){
            coords[counter] = _coords;
        }
        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        template<typename COORDS, typename... REST>
        void lookup_table<DIM, COORD_TYPE, DATA_TYPE>::init_lookup_table(size_t counter, COORDS _coords, REST... _other_coords){
            coords[counter] = (_coords);
            init_lookup_table(++counter, _other_coords...);
        }

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        size_t lookup_table<DIM, COORD_TYPE, DATA_TYPE>::get_dims(){
            return coords.size();
        }

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        size_t lookup_table<DIM, COORD_TYPE, DATA_TYPE>::get_data_size(){
            return raw_data.size();
        }

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        std::vector<COORD_TYPE> lookup_table<DIM, COORD_TYPE, DATA_TYPE>::get_coords(size_t dim){
            return coords.at(dim);
        }

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        void lookup_table<DIM, COORD_TYPE, DATA_TYPE>::do_coords_checks(){
            for(auto const& clist: coords){
                if (!std::is_sorted(clist.begin(), clist.end()))
                    throw std::logic_error(std::string("Coordinates must be sorted!"));
            }

            for (auto const& clist: coords){
                for (auto it = clist.begin() + 1; it != clist.end(); ++it){
                    if (*it == *(it - 1))
                        throw std::logic_error(std::string("Coordinates contains duplicates!"));
                }
            }
        }

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        void lookup_table<DIM, COORD_TYPE, DATA_TYPE>::compute_total_size(){
            total_size = 1;
            for(auto const& clist: coords)
                total_size *= clist.size();
        }


        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        size_t lookup_table<DIM, COORD_TYPE, DATA_TYPE>::compute_index(std::array<size_t, DIM> pos){
            size_t index = 0;
            size_t mul = 1;
            for(auto it = pos.rbegin(); it != pos.rend();  it++){
                index += mul*(*it);
                mul *= coords[ std::distance(pos.begin(),it.base())-1].size();
            }
            return index;
        }

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        bool lookup_table<DIM, COORD_TYPE, DATA_TYPE>::where_to_pos(const std::array<COORD_TYPE, DIM>& where, std::array<size_t, DIM>& pos){
            for(size_t i = 0; i < where.size(); i++){
                auto lower = std::lower_bound(coords[i].begin(), coords[i].end(), where[i]);
                if(lower != coords[i].end() && *lower == where[i]){
                    pos[i] = std::distance(coords[i].begin(),lower);
                }
                else{
                    return false;
                }
            }
            return true;
        }

        template<size_t DIM, typename COORD_TYPE, typename DATA_TYPE>
        void lookup_table<DIM, COORD_TYPE, DATA_TYPE>::print_on_disk(std::string file_name){
            std::ofstream of;
            of.open(file_name.c_str());

            std::array<size_t, DIM> base_vec;

            size_t m = 1;
            for(auto it = base_vec.rbegin(); it != base_vec.rend();  it++){
                *it = m;
                m *= coords[ std::distance(base_vec.begin(),it.base())-1].size();
            }

            for(auto raw = raw_data.begin(); raw != raw_data.end(); raw++){
                size_t indx = std::distance(raw_data.begin(), raw);
                for(auto c_list = coords.begin(); c_list != coords.end(); c_list++){
                    auto c_index = std::distance(coords.begin(), c_list);
                    auto cc_index = indx/base_vec[c_index];
                    of << (*c_list)[cc_index] << " " ;
                    indx -= cc_index*base_vec[c_index];
                }
                of << *raw << std::endl;
            }

            of.close();
        }

    }
}



#endif
