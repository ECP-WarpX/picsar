#ifndef __PMP_LOOKUP__
#define __PMP_LOOKUP__

#include <stdexcept>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <array>
#include <utility>
#include <iostream>

//This source file contains a templated lookup table (work in progress)
//It allows to do this:
//vector<double> c1{1,2,3,4,5,6,7,8,9,10};
//vector<double> c2{-1.0,0,1.0};
//vector<double> c3{-100,0.0,100};
//picsar::multi_physics::lookup_table<3, double, double> table{c1, c2, c3};

namespace picsar{
    namespace multi_physics{

        template<size_t DIM, typename DATA_TYPE>
        class lookup_table{
        public:
            lookup_table();

            template <typename... COORDS>
            lookup_table(COORDS... _coords);

            template <typename... COORDS>
            lookup_table(COORDS... _coords, std::vector<DATA_TYPE>&& _raw_data);

            bool fill_at(const std::array<DATA_TYPE, DIM>& where, DATA_TYPE what);
            bool get_at(const std::array<DATA_TYPE, DIM>& where, DATA_TYPE& what);

            bool interp_at(const std::array<DATA_TYPE, DIM>& where, DATA_TYPE& what);
            bool interp_at_one_dim(size_t which_dim, DATA_TYPE where, const std::array<size_t, DIM>& pos, DATA_TYPE& what);

            size_t get_dims();
            size_t get_data_size();
            std::vector<DATA_TYPE> get_coords(size_t dim);

            void print_on_disk(std::string file_name, bool reloadable);
            void read_from_disk(std::string file_name);

        private:
            std::array<std::vector<DATA_TYPE>, DIM> coords;
            size_t total_size = 0;
            std::vector<DATA_TYPE> raw_data;

            static const size_t SIZE_2N = (1 << DIM);

            template <typename COORDS>
            void init_lookup_table(size_t counter, COORDS _coords);

            template <typename COORDS, typename... REST>
            void init_lookup_table(size_t counter, COORDS _coords, REST... _other_coords);

            void do_coords_checks();
            bool where_to_pos(const std::array<DATA_TYPE, DIM>& where, std::array<size_t, DIM>& pos);
            void compute_total_size();

            PXR_FORCE_INLINE
            size_t compute_index(std::array<size_t, DIM> pos);
        };

        //##############################################################Implementation

        template<size_t DIM, typename DATA_TYPE>
        lookup_table<DIM, DATA_TYPE>::lookup_table(){}

        template<size_t DIM, typename DATA_TYPE>
        template <typename... COORDS>
        lookup_table<DIM, DATA_TYPE>::lookup_table(COORDS... _coords){
            static_assert(sizeof...(_coords) == DIM, "Wrong number of arguments.");
            init_lookup_table(0, _coords...);
            do_coords_checks();
            compute_total_size();

            total_size = 1;
            for(auto const& clist: coords)
                total_size *= clist.size();

            raw_data.resize(total_size);
        }

        template<size_t DIM, typename DATA_TYPE>
        template <typename... COORDS>
        lookup_table<DIM, DATA_TYPE>::lookup_table(COORDS... _coords, std::vector<DATA_TYPE>&& _raw_data){
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

        template<size_t DIM, typename DATA_TYPE>
        bool lookup_table<DIM, DATA_TYPE>::fill_at(const std::array<DATA_TYPE, DIM>& where, DATA_TYPE what){
            std::array<size_t, DIM> pos;
            bool success = where_to_pos(where, pos);
            if(success)
                raw_data[compute_index(pos)] = what;

            return success;
        }

        template<size_t DIM, typename DATA_TYPE>
        bool lookup_table<DIM, DATA_TYPE>::get_at(const std::array<DATA_TYPE, DIM>& where, DATA_TYPE& what){
            std::array<size_t, DIM> pos;
            bool success = where_to_pos(where, pos);
            if(success)
                what = raw_data[compute_index(pos)];

            return success;
        }

        template<size_t DIM, typename DATA_TYPE>
        template<typename COORDS>
        void lookup_table<DIM, DATA_TYPE>::init_lookup_table(size_t counter, COORDS _coords){
            coords[counter] = _coords;
        }
        template<size_t DIM, typename DATA_TYPE>
        template<typename COORDS, typename... REST>
        void lookup_table<DIM, DATA_TYPE>::init_lookup_table(size_t counter, COORDS _coords, REST... _other_coords){
            coords[counter] = (_coords);
            init_lookup_table(++counter, _other_coords...);
        }

        template<size_t DIM, typename DATA_TYPE>
        size_t lookup_table<DIM, DATA_TYPE>::get_dims(){
            return coords.size();
        }

        template<size_t DIM, typename DATA_TYPE>
        size_t lookup_table<DIM, DATA_TYPE>::get_data_size(){
            return raw_data.size();
        }

        template<size_t DIM, typename DATA_TYPE>
        std::vector<DATA_TYPE> lookup_table<DIM, DATA_TYPE>::get_coords(size_t dim){
            return coords.at(dim);
        }

        template<size_t DIM, typename DATA_TYPE>
        void lookup_table<DIM, DATA_TYPE>::do_coords_checks(){
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

        template<size_t DIM, typename DATA_TYPE>
        void lookup_table<DIM, DATA_TYPE>::compute_total_size(){
            total_size = 1;
            for(auto const& clist: coords)
                total_size *= clist.size();
        }


        template<size_t DIM, typename DATA_TYPE>
        size_t lookup_table<DIM, DATA_TYPE>::compute_index(std::array<size_t, DIM> pos){
            size_t index = 0;
            size_t mul = 1;
            for(auto it = pos.rbegin(); it != pos.rend();  it++){
                index += mul*(*it);
                mul *= coords[ std::distance(pos.begin(),it.base())-1].size();
            }
            return index;
        }

        template<size_t DIM, typename DATA_TYPE>
        bool lookup_table<DIM, DATA_TYPE>::where_to_pos(const std::array<DATA_TYPE, DIM>& where, std::array<size_t, DIM>& pos){
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

        template<size_t DIM, typename DATA_TYPE>
        void lookup_table<DIM, DATA_TYPE>::print_on_disk(std::string file_name, bool reloadable){
            std::ofstream of;
            of.open(file_name.c_str());

            if(!reloadable){

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
            }
            else{
                of << DIM << std::endl;
                for(auto& cc: coords){
                    of << cc.size() << std::endl;
                    for(auto ee: cc){
                        of << ee << " ";
                    }
                    of << std::endl;
                }
                for (auto dd: raw_data){
                    of << dd << " ";
                }
                of << std::endl;
            }

            of.close();
        }

        template<size_t DIM, typename DATA_TYPE>
        void lookup_table<DIM, DATA_TYPE>::read_from_disk(std::string file_name){
            std::ifstream ifile;
            ifile.open(file_name.c_str());

            size_t dim;
            ifile >> dim;
            if(dim != DIM)
                throw std::logic_error(std::string("Wrong file!"));

            size_t totsize = 1;
            for(size_t idim = 0; idim < DIM; idim++){
                size_t dimsize;
                ifile >> dimsize;
                totsize *= dimsize;
                coords[idim].resize(dimsize);
                for(size_t icc = 0; icc < dimsize; icc++){
                    DATA_TYPE dummy;
                    ifile >> dummy;
                    coords[idim][icc] = dummy;
                }

            }

            raw_data.resize(totsize);
            for(size_t iraw = 0; iraw < totsize; iraw++){
                DATA_TYPE dummy;
                ifile >> dummy;
                raw_data[iraw] = dummy;
            }

            total_size = totsize;

            ifile.close();
        }

        template<size_t DIM, typename DATA_TYPE>
        bool lookup_table<DIM, DATA_TYPE>::interp_at(const std::array<DATA_TYPE, DIM>& where, DATA_TYPE& what){
            for(size_t i = 0; i < DIM; i++){
                if(where[i] < coords[i].front() || where[i] > coords[i].back()){
                    err("Warning: coord " + std::to_string(where[i]) +
                        " outside table limits [" + std::to_string(coords[i].front()) + ":" + std::to_string(coords[i].back())  + "]!");
                    return false;
                }
            }
            std::array<DATA_TYPE, DIM> left_coord;
            std::array<DATA_TYPE, DIM> right_coord;
            std::array<size_t, DIM> left_idx;
            std::array<size_t, DIM> right_idx;
            std::array<DATA_TYPE, DIM> dist;

            for(size_t i = 0; i < DIM; i++){
                auto right = std::upper_bound(coords[i].begin(), coords[i].end(), where[i]);
                auto left = right - 1;
                left_coord[i] = *left;
                right_coord[i] = *right;
                left_idx[i] = std::distance(coords[i].begin(), left);
                right_idx[i] = std::distance(coords[i].begin(), right);
                dist[i] = right_coord[i] - left_coord[i];
            }

            DATA_TYPE hypercube_vol = 1.0;
            for(auto cc: dist){
                hypercube_vol*= cc;
            }
            DATA_TYPE inverse_hypercube_vol = 1.0/hypercube_vol;


            what = 0.0;

            std::array<size_t, DIM> pos;
            for (size_t idx = 0; idx < SIZE_2N; idx++){
                DATA_TYPE temp_mul = 1.0;
                for(size_t dim = 0; dim < DIM; dim++){
                    if(idx & (1 << dim)){
                        temp_mul *= (where[dim] - left_coord[dim]);
                        pos[dim] = right_idx[dim];
                    }
                    else{
                        temp_mul *= (right_coord[dim] - where[dim]);
                        pos[dim] = left_idx[dim];
                    }
                }
                what += raw_data[compute_index(pos)]*temp_mul;
            }

            what *= inverse_hypercube_vol;

            return true;
        }

        template<size_t DIM, typename DATA_TYPE>
        bool lookup_table<DIM, DATA_TYPE>::interp_at_one_dim(size_t which_dim, DATA_TYPE where, const std::array<size_t, DIM>& pos, DATA_TYPE& what){
            if(where < coords[which_dim].front() || where > coords[which_dim].back())
                return false;

            auto right = std::upper_bound(coords[which_dim].begin(), coords[which_dim].end(), where);
            auto left = right - 1;
            double left_coord = *left;
            double right_coord = *right;
            double dist = right_coord - left_coord;
            size_t left_idx = std::distance(coords[which_dim].begin(), left);
            size_t right_idx = std::distance(coords[which_dim].begin(), right);

            double left_coeff = (right_coord - where);
            double right_coeff = (where - left_coord);
            std::array<size_t, DIM> left_pos = pos;
            std::array<size_t, DIM> right_pos = pos;
            left_pos[which_dim] = left_idx;
            right_pos[which_dim] = right_idx;

            what = (raw_data[compute_index(left_pos)]*left_coeff + raw_data[compute_index(right_pos)]*right_coeff)/dist;

            return true;
        }
    }
}



#endif
