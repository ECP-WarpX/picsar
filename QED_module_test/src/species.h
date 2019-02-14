#ifndef __SPECIES__
#define __SPECIES__

#include <string>
#include <iostream>
#include <fstream>
#include <tuple>

#include "commons.h"

namespace testbed{
    class species{
    public:
        species(std::string name);
        virtual ~species();
        std::string get_name() const;
        double get_mass() const;
        double get_charge() const;
        int get_num_particles() const;

        void add_particle(picsar::position part_pos, picsar::momentum part_mom);

        virtual void push_positions(picsar::ttime dt) = 0;
        virtual void push_momenta(picsar::ttime dt) = 0;
        void calc_fields(picsar::em_field_function em_function, picsar::ttime tt);

        void print_on_disk(std::string prefix, int step_num) const;

        std::tuple<picsar::positions_list, picsar::momenta_list,
        picsar::em_field_list,  std::vector<double>> get_copy_of_all_data();
        picsar::positions_list get_copy_of_positions();
        picsar::momenta_list get_copy_of_momenta();
    protected:
        std::string name;
        double mass;
        double charge;

        picsar::positions_list pos;
        picsar::momenta_list mom;
        picsar::em_field_list fields;
        std::vector<double> optical_depth;

        std::string header = "species";
    };
}

#endif
