#ifndef __SPECIES__
#define __SPECIES__

#include <string>
#include <iostream>
#include <fstream>

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

        void add_particle(position part_pos, momentum part_mom);

        virtual void push_positions(time dt) = 0;
        virtual void push_momenta(time dt) = 0;
        void calc_fields(em_field_function em_function, time tt);

        void print_on_disk(std::string prefix, int step_num) const;
    protected:
        std::string name;
        double mass;
        double charge;

        positions_list pos;
        momenta_list mom;
        em_field_list fields;

        std::string header = "species";
    };
}

#endif
