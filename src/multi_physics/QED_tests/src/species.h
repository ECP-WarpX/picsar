#ifndef __SPECIES__
#define __SPECIES__

#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
#include <map>
#include <algorithm>

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
        void swap_particles(size_t index1, size_t index2);
        void remove_particle(size_t particle_index);
        void remove_particles_from_toend(size_t particle_index);
        void remove_particles_given_sorted_indices(const std::vector<size_t>& kill_index);

        virtual void push_positions(ttime dt) = 0;
        virtual void push_momenta(ttime dt) = 0;
        void calc_fields(em_field_function em_function, ttime tt);

        void add_simple_process(simple_process proc, int ID);
        void do_simple_process(int ID, ttime dt);
        void do_all_simple_processes(ttime dt);

        void add_process_with_destruction(process_with_destruction proc, int ID);
        void do_process_with_destruction(int ID, ttime dt);

        void print_on_disk(std::string prefix, int step_num) const;

        std::tuple<positions_list, momenta_list,
        em_field_list,  std::vector<double>> get_copy_of_all_data();
        positions_list get_copy_of_positions();
        momenta_list get_copy_of_momenta();

        std::tuple<positions_list&, momenta_list&,
        em_field_list&,  std::vector<double>&> get_ref_of_all_data();
        positions_list& get_ref_of_positions();
        momenta_list& get_ref_of_momenta();
        std::vector<double>& get_ref_of_optical_depth();

    protected:
        std::string name;
        double mass;
        double charge;

        positions_list pos;
        momenta_list mom;
        em_field_list fields;
        std::vector<double> optical_depth;

        std::map<int,simple_process> simple_processes;
        std::map<int,process_with_destruction> processes_with_destruction;

        std::string header = "species";

        void overwrite_particle(size_t isource, size_t idest);
    };
}

#endif
