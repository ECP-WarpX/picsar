#include "species.h"

using namespace testbed;

species::species(std::string name):name(name){}

species::~species(){}

std::string species::get_name() const{
    return name;
}

double species::get_mass() const{
    return mass;
}

double species::get_charge() const{
    return charge;
}

int species::get_num_particles() const{
    return pos[0].size();
}

void species::add_simple_process(simple_process proc, int ID){
    if (simple_processes.find(ID) != simple_processes.end())
        throw std::logic_error("ID already assigned");
    else
        simple_processes.insert(std::pair<int, simple_process>(ID, proc));
}

void species::do_all_simple_processes(ttime dt){
    for(auto& proc : simple_processes){
        proc.second(pos, mom,  fields, optical_depth, flag, mass, charge, dt);
    }
}

void species::do_simple_process(int ID, ttime dt){
    simple_processes[ID](pos, mom,  fields, optical_depth, flag, mass, charge, dt);
}

void species::add_process_with_destruction(process_with_destruction proc, int ID){
    if (processes_with_destruction.find(ID) != processes_with_destruction.end())
        throw std::logic_error("ID already assigned");
    else
        processes_with_destruction.insert(std::pair<int, process_with_destruction>(ID, proc));
}
void species::do_process_with_destruction(int ID, ttime dt){
    std::vector<size_t> to_be_destroyed;
    for(size_t i = 0; i < pos[0].size(); i++){
        if(processes_with_destruction[ID](
            {pos[0][i], pos[1][i], pos[2][i]},
            {mom[0][i], mom[1][i], mom[2][i]},
            {fields[0][i], fields[1][i], fields[2][i],
            fields[3][i], fields[4][i], fields[5][i]},
            optical_depth[i], flag[i], mass, charge, dt))
                to_be_destroyed.push_back(i);
    }
    remove_particles_given_sorted_indices(to_be_destroyed);


}


void species::add_particle(position part_pos, momentum part_mom){
    for(int i = 0; i < 3; i++){
        pos[i].push_back(part_pos[i]);
        mom[i].push_back(part_mom[i]);
    }
    for(int i = 0; i < 6; i++){
        fields[i].push_back(0);
    }
    optical_depth.push_back(0.0);
    flag.push_back(false);

}

void species::swap_particles(size_t index1, size_t index2){
    for(int i = 0; i < 3; i++){
        std::swap(pos[index1], pos[index2]);
        std::swap(mom[index1], mom[index2]);
    }
    for(int i = 0; i < 6; i++){
        std::swap(fields[index1], fields[index2]);
    }
    std::swap(optical_depth[index1], optical_depth[index2]);
    std::swap(flag[index1], flag[index2]);
}

void species::remove_particle(size_t particle_index){
    for (auto& p : pos)
        p.erase(p.begin() + particle_index);

    for (auto& m : mom)
        m.erase(m.begin() + particle_index);

    for (auto& f: fields)
        f.erase(f.begin() + particle_index);

    optical_depth.erase(optical_depth.begin() + particle_index);
    flag.erase(flag.begin() + particle_index);
}

void species::remove_particles_from_toend(size_t particle_index){
    for (auto& p : pos)
        p.erase(p.begin() + particle_index, p.end());

    for (auto& m : mom)
        m.erase(m.begin() + particle_index, m.end());

    for (auto& f: fields)
        f.erase(f.begin() + particle_index, f.end());

    optical_depth.erase(optical_depth.begin() + particle_index, optical_depth.end());
    flag.erase(flag.begin() + particle_index, flag.end());
}

void species::overwrite_particle(size_t isource, size_t idest){
    for(int i = 0; i < 3; i++){
        pos[i][idest] = pos[i][isource];
        mom[i][idest] = mom[i][isource];
    }
    for(int i = 0; i < 6; i++){
         fields[i][idest] = fields[i][isource];
    }
     optical_depth[idest] = optical_depth[isource];
     flag[idest] = flag[isource];
}

void species::remove_particles_given_sorted_indices(const std::vector<size_t>& kill_index){
    size_t last = pos[0].size()-1;
    for(auto ki = kill_index.begin(); ki != kill_index.end() ; ki++){
        size_t ll = last;
        while(ll > *ki){
            if(!binary_search( kill_index.begin(), kill_index.end(), ll)){
                overwrite_particle(ll, *ki);
                ll--;
                break;
            }
            ll--;
        }
        last = ll;
    }
    remove_particles_from_toend(pos[0].size()-kill_index.size());
}

void species::calc_fields(em_field_function em_function, ttime tt){
    int num_particles = pos[0].size();

    for (int i = 0; i < num_particles; i++){
        em_field EB = em_function({pos[0][i], pos[1][i], pos[2][i]}, tt);
        for (size_t j = 0; j < EB.size(); j++){
            fields[j][i] = EB[j];
        }
    }
}

void species::print_on_disk(std::string prefix, int step_num) const{
    std::ofstream out_file(prefix + "_" + name + "_" + std::to_string(step_num) +".dat");

    out_file << "# " << header << std::endl;

    int num_particles = pos[0].size();
    out_file << "#N = " << num_particles << std::endl;

    out_file << "#x\ty\tz\tpx\tpy\tpz\tEx\tEy\tEz\tBx\tBy\tBz\tOptDepth" << std::endl;

    for (int i = 0; i < num_particles; i++){
        out_file << pos[0][i] << " " << pos[1][i] << " " << pos[2][i] << " ";
        out_file << mom[0][i] << " " << mom[1][i] << " " << mom[2][i] << " ";
        out_file << fields[0][i] << " " << fields[1][i] << " " << fields[2][i] << " ";
        out_file << fields[3][i] << " " << fields[4][i] << " " << fields[5][i] << " ";
        out_file <<  optical_depth[i] << " ";
        out_file <<  flag[i] << std::endl;
    }

    out_file.close();
}

std::tuple<positions_list, momenta_list,
em_field_list,  std::vector<double>, std::vector<bool>>
species::get_copy_of_all_data(){
    return std::tuple<positions_list, momenta_list,
    em_field_list,  std::vector<double>, std::vector<bool> >
        {pos, mom, fields, optical_depth, flag};
}

positions_list species::get_copy_of_positions(){
    return pos;
}

momenta_list species::get_copy_of_momenta(){
    return mom;
}

std::tuple<positions_list&, momenta_list&,
em_field_list&,  std::vector<double>&, std::vector<bool>&>
species::get_ref_of_all_data(){
    return std::tuple<positions_list&, momenta_list&,
    em_field_list&,  std::vector<double>&, std::vector<bool>&>
        {pos, mom, fields, optical_depth, flag};
}

positions_list& species::get_ref_of_positions(){
    return pos;
}
momenta_list& species::get_ref_of_momenta(){
    return mom;
}

std::vector<double>& species::get_ref_of_optical_depth(){
    return optical_depth;
}
