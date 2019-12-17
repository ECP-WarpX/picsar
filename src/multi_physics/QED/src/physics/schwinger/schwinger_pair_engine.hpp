#ifndef PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE
#define PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE

//This .hpp file contais the implementation of the
//Schwinger pair engine (as described in Gonoskov et al. PRE 92, 023305 2015
// and Banerjee et al. PRA 98, 032121 2018)

//Should be included by all the src files of the library
#include "../../qed_commons.h"

#include "../../containers/picsar_array.hpp"


#include "schwinger_pair_engine_core.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{

    template<
        typename RealType, class RandWrap,
        unit_system UnitSystem = unit_system::SI>
    class schwinger_pair_engine
    {
    public:
        //A random number generatator has to be passed by move.
        //The RNG can be ANY object implementing the functions
        //RealType unf (RealType a, RealType b)
        //and
        //RealType exp (RealType l)
        //The constructor can accept a lambda parameter.
        //It is ignored if the SI units option is selected
        schwinger_pair_engine(
            RandWrap&& t_rng, RealType t_lambda = static_cast<RealType>(1.0)):
            rng{t_rng}, lambda{t_lambda}{}

        //Setter for lambda
        RealType set_lambda(RealType t_lambda)
        {
            m_lambda = t_lambda;
        }

        //Getter for lambda
        RealType get_lambda() const
        {
            return m_lambda;
        }

        //This function determines if a pair has been generated in a given
        //cell
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        int generate_pairs_single(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType dx, const RealType dy, const RealType dz,
        const RealType dt)
        {
            return internal_generate_pairs_single<RealType, UnitSystem>(
                ex, ey, ez, bx, by, bz, dx, dy, dz, dt
                m_rng.unf(zero, one), m_lambda);
        }

        //This function determines how many pairs have been generated in a given
        //cell.
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        int generate_pairs_multiple(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType dx, const RealType dy, const RealType dz,
        const RealType dt)
        {
            return generate_pairs_multiple<RealType, UnitSystem>(
                ex, ey, ez, bx, by, bz, dx, dy, dz, dt,
                m_rng.unf(zero, one), m_lambda);
        }

        //This function computes 3 random numbers between 0 and 1.
        //They should be used to initialize the position of the pair.
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        containers::picsar_array<RealType, 3> get_new_pair_position()
        {
            auto res = picsar_array<RealType, 3> res;
            for(auto& val: res) val = rng.unf(zero, one);
            return res;
        }

        //Computes the pair production rate per unit time per unit volume
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType compute_schwinger_pair_production_rate(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz
        ) const
        {
            return compute_schwinger_pair_production_rate(
                ex, ey, ez, bx, by, bz, m_lambda);
        }

    private:

        RealType m_lambda;
        //The only requrement for the RNG is to be able to provide unf(a,b) and
        //exp(l)
        RandWrap m_rng;

        //Some handy constants
        constexpr const auto zero = static_cast<RealType>(0.0);
        constexpr const auto one = static_cast<RealType>(1.0);
        constexpr const auto two = static_cast<RealType>(2.0);

}
}
}

#endif //__PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE__
