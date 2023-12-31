//Copyright (c) 2018 Ultimaker B.V.
//CuraEngine is released under the terms of the AGPLv3 or higher.

#ifndef CXUTIL_VELOCITY_H
#define CXUTIL_VELOCITY_H

namespace cxutil
{

    /*
     * Represents a velocity in millimetres per second.
     *
     * This is a facade. It behaves like a double, only it can't be negative.
     */
    struct Velocity
    {
        /*
         * \brief Default constructor setting velocity to 0.
         */
        constexpr Velocity() : value(0.0) {};

        /*
         * \brief Casts a double to a Velocity instance.
         */
        constexpr Velocity(double value) : value(value) {};

        /*
         * \brief Casts the Temperature instance to a double.
         */
        constexpr operator double() const
        {
            return value;
        }

        /*
         * Some operators for arithmetic on velocities.
         */
        Velocity operator *(const Velocity& other) const
        {
            return Velocity(value * other.value);
        }
        template<typename E> Velocity operator *(const E& other) const
        {
            return Velocity(value * other);
        }
        Velocity operator /(const Velocity& other) const
        {
            return Velocity(value / other.value);
        }
        template<typename E> Velocity operator /(const E& other) const
        {
            return Velocity(value / other);
        }
        Velocity& operator *=(const Velocity& other)
        {
            value *= other.value;
            return *this;
        }
        template<typename E> Velocity& operator *=(const E& other)
        {
            value *= other;
            return *this;
        }
        Velocity& operator /=(const Velocity& other)
        {
            value /= other.value;
            return *this;
        }
        template<typename E> Velocity& operator /=(const E& other)
        {
            value /= other;
            return *this;
        }

        /*
         * \brief The actual temperature, as a double.
         */
        double value = 0;
    };

    using Acceleration = Velocity; //Use the same logic for acceleration variables.

}

#endif //CXUTIL_VELOCITY_H