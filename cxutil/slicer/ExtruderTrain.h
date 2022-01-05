//Copyright (c) 2018 Ultimaker B.V.
//CuraEngine is released under the terms of the AGPLv3 or higher.
#include<iostream>
#ifndef EXTRUDER_TRAIN_H
#define EXTRUDER_TRAIN_H

namespace cxutil
{
    class Settings;
    class ExtruderTrain
    {
    public:
        /*
         * \brief Creates a new extruder. This extruder will have no settings but
         * inherits everything from the global settings.
         */
        ExtruderTrain(const size_t extruder_nr, Settings* parent_settings);
        ExtruderTrain(const ExtruderTrain& train);
        ~ExtruderTrain();
        /*
         * \brief The settings that this extruder overwrites.
         */
        Settings* settings;

        /*
         * \brief The position of this extruder.
         *
         * This may be used by g-code commands such as T to indicate to which
         * tool we must switch.
         */
        const size_t extruder_nr;
    };
}//namespace cxutil
#endif // EXTRUDER_TRAIN_H