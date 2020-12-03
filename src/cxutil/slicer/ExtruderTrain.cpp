/** Copyright (C) 2016 Ultimaker - Released under terms of the AGPLv3 License */
#include "ExtruderTrain.h"
#include "cxutil/settings/Settings.h"

namespace cxutil
{
    ExtruderTrain::ExtruderTrain(const size_t extruder_nr, Settings* parent_settings)
        : extruder_nr(extruder_nr)
    {
        settings = new Settings(parent_settings);
    }

    ExtruderTrain::~ExtruderTrain()
    {
        delete settings;
    }

}//namespace cxutil