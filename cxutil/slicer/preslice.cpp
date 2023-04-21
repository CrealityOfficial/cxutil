#include "preslice.h"
#include "cxutil/input/groupinput.h"
#include "cxutil/input/dlpinput.h"
#include "trimesh2/TriMesh.h"
#include <float.h>
#include "cxutil/settings/AdaptiveLayerHeights.h"
namespace cxutil
{
    void buildSliceInfos(const DLPInput& input, std::vector<int>& z)
    {
        const std::vector<MeshObjectPtr>& meshes = input.Meshes;
        const DLPParam& param = input.Param;

        int slice_layer_count = 0;
        const coord_t initial_layer_thickness = param.initial_layer_thickness;
        const coord_t layer_thickness = param.layer_thickness;

        AABB3D box = input.box();
        slice_layer_count = (box.max.z - initial_layer_thickness) / layer_thickness + 1;

        if (slice_layer_count > 0)
        {
            z.resize(slice_layer_count, 0);

            z[0] = std::max(0LL, initial_layer_thickness - layer_thickness);
            coord_t adjusted_layer_offset = initial_layer_thickness;

            z[0] = initial_layer_thickness / 2;
            adjusted_layer_offset = initial_layer_thickness + (layer_thickness / 2);

            // define all layer z positions (depending on slicing mode, see above)
            for (int layer_nr = 1; layer_nr < slice_layer_count; layer_nr++)
            {
                z[layer_nr] = adjusted_layer_offset + (layer_thickness * (layer_nr - 1));
            }
        }
    }
}
