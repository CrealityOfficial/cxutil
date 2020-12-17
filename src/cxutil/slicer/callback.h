#ifndef CXUTIL_CALLBACK_1607502965747_H
#define CXUTIL_CALLBACK_1607502965747_H
#include "cxutil/math/AABB3D.h"

namespace cxutil
{
	class SliceCallback
	{
	public:
		virtual ~SliceCallback() {}

		virtual void onSceneBox(const AABB3D& box3) = 0;
		
		virtual void onLayerCount(int layer) = 0;

		virtual void onFilamentLen(double len) = 0;
		virtual void onPrintTime(int time) = 0;
	};
}

#endif // CXUTIL_CALLBACK_1607502965747_H