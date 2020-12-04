#include "cxutil/input/dlpdata.h"

namespace cxutil
{
	DLPData::DLPData()
	{
	}
	DLPData::~DLPData()
	{
		for (std::vector<DLPmeshs>::iterator it= m_dlpmeshsgroup.dlpmeshsgroup.begin(); it!= m_dlpmeshsgroup.dlpmeshsgroup.end();it++)
		{
			for (std::vector<DLPmesh>::iterator it2 = it->dlpmeshs.begin();it2!= it->dlpmeshs.end();it2++)
			{
				for (std::vector<DLPLayer>::iterator it3= it2->layers.begin(); it3!=it2->layers.end();it3++)
				{
					for (std::vector< ClipperLib::PolyTree*>::iterator it4=it3->parts.begin();it4!=it3->parts.end();it4++)
					{
						delete *it4;
					}
				}
			}
		}
		
	}
}
