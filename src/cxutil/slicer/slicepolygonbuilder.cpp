#include "cxutil/slicer/slicepolygonbuilder.h"
#include "cxutil/math/utilconst.h"

#include <fstream>

namespace cxutil
{
	SlicePolygonBuilder::SlicePolygonBuilder()
	{

	}

	SlicePolygonBuilder::~SlicePolygonBuilder()
	{

	}

	void SlicePolygonBuilder::makePolygon(Polygons* polygons, Polygons* openPolygons)
	{
		size_t segment_size = segments.size();
		if (segment_size > 0) visitFlags.resize(segment_size, false);
		for (size_t start_segment_idx = 0; start_segment_idx < segment_size; start_segment_idx++)
		{
			if (!visitFlags[start_segment_idx])
			{
				ClipperLib::Path poly;
				poly.push_back(segments[start_segment_idx].start);

				bool closed = false;
				for (int segment_idx = start_segment_idx; segment_idx != -1; )
				{
					const SlicerSegment& segment = segments[segment_idx];
					poly.push_back(segment.end);

					visitFlags[segment_idx] = true;
					segment_idx = getNextSegmentIdx(segment, start_segment_idx);
					if (segment_idx == static_cast<int>(start_segment_idx))
					{ // polyon is closed
						closed = true;
						break;
					}
				}

				if (closed)
				{
					polygons->add(poly);
				}
				else
				{
					// polygon couldn't be closed
					openPolygons->add(poly);
				}
			}
		}
	}

	int SlicePolygonBuilder::tryFaceNextSegmentIdx(const SlicerSegment& segment,
		const int face_idx, const size_t start_segment_idx)
	{
		decltype(face_idx_to_segment_idx.begin()) it;
		auto it_end = face_idx_to_segment_idx.end();
		it = face_idx_to_segment_idx.find(face_idx);
		if (it != it_end)
		{
			const int segment_idx = (*it).second;
			Point p1 = segments[segment_idx].start;
			Point diff = segment.end - p1;
			if (shorterThen(diff, largest_neglected_gap_first_phase))
			{
				if (segment_idx == static_cast<int>(start_segment_idx))
				{
					return start_segment_idx;
				}
				if (visitFlags[segment_idx])
				{
					return -1;
				}
				return segment_idx;
			}
		}

		return -1;
	}

	int SlicePolygonBuilder::getNextSegmentIdx(const SlicerSegment& segment, const size_t start_segment_idx)
	{
		int next_segment_idx = -1;

		const bool segment_ended_at_edge = segment.endVertex == nullptr;
		if (segment_ended_at_edge)
		{
			const int face_to_try = segment.endOtherFaceIdx;
			if (face_to_try == -1)
			{
				return -1;
			}
			return tryFaceNextSegmentIdx(segment, face_to_try, start_segment_idx);
		}
		else
		{
			// segment ended at vertex

			const std::vector<uint32_t>& faces_to_try = segment.endVertex->connected_faces;
			for (int face_to_try : faces_to_try)
			{
				const int result_segment_idx =
					tryFaceNextSegmentIdx(segment, face_to_try, start_segment_idx);
				if (result_segment_idx == static_cast<int>(start_segment_idx))
				{
					return start_segment_idx;
				}
				else if (result_segment_idx != -1)
				{
					// not immediately returned since we might still encounter the start_segment_idx
					next_segment_idx = result_segment_idx;
				}
			}
		}

		return next_segment_idx;
	}

	void SlicePolygonBuilder::write(std::fstream& out)
	{
		for (size_t j = 0; j < segments.size(); ++j)
		{
			const SlicerSegment& segment = segments.at(j);
			out << segment.start << " " << segment.endOtherFaceIdx << std::endl;
			out << segment.start << " " << segment.end << std::endl;
		}
		for (const std::pair<int, int>& p : face_idx_to_segment_idx)
		{
			out << p.first << " " << p.second << std::endl;
		}
	}
}