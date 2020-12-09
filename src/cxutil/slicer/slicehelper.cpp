#include "cxutil/slicer/slicehelper.h"
#include "cxutil/math/linearAlg2D.h"

namespace cxutil
{
	SliceHelper::SliceHelper()
		:mesh(nullptr)
	{

	}

	SliceHelper::~SliceHelper()
	{

	}

	void SliceHelper::prepare(MeshObject* _mesh)
	{
		mesh = _mesh;
		buildMeshFaceHeightsRange(mesh, faceRanges);
	}

	SlicerSegment SliceHelper::project2D(const Point3& p0, const Point3& p1, const Point3& p2, const coord_t z)
	{
		SlicerSegment seg;

		seg.start.X = interpolate(z, p0.z, p1.z, p0.x, p1.x);
		seg.start.Y = interpolate(z, p0.z, p1.z, p0.y, p1.y);
		seg.end.X = interpolate(z, p0.z, p2.z, p0.x, p2.x);
		seg.end.Y = interpolate(z, p0.z, p2.z, p0.y, p2.y);

		return seg;
	}

	void SliceHelper::buildMeshFaceHeightsRange(const MeshObject* mesh, std::vector<Point2>& heightRanges)
	{
		int faceSize = (int)mesh->faces.size();
		if (faceSize > 0) heightRanges.resize(faceSize);
		for (int i = 0; i < faceSize; ++i)
		{
			const MeshFace& face = mesh->faces[i];
			const MeshVertex& v0 = mesh->vertices[face.vertex_index[0]];
			const MeshVertex& v1 = mesh->vertices[face.vertex_index[1]];
			const MeshVertex& v2 = mesh->vertices[face.vertex_index[2]];

			// get all vertices represented as 3D point
			Point3 p0 = v0.p;
			Point3 p1 = v1.p;
			Point3 p2 = v2.p;

			// find the minimum and maximum z point		
			int32_t minZ = p0.z;
			if (p1.z < minZ)
			{
				minZ = p1.z;
			}
			if (p2.z < minZ)
			{
				minZ = p2.z;
			}

			int32_t maxZ = p0.z;
			if (p1.z > maxZ)
			{
				maxZ = p1.z;
			}
			if (p2.z > maxZ)
			{
				maxZ = p2.z;
			}

			heightRanges.at(i) = Point2(minZ, maxZ);
		}
	}

	void SliceHelper::sliceOneLayer(int z,
		std::vector<SlicerSegment>& segments, std::unordered_map<int, int>& face_idx_to_segment_idx)
	{
		segments.reserve(100);

		int faceNum = mesh->faces.size();
		// loop over all mesh faces
		for (int faceIdx = 0; faceIdx < faceNum; ++faceIdx)
		{
			if ((z < faceRanges[faceIdx].x) || (z > faceRanges[faceIdx].y))
				continue;

			// get all vertices per face
			const MeshFace& face = mesh->faces[faceIdx];
			const MeshVertex& v0 = mesh->vertices[face.vertex_index[0]];
			const MeshVertex& v1 = mesh->vertices[face.vertex_index[1]];
			const MeshVertex& v2 = mesh->vertices[face.vertex_index[2]];

			// get all vertices represented as 3D point
			Point3 p0 = v0.p;
			Point3 p1 = v1.p;
			Point3 p2 = v2.p;

			SlicerSegment s;
			s.endVertex = nullptr;
			int end_edge_idx = -1;

			if (p0.z < z && p1.z >= z && p2.z >= z)
			{
				s = project2D(p0, p2, p1, z);
				end_edge_idx = 0;
				if (p1.z == z)
				{
					s.endVertex = &v1;
				}
			}
			else if (p0.z > z && p1.z < z && p2.z < z)
			{
				s = project2D(p0, p1, p2, z);
				end_edge_idx = 2;
			}
			else if (p1.z < z && p0.z >= z && p2.z >= z)
			{
				s = project2D(p1, p0, p2, z);
				end_edge_idx = 1;
				if (p2.z == z)
				{
					s.endVertex = &v2;
				}
			}
			else if (p1.z > z && p0.z < z && p2.z < z)
			{
				s = project2D(p1, p2, p0, z);
				end_edge_idx = 0;
			}
			else if (p2.z < z && p1.z >= z && p0.z >= z)
			{
				s = project2D(p2, p1, p0, z);
				end_edge_idx = 2;
				if (p0.z == z)
				{
					s.endVertex = &v0;
				}
			}
			else if (p2.z > z && p1.z < z && p0.z < z)
			{
				s = project2D(p2, p0, p1, z);
				end_edge_idx = 1;
			}
			else
			{
				//Not all cases create a segment, because a point of a face could create just a dot, and two touching faces
				//  on the slice would create two segments
				continue;
			}

			// store the segments per layer
			face_idx_to_segment_idx.insert(std::make_pair(faceIdx, segments.size()));
			s.faceIndex = faceIdx;
			s.endOtherFaceIdx = face.connected_face_index[end_edge_idx];
			segments.push_back(s);
		}
	}
}