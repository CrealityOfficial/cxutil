#include "cxutil/slicer/slicepolygonbuilder.h"
#include "cxutil/math/utilconst.h"
#include "cxutil/slicer/slicehelper.h"
#include "trimesh2/TriMesh.h"

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

	coord_t interpolate_my(const coord_t x, const coord_t x0, const coord_t x1, const coord_t y0, const coord_t y1)
	{
		const coord_t dx_01 = x1 - x0;
		coord_t num = (y1 - y0) * (x - x0);
		num += num > 0 ? dx_01 / 2 : -dx_01 / 2; // add in offset to round result
		return y0 + num / dx_01;
	}

	void project2D(const Point3& p0, const Point3& p1, const Point3& p2, const coord_t z, SlicerSegment& seg)
	{
		if (seg.faceIndex == -1)
		{
			seg.start.X = interpolate_my(z, p0.z, p1.z, p0.x, p1.x);
			seg.start.Y = interpolate_my(z, p0.z, p1.z, p0.y, p1.y);
		}
		else
		{
			seg.start.X = seg.end.X;
			seg.start.Y = seg.end.Y;
		}
		seg.end.X = interpolate_my(z, p0.z, p2.z, p0.x, p2.x);
		seg.end.Y = interpolate_my(z, p0.z, p2.z, p0.y, p2.y);

	}

	bool zCrossFace(trimesh::TriMesh* mesh, cxutil::SliceHelper* helper, int faceId, int z, SlicerSegment& s)
	{
		const trimesh::TriMesh::Face& face = mesh->faces[faceId];

		// get all vertices represented as 3D point
		Point3 p0 = Point3(MM2INT(mesh->vertices[face[0]].x), MM2INT(mesh->vertices[face[0]].y), MM2INT(mesh->vertices[face[0]].z));
		Point3 p1 = Point3(MM2INT(mesh->vertices[face[1]].x), MM2INT(mesh->vertices[face[1]].y), MM2INT(mesh->vertices[face[1]].z));
		Point3 p2 = Point3(MM2INT(mesh->vertices[face[2]].x), MM2INT(mesh->vertices[face[2]].y), MM2INT(mesh->vertices[face[2]].z));

		int end_edge_idx = -1;

		if (p0.z < z && p1.z >= z && p2.z >= z)
		{
			project2D(p0, p2, p1, z, s);
			if (p1.z == z)
			{
				if (s.end != s.start)
				{
					s.end.X = p1.x;
					s.end.Y = p1.y;
					s.endOtherFaceIdx = -1;
					return true;
				}
				else
				{
					return false;
				}
			}
			else
			{
				end_edge_idx = 0;
			}
		}
		else if (p0.z > z && p1.z < z && p2.z < z)
		{
			project2D(p0, p1, p2, z, s);
			end_edge_idx = 2;
		}
		else if (p1.z < z && p0.z >= z && p2.z >= z)
		{
			project2D(p1, p0, p2, z, s);
			if (p2.z == z)
			{
				if (s.end != s.start)
				{
					s.end.X = p2.x;
					s.end.Y = p2.y;
					s.endOtherFaceIdx = -1;
					return true;
				}
				else
				{
					return false;
				}
			}
			else
			{
				end_edge_idx = 1;
			}
		}
		else if (p1.z > z && p0.z < z && p2.z < z)
		{
			project2D(p1, p2, p0, z, s);
			end_edge_idx = 0;
		}
		else if (p2.z < z && p1.z >= z && p0.z >= z)
		{
			project2D(p2, p1, p0, z, s);
			if (p0.z == z)
			{
				if (s.end != s.start)
				{
					s.end.X = p0.x;
					s.end.Y = p0.y;
					s.endOtherFaceIdx = -1;
					return true;
				}
				else
				{
					return false;
				}
			}
			else
			{
				end_edge_idx = 2;
			}
		}
		else if (p2.z > z && p1.z < z && p0.z < z)
		{
			project2D(p2, p0, p1, z, s);
			end_edge_idx = 1;
		}
		else
		{
			return false;
			//Not all cases create a segment, because a point of a face could create just a dot, and two touching faces
			//  on the slice would create two segments
		}

		s.faceIndex = faceId;
		if (end_edge_idx != -1)
		{
			s.endOtherFaceIdx = helper->faces[faceId].connected_face_index[end_edge_idx];
		}
		return true;
	}

	void SlicePolygonBuilder::sliceOneLayer_dst(cxutil::SliceHelper* helper, int z, Polygons* polygons, Polygons* openPolygons)
	{
		trimesh::TriMesh* meshSrc = helper->getMeshSrc();
		std::vector<Point2>* faceRanges = helper->getFaceRanges();
		int faceNum = meshSrc->faces.size();

		if (faceNum > 0) visitFlags.resize(faceNum, false);
		for (int i = 0; i < faceNum; ++i)
		{
			if ((z < faceRanges->at(i).x) || (z > faceRanges->at(i).y))
				continue;
			if (visitFlags[i]) continue;
			bool bFirstSegment = true;
			bool bClosed = false;
			SlicerSegment s;
			ClipperLib::Path poly;
			poly.reserve(100);
			for (int faceid = i; faceid != -1; )
			{
				if (visitFlags[faceid]) //开轮廓
				{
					break;
				}
				if (zCrossFace(meshSrc, helper, faceid, z, s))
				{
					if (bFirstSegment) poly.push_back(s.start);
					poly.push_back(s.end);
					visitFlags[faceid] = true;

					faceid = s.endOtherFaceIdx;
					if (bFirstSegment) bFirstSegment = false;
				}
				else
				{
					break; //开轮廓
				}

				if (faceid == i) //闭轮廓
				{
					bClosed = true;
					break;
				}
			}
			if (bClosed)
			{
				ClipperLib::CleanPolygon(poly);
				polygons->add(poly);
			}
			else if (!poly.empty())
			{////开轮廓不能在这里清理轮廓，否则在后续开轮廓闭合时会出错
				openPolygons->add(poly);
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

    void tryConnectPath(Polygons& open_polylines, int cell_size, bool needReverse = false)
    {
        //获取开多边形的起始点
        std::vector<SlicerSegment> openSegs;
        for (const ClipperLib::Path& path : open_polylines.paths)
        {
            if (path.size() < 2)
                continue;
            SlicerSegment seg;
            seg.start = path.front();
            seg.end = path.back();
            openSegs.push_back(seg);
        }

        for (size_t i = 0; i < openSegs.size(); i++)
        {
            SlicerSegment& seg1 = openSegs.at(i);

            for (size_t j = 0; j < openSegs.size(); j++)
            {
                SlicerSegment& seg2 = openSegs.at(j);
                if (i == j || seg2.addedToPolygon || open_polylines.paths[i].empty())
                    continue;

                Point diff = seg1.end - seg2.start;
                bool bdiff = shorterThen(diff, cell_size);
                if (bdiff)
                {
                    open_polylines.paths[i].insert(open_polylines.paths[i].end(), open_polylines.paths[j].begin(), open_polylines.paths[j].end());
                    open_polylines.paths[j].clear();
                    seg2.addedToPolygon = true;
                    seg1.end = open_polylines.paths[i].back();
                    break;
                }
                else if (needReverse) //处理反序
                {
                    Point diffR1 = seg1.start - seg2.end;
                    Point diffR2 = seg1.end - seg2.end;
                    Point diffR3 = seg1.start - seg2.start;
                    bool bdiffR1 = shorterThen(diffR1, cell_size);
                    bool bdiffR2 = shorterThen(diffR2, cell_size);
                    bool bdiffR3 = shorterThen(diffR3, cell_size);
                    if (bdiffR1)
                    {
                        open_polylines.paths[i].insert(open_polylines.paths[i].begin(), open_polylines.paths[j].begin(), open_polylines.paths[j].end());
                        open_polylines.paths[j].clear();
                        seg2.addedToPolygon = true;
                        seg1.end = open_polylines.paths[i].back();
                        break;
                    }
                    else if (bdiffR2)
                    {
                        if (open_polylines.paths[i].size() > open_polylines.paths[j].size())
                        {
                            std::reverse(open_polylines.paths[j].begin(), open_polylines.paths[j].end());
                            open_polylines.paths[i].insert(open_polylines.paths[i].end(), open_polylines.paths[j].begin(), open_polylines.paths[j].end());
                        }
                        else
                        {
                            std::reverse(open_polylines.paths[i].begin(), open_polylines.paths[i].end());
                            open_polylines.paths[i].insert(open_polylines.paths[i].begin(), open_polylines.paths[j].begin(), open_polylines.paths[j].end());
                        }

                        open_polylines.paths[j].clear();
                        seg2.addedToPolygon = true;
                        seg1.end = open_polylines.paths[i].back();
                        break;
                    }
                    else if (bdiffR3)
                    {
                        if (open_polylines.paths[i].size() > open_polylines.paths[j].size())
                        {
                            std::reverse(open_polylines.paths[j].begin(), open_polylines.paths[j].end());
                            open_polylines.paths[i].insert(open_polylines.paths[i].begin(), open_polylines.paths[j].begin(), open_polylines.paths[j].end());
                        }
                        else
                        {
                            std::reverse(open_polylines.paths[i].begin(), open_polylines.paths[i].end());
                            open_polylines.paths[i].insert(open_polylines.paths[i].end(), open_polylines.paths[j].begin(), open_polylines.paths[j].end());
                        }

                        open_polylines.paths[j].clear();
                        seg2.addedToPolygon = true;
                        seg1.end = open_polylines.paths[i].back();
                        break;
                    }
                }
            }
        }

        for (auto it = open_polylines.paths.begin(); it != open_polylines.paths.end(); )
        {
            if (it->empty())
            {
                it = open_polylines.paths.erase(it);
            }
            else
                it++;
        }
    }


    ClipperLib::IntPoint Vector(ClipperLib::IntPoint a, ClipperLib::IntPoint b)   //向量ab
    {
        return{ b.X - a.X, b.Y - a.Y };
    }

    double dot(ClipperLib::IntPoint A, ClipperLib::IntPoint B, ClipperLib::IntPoint P)      //向量的内积
    {
        ClipperLib::IntPoint AB = Vector(A, B);
        ClipperLib::IntPoint AP = Vector(A, P);
        return AB.X* AP.X + AB.Y * AP.Y;
    }

    double cross(ClipperLib::IntPoint A, ClipperLib::IntPoint B, ClipperLib::IntPoint P)  //向量的外积
    {
        ClipperLib::IntPoint AB = Vector(A, B);
        ClipperLib::IntPoint AP = Vector(A, P);
        return AB.X* AP.Y - AB.Y * AP.X;
    }

    double dis2(ClipperLib::IntPoint a, ClipperLib::IntPoint b)           //两点间的距离的平方
    {
        return (b.X - a.X) * (b.X - a.X) + (b.Y - a.Y) * (b.Y - a.Y);

    }

    int dir(ClipperLib::IntPoint A, ClipperLib::IntPoint B, ClipperLib::IntPoint P)     //点与线段方位判定
    {
        if (cross(A, B, P) > 0)  return -1;
        else if (cross(A, B, P) < 0) return 1;
        else if (dot(A, B, P) < 0) return -2;
        else if (dot(A, B, P) >= 0)
        {
            if (dis2(A, B) < dis2(A, P)) return 2;
            else return 0;
        }
    }

    double disLine(ClipperLib::IntPoint A, ClipperLib::IntPoint B, ClipperLib::IntPoint P)    //点P到直线AB的距离
    {
        return fabs(cross(A, B, P)) / sqrt(dis2(A, B));
    }

    bool pointInSegment(ClipperLib::IntPoint & A1, ClipperLib::IntPoint & A2, ClipperLib::IntPoint & p)
    {
        //判断交点有效
        if ((p.X >= A1.X && p.X <= A2.X) || (p.X <= A1.X && p.X >= A2.X))
        {
            if ((p.Y >= A1.Y && p.Y <= A2.Y) || (p.Y <= A1.Y && p.Y >= A2.Y))
                return true;
        }

        return false;
    }

    void insertPoint(ClipperLib::Path & path, ClipperLib::IntPoint & p)
    {
        for (size_t i = 0; i < path.size() - 1; i++)
        {
            if (pointInSegment(path[i], path[i + 1], p))
            {
                path.insert(path.begin() + i + 1, p);
                return;
            }
        }
    }

    struct InsertData
    {
        int index1;
        int index2;
        ClipperLib::IntPoint p;
    };

    void tryConnectPaths(Polygons & open_polylines, int cell_size, ClipperLib::Path& intersectPoints,bool needReverse = false)
    {

        //Polygons result = open_polylines.intersection_open();
        int test = 0;


        if (open_polylines.paths.size() > 1)
        {
            for (size_t i = 0; i < open_polylines.paths.size() - 1; i++, i++)
            {
                ClipperLib::Path& path1 = open_polylines.paths[i];
                ClipperLib::Path& path2 = open_polylines.paths[i + 1];

                //插入交点
                ClipperLib::Paths paths1;
                ClipperLib::Paths paths2;
                paths1.reserve(path1.size());
                paths2.reserve(path2.size());

                std::vector<InsertData> insertDataS;
                for (size_t j = 0; j < path1.size() - 1; j++)
                {
                    ClipperLib::IntPoint& A1 = path1[j];
                    ClipperLib::IntPoint& A2 = path1[j + 1];

                    ClipperLib::Path path;
                    path.push_back(A1);
                    path.push_back(A2);
                    paths1.push_back(path);

                    ClipperLib::Path::iterator iter = path2.begin();
                    for (size_t k = 0; k < path2.size() - 1; k++)
                    {
                        ClipperLib::IntPoint& B1 = path2[k];
                        ClipperLib::IntPoint& B2 = path2[k + 1];
                        if (dir(A1, A2, B1) * dir(A1, A2, B2) <= 0 && dir(B1, B2, A1) * dir(B1, B2, A2) <= 0)
                        {//判断有无交点
                            double t = disLine(A1, A2, B1) / (disLine(A1, A2, B1) + disLine(A1, A2, B2));

                            ClipperLib::IntPoint B1B2 = Vector(B1, B2);
                            ClipperLib::IntPoint inter = { B1.X + (ClipperLib::cInt)(B1B2.X * t), B1.Y + (ClipperLib::cInt)(B1B2.Y * t) };

                            if (pointInSegment(A1, A2, inter) && pointInSegment(B1, B2, inter))
                            {
                                InsertData insertData;
                                insertData.index1 = j;
                                insertData.index2 = k;
                                insertData.p = inter;
                                insertDataS.push_back(insertData);
                            }
                        }
                    }
                }

                for (size_t k = 0; k < path2.size() - 1; k++)
                {
                    ClipperLib::IntPoint& B1 = path2[k];
                    ClipperLib::IntPoint& B2 = path2[k + 1];

                    ClipperLib::Path path;
                    path.push_back(B1);
                    path.push_back(B2);
                    paths2.push_back(path);
                }

                for (size_t i = 0; i < insertDataS.size(); i++)
                {
                    ClipperLib::Path& path1 = paths1[insertDataS[i].index1];
                    ClipperLib::Path& path2 = paths2[insertDataS[i].index2];
                    ClipperLib::IntPoint& p = insertDataS[i].p;
                    insertPoint(path1, p);
                    insertPoint(path2, p);

                    intersectPoints.push_back(p);
                }

                ClipperLib::Path path11;
                ClipperLib::Path path22;


            }
        }
    }

    void SlicePolygonBuilder::connectOpenPolylines(Polygons & open_polylines,ClipperLib::Path& intersectPoints)
    {
        if (open_polylines.empty())
        {
            return;
        }

        constexpr bool allow_reverse = false;
        // Search a bit fewer cells but at cost of covering more area.
        // Since acceptance area is small to start with, the extra is unlikely to hurt much.
        constexpr coord_t cell_size = largest_neglected_gap_first_phase * 2;

        tryConnectPath(open_polylines, largest_neglected_gap_first_phase);

        //再次检测
        if (open_polylines.paths.size())
        {
            tryConnectPath(open_polylines, cell_size, true);
        }

        tryConnectPaths(open_polylines, cell_size, intersectPoints, true);
    }
}