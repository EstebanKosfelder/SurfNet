using System.ComponentModel.DataAnnotations;
using System.Formats.Asn1;
using System.Net.Http.Headers;
using System.Runtime;
using System.Security.Cryptography.X509Certificates;
using TriangleNet;
using TriangleNet.Geometry;
using SurfNet.Meshing;
using SurfNet.Topology;
using SurfNet.Topology.DCEL;
using TriangleNet.Meshing;

namespace SurfNet
{
    using static DebugLog;

    using HVertex = TriangleNet.Topology.DCEL.Vertex;
    using TVertex = TriangleNet.Geometry.Vertex;
    using FT = double;
    public class BasicInput
    {
        public List<BasicFace> InnerFaces = new List<BasicFace>();

        private Dictionary<(int a, int b), TrianglesOnEdge> halfedges = new Dictionary<(int a, int b), TrianglesOnEdge>();

        public BasicInput()
        { }

        public List<Point2> Points = new List<Point2>();
        public Dictionary<Point2, int> PointsMap = new Dictionary<Point2, int>();
        public List<BasicHalfedge> Halfedges { get; private set; } = new List<BasicHalfedge>();

        public int NumExtraBevelingVertices { get; private set; }

        public List<BasicContour> BasicContours { get; private set; } = new List<BasicContour>();
        public List<BasicPolygon> Polygons { get; private set; } = new List<BasicPolygon>();

        public List<BasicTriangle> Triangles { get; private set; }

        public List<BasicVertex> Vertices { get; private set; } = new List<BasicVertex>();
        public List<BasicPolygonFace> PolygonFaces { get; private set; } = new List<BasicPolygonFace>();
        public DcelMesh Dcel { get; internal set; }

        public void _Inicializar(IMesh mesh, List<ContourVertex> contourVertices)
        {
            /*

                        var dcel = ToDCEL((Mesh)mesh);
                        Vertices = new List<BasicVertex>(dcel.Vertices.Count);
                        foreach (var v in dcel.Vertices)
                        {
                            assert(new Point2(v.X, v.Y) == contourVertices[v.ID].Point);
                            var bv = new BasicVertex(v.ID, new Point2(v.X, v.Y));
                            assert(bv.Id == Vertices.Count);
                            Vertices.Add(bv);
                        }
                        Triangles = new List<BasicTriangle>(dcel.Faces.Count);

                        foreach (var f in dcel.Faces)
                        {
                            var t = new BasicTriangle(f.ID);
                            assert(t.Id == Triangles.Count);
                            Triangles.Add(t);
                            var he = f.Edge;
                            for (var i = 0; i < 3; i++)
                            {
                                t.Vertices[i] = Vertices[he.Origin.ID];
                                he = he.Next;
                            }

                            Log($"area t:{t.Id} = {Mathex.TriangleArea(t.Vertices.Select(v => v.Point).ToArray())}");
                        }
                        foreach (var f in dcel.Faces)
                        {
                            var t = Triangles[f.ID];
                            var he = f.Edge;
                            for (var i = 0; i < 3; i++)
                            {
                                var nt = he.Next.Opposite.Face;
                                if (nt.ID >= 0)
                                {
                                    t.Neighbors[i] = Triangles[he.Next.Opposite.Face.ID];
                                }
                                else
                                {
                                    t.Neighbors[i] = null;
                                }
                                he = he.Next;
                            }
                        }
            */
        }

        public void DeleteBasicContour(BasicContour contour)
        {
            assert(contour != null, "Contour should not be null.");
            assert(BasicContours[contour.Id] == contour, "Contour should be in BasicContours.");

            BasicContours[contour.Id] = BasicContours[BasicContours.Count - 1];
            BasicContours[contour.Id].Id = contour.Id;
            BasicContours.RemoveAt(BasicContours.Count - 1);
        }
        public void DeleteBasicPolygon(BasicPolygon polygon)
        {
            assert(polygon != null, "Polygon should not be null.");
            assert(Polygons[polygon.Id] == polygon, "Polygon should be in Polygons.");
            Polygons[polygon.Id] = Polygons[Polygons.Count - 1];
            Polygons[polygon.Id].Id = polygon.Id;
            Polygons.RemoveAt(Polygons.Count - 1);
        }
        public void DeleteEdge(BasicHalfedge halfedge)
        {
            if (halfedge.Id > halfedge.Opposite.Id)
            {
                halfedge = halfedge.Opposite;
            }
            assert(halfedge.Id % 2 == 0, "Halfedge ID should be even.");
            assert(halfedge.Opposite.Id - halfedge.Id == 1, "Opposite halfedge ID should be one greater than the original halfedge ID.");
            assert(halfedge != null, "Halfedge should not be null.");
            assert(Halfedges[halfedge.Id] == halfedge, "Halfedge should be in Halfedges.");
            halfedge.Type = BasicHalfedgeType.Deleted;
            halfedge.Opposite.Type = BasicHalfedgeType.Deleted;

            Halfedges[halfedge.Id] = Halfedges[Halfedges.Count - 2];
            Halfedges[halfedge.Id + 1] = Halfedges[Halfedges.Count - 1];
            Halfedges[halfedge.Id].Id = halfedge.Id;
            Halfedges[halfedge.Id + 1].Id = halfedge.Id + 1;

            halfedge.Id = Halfedges.Count - 2;
            halfedge.Opposite.Id = Halfedges.Count - 1;

            Halfedges.RemoveAt(Halfedges.Count - 1);
            Halfedges.RemoveAt(Halfedges.Count - 1);
        }
        public void SwapHalfedges(int firstIdx, int seconIdx)
        {
            var he1 = Halfedges[firstIdx];
            Halfedges[firstIdx] = Halfedges[seconIdx];
            Halfedges[seconIdx] = he1;
            Halfedges[firstIdx].Id = firstIdx;
            Halfedges[seconIdx].Id = seconIdx;
        }

        public void CollapseEdgeContour(BasicHalfedge hc1, BasicHalfedge hc2)
        {
            BasicHalfedge ho1 = hc1.Opposite;
            BasicHalfedge ho2 = hc2.Opposite;

            assert(hc1.IsContour && !hc1.IsOutside, $"{hc1} no es Contorno .");
            assert(hc2.IsContour && !hc2.IsOutside, $"{hc2} no es Contorno .");
            assert(ho1.IsOutside && ho2.IsOutside, $"{ho1} and {ho2} no Son OutSide.");

            SwapHalfedges(hc2.Id, ho1.Id);
            hc1.LinkOpposite(hc2);
            ho1.LinkOpposite(ho2);
            Log($"  Van a Colapra ");
            Log($"    {ho1}");
            Log($"    {ho2}");
            Log($"  Quedan como E ");
            Log($"    {hc1}");
            Log($"    {hc2}");
            //Log($"");

            ho1.Prev.LinkNext(ho2.Next);
            ho2.Prev.LinkNext(ho1.Next);

            //if ((ho1.Next == ho2 && ho2.Next == ho1))
            //{
            //    assert(ho1.Face == ho2.Face, $"{ho1.Face?.Id} should be the same as {ho2.Face?.Id}.");
            //    var contour = (BasicContour)ho1.Face;
            //    DeleteBasicContour(contour);
            //    ho1.Face = ho2.Face = null;
            //}
            //else if (ho1.Next == ho2)
            //{
            //    assert(ho1.Face == ho2.Face, $"{ho1.Face?.Id} should be the same as {ho2.Face?.Id}.");
            //    if (ho1.Face.Halfedge == ho1 || ho1.Face.Halfedge == ho2)
            //    {
            //        ho1.Face.Halfedge = ho2.Next;
            //    }
            //    ho1.Prev.LinkNext(ho2.Next);
            //}
            //else if (ho2.Next == ho1)
            //{
            //    assert(ho1.Face == ho2.Face, $"{ho1.Face?.Id} should be the same as {ho2.Face?.Id}.");
            //    if (ho2.Face.Halfedge == ho1 || ho2.Face.Halfedge == ho2)
            //    {
            //        ho2.Face.Halfedge = ho1.Prev;
            //    }
            //    ho1.Prev.LinkNext(ho2.Next);
            //}
            //else if (ho1.Face == ho2.Face)
            //{
            //    // se separan los contornos
            //    var ho1Prev = ho1.Prev;
            //    var ho2Prev = ho2.Prev;
            //    var ho1Next = ho1.Next;
            //    var ho2Next = ho2.Next;

            //    ho1Prev.LinkNext(ho2Next);
            //    ho1Prev.Face = ho1.Face;
            //    ho1Prev.Face.Halfedge = ho1Prev;
            //    ho1Prev.Face.Assert();

            //    var newContour = NewContour();

            //    ho2Prev.LinkNext(ho1.Next);
            //    ho2Prev.Face = newContour;
            //    newContour.Halfedge = ho2Prev;

            //    assert(ho1.Face.Halfedge != null, $"{ho1.Face?.Id} should have a Halfedge.");
            //    assert(ho1.Face.Halfedge == ho2.Face.Halfedge, $"{ho1.Face?.Id} and {ho2.Face?.Id} should have the same Halfedge.");
            //    ho1.Prev.LinkNext(ho2.Next);
            //}
            //else
            //{
            //    throw new NotImplementedException($"CollapseEdgeContour not implemented for {ho1.Face?.Id} and {ho2.Face?.Id}.");
            //}
            DeleteEdge(ho1);
        }

        public BasicContour NewContour()
        {
            BasicContour outerFace = new BasicContour(BasicContours.Count);
            BasicContours.Add(outerFace);
            return outerFace;
        }
        public BasicPolygon NewPolygon(bool isHole = false)
        {
            BasicPolygon polygon = new BasicPolygon(Polygons.Count, isHole);
            Polygons.Add(polygon);
            return polygon;
        }
        public BasicPolygonFace NewPolygonFace(BasicPolygon polygon)
        {
            BasicPolygonFace polygonface = new BasicPolygonFace(PolygonFaces.Count, polygon);
            PolygonFaces.Add(polygonface);
            return polygonface;
        }

        public void CreateContourBasicHalfedges(IEnumerable<Point2> contourVertices, bool isHole = false)
        {
            BasicHalfedge first = null;
            BasicHalfedge last = null;

            BasicPolygon polygon = NewPolygon(isHole);

            List<BasicVertex> vertices = new List<BasicVertex>();
            foreach (Point2 p in contourVertices)
            {
                if (last != null && last.Vertex.Point.AreNear(p))
                {
                    Warning($"P{polygon.Id} Skip duplicate point {p}");
                    continue;
                }

                var vertex = NewVertex(p);
                var he = NewCountorEdges(vertex, polygon);

                if (last != null)
                {
                    last.FullLinkNext(he);
                }
                if (first == null)
                {
                    first = he;
                }
                last = he;
            }
            if (last != null)
            {
                last.FullLinkNext(first);
            }
            polygon.Halfedge = first;

            var i = new SurfNetPolyBool.Intersecter.RegionIntersecter(this, first);
            i.Calculate(false);

            foreach (var h in polygon.Halfedge.Circulation())
            {
                if (h.Vertex.Halfedge == null)
                    h.Vertex.Halfedge = h;
            }
            // var prueba = polygon.PolygonFace.Halfedge.Opposite.Circulation().Select(v => (v.Vertex.PointId)).ToArray();
        }

        public void Inicializar(IMesh mesh, Dictionary<(int start, int end), BasicHalfedge> constrains)
        {
            Dcel = SurfNet.Meshing.Converter.Instance.ToDCEL(this,(Mesh)mesh);

          //  ToDCEL((Mesh)mesh, constrains);

            //foreach (var v in dcel.Vertices)
            //{
            //    assert(new Point2(v.X, v.Y) == contourVertices[v.ID].Point);
            //    var bv = new BasicVertex(v.ID, new Point2(v.X, v.Y));
            //    assert(bv.Id == Vertices.Count);
            //    Vertices.Add(bv);

            //}
            Triangles = new List<BasicTriangle>(InnerFaces.Count);

            //         assert(InnerFaces.Count == Vertices.Count, "The number of inner faces should match the number of vertices.");
            var verticesMap = new Dictionary<(Point2, Point2), BasicHalfedge>();

            foreach (var fi in InnerFaces)
            {
                //  assert(fi.Id == fi.Halfedge.Vertex.Id, "The face ID should match the vertex ID of the halfedge.");
                (Point2, Point2) key = (fi.Halfedge.Vertex.Point, fi.Halfedge.Next.Vertex.Point);
                if (!verticesMap.ContainsKey(key))
                {
                    verticesMap[key] = fi.Halfedge;
                }
                else
                {
                    assert(verticesMap[key] == fi.Halfedge, $"The halfedge {fi.Halfedge} should match the existing {verticesMap[key]} for the vertex pair.");
                }
            }

            foreach (var t in Triangles)
            {
                //var t = new BasicTriangle(f.Id);
                //assert(t.Id == Triangles.Count);
                //Triangles.Add(t);
                //var he = f.Halfedge;
                //for (var i = 0; i < 3; i++)
                //{
                //    t.Vertices[i] = Vertices[he.Vertex.Id];
                //    he = he.Next;
                //}
                Log($"area t:{t.Id} = {Mathex.TriangleArea(t.Vertices.Select(v => v.Point).ToArray())}");
            }

            foreach (var f in Triangles)
            {
                var t = Triangles[f.Id];
                var he = f.Halfedge;
                for (var i = 0; i < 3; i++)
                {
                    var nt = he.Next.Opposite.Face;
                    if (nt.Id >= 0)
                    {
                        t.Neighbors[i] = Triangles[he.Next.Opposite.Face.Id];
                    }
                    else
                    {
                        t.Neighbors[i] = null;
                    }
                    he = he.Next;
                }
            }
        }

        public BasicHalfedge NewCountorEdges(BasicVertex vertex, BasicPolygon polygon)
        {
            var he = new BasicHalfedge(Halfedges.Count, BasicHalfedgeType.Countour | BasicHalfedgeType.Constraint);
            var other = new BasicHalfedge(Halfedges.Count + 1, BasicHalfedgeType.IsOutside | BasicHalfedgeType.Constraint);
            he.IsHole = polygon.IsHole;
            other.IsHole = polygon.IsHole;
            he.Vertex = vertex;
            Halfedges.Add(he);
            Halfedges.Add(other);
            he.LinkOpposite(other);
            return he;
        }

        public (BasicHalfedge, BasicHalfedge) NewPairHalfedges()
        {
            var he = new BasicHalfedge(Halfedges.Count, BasicHalfedgeType.Interior);
            var other = new BasicHalfedge(Halfedges.Count + 1, BasicHalfedgeType.Interior);
            Halfedges.Add(he);
            Halfedges.Add(other);
            return (he, other);
        }

        public BasicHalfedge NewEdges(BasicVertex? vertex)
        {
            var he = new BasicHalfedge(Halfedges.Count, BasicHalfedgeType.Interior);
            var other = new BasicHalfedge(Halfedges.Count + 1, BasicHalfedgeType.Interior);
            Halfedges.Add(he);
            Halfedges.Add(other);
            he.LinkOpposite(other);
            he.Vertex = vertex;
            return he;
        }

        public BasicVertex NewVertex(Point2 point)
        {
            if (!PointsMap.TryGetValue(point, out var pointId))
            {
                PointsMap[point] = pointId = Points.Count;
                Points.Add(point);
            }
            else
            {
            }

            var result = new BasicVertex(Vertices.Count, point, pointId);
            Vertices.Add(result);
            return result;
        }

        public BasicVertex NewVertex(BasicVertex vertex)
        {
            var result = new BasicVertex(Vertices.Count, vertex.Point, vertex.OriginalId, vertex.PointId);
            Vertices.Add(result);
            return result;
        }

        /// <summary>
        /// Convert the triangle mesh topology to DCEL.
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        //public void ToDCEL(Mesh mesh, Dictionary<(int start, int end), BasicHalfedge> constrainsMap)
        //{
        //    Triangles = new List<BasicTriangle>();

        //    if (mesh == null)
        //    {
        //        throw new ArgumentNullException(nameof(mesh), "Mesh cannot be null.");
        //    }
        //    mesh.Renumber();
        //    var pointIdMap = new int[mesh.vertices.Count];
        //    var vertexIdMap = new int[mesh.vertices.Count];

        //    Array.Fill(pointIdMap, -1);
        //    foreach (var kv in mesh.vertices)
        //    {
        //        var point = new Point2(kv.Value.X, kv.Value.Y);
        //        if (this.PointsMap.TryGetValue(point, out var pointId))
        //        {
        //            pointIdMap[kv.Key] = pointId;
        //            vertexIdMap[pointId] = kv.Key;
        //        }
        //        else
        //        {
        //            Warning($" no se encontro punto {point} id {kv.Key} ");
        //        }
        //    }

        //    // Maps a triangle to its 3 edges (used to set next pointers).
        //    var map = new List<BasicHalfedge>[mesh.triangles.Count];
        //    var triangles = new BasicTriangle[mesh.triangles.Count];
        //    //// Maps a triangle to its 3 edges (used to set next pointers).
        //    //var map = new List<BasicHalfedge>[mesh.triangles.Count];

        //    foreach (var t in mesh.triangles)
        //    {
        //        var bt = new BasicTriangle(t.id);
        //        triangles[t.id] = bt;
        //        map[t.id] = new List<BasicHalfedge>(3);
        //    }

        //    BasicTriangle face;

        //    foreach (var t in mesh.triangles)
        //    {
        //        face = triangles[t.ID];
        //        face.Id = t.id;

        //        //   faces[t.id] = face;

        //        map[t.id] = new List<BasicHalfedge>(3);
        //    }

        //    Otri tri = default(Otri), neighbor = default(Otri);
        //    TVertex org, dest;

        //    int id, nid, count = mesh.triangles.Count;

        //    BasicHalfedge edge, twin, next;

        //    // Count half-edges (edge ids).
        //    int k = 0;

        //    // Maps a vertex to its leaving boundary edge.
        //    var boundary = new List<BasicHalfedge>();
        //    int constraintsFounded = 0;
        //    foreach (var t in mesh.triangles)
        //    {
        //        id = t.id;

        //        tri.tri = t;

        //        for (int i = 0; i < 3; i++)
        //        {
        //            tri.orient = i;
        //            tri.Sym(ref neighbor);

        //            nid = neighbor.tri.id;

        //            if (id < nid || nid < 0)
        //            {
        //                face = triangles[id];

        //                // Get the endpoints of the current triangle edge.
        //                org = tri.Org();
        //                dest = tri.Dest();

        //                if (org.id > dest.id)
        //                {
        //                    // Swap the endpoints to maintain a consistent ordering.
        //                    (org, dest) = (dest, org);
        //                }
        //                var poId = pointIdMap[org.id];

        //                assert(Points[poId].X == mesh.vertices[org.id].X && Points[poId].Y == mesh.vertices[org.id].Y);

        //                var pdId = pointIdMap[dest.id];
        //                assert(Points[pdId].X == mesh.vertices[dest.id].X && Points[pdId].Y == mesh.vertices[dest.id].Y);
        //                bool found = false;
        //                if (constrainsMap.TryGetValue((poId, pdId), out var value))
        //                {
        //                    found = true;
        //                }
        //                else if (constrainsMap.TryGetValue((pdId, poId), out value))
        //                {
        //                    found = true;
        //                }
        //                if (found)
        //                {
        //                    edge = value;
        //                    twin = value.Opposite;
        //                    if (twin.IsOutside)
        //                    {
        //                        boundary.Add(edge);
        //                    }
        //                    constraintsFounded++;
        //                    Log($"fount {edge}");
        //                }
        //                else
        //                {
        //                    // Create half-edges.
        //                    edge = NewEdges(null);
        //                    twin = edge.Opposite;
        //                }
        //                edge.Face = face;
        //                edge.TagPoint = poId;
        //                if (nid >= 0)
        //                    twin.Face = triangles[nid];
        //                else
        //                {
        //                    Warning($"No se encontro el triangulo {nid} en el mapa de triangulos");
        //                }
        //                twin.TagPoint = pdId;

        //                map[id].Add(edge);

        //                if (nid >= 0)
        //                {
        //                    map[nid].Add(twin);
        //                }

        //                // Set leaving edges.
        //                //edge.origin.leaving = edge;
        //                //twin.origin.leaving = twin;

        //                // Set twin edges.
        //                //edge.twin = twin;
        //                //twin.twin = edge;

        //                //edge.id = k++;
        //                //twin.id = k++;

        //                //edges.Add(edge);
        //                //edges.Add(twin);
        //            }
        //        }
        //    }

        //    // Set next pointers for each triangle face.
        //    foreach (var t in map)
        //    {
        //        edge = t[0];
        //        next = t[1];
        //        if (edge.Face.Halfedge == null) edge.Face.Halfedge = edge;
        //        var eoi = edge.Opposite.Vertex == null ? edge.Opposite.TagPoint : edge.Opposite.Vertex.PointId;
        //        var nexti = next.Vertex == null ? next.TagPoint : next.Vertex.PointId;

        //        if (eoi == nexti)
        //        {
        //            edge.LinkNext(next);
        //            next.LinkNext(t[2]);
        //            t[2].LinkNext(edge);
        //        }
        //        else
        //        {
        //            edge.LinkNext(t[2]);
        //            next.LinkNext(edge);
        //            t[2].LinkNext(next);
        //        }

        //        Log($" ");
        //    }
        //    foreach (var t in triangles)
        //    {
        //        t.Halfedge.Circulation();
        //        assert(t.Halfedge.Circulation().Count() == 3, $" {t} no tiene tres lados");
        //        assert(Mathex.TriangleArea(t.Halfedge.Circulation().Select(h => Points[h.TagPoint]).ToArray()) > 0, $" {t} tiene area negativa");
        //    }

        //    foreach (var bv in Vertices)
        //    {
        //        var cv = bv;
        //        foreach (var he in bv.Halfedge.AroundVerticeCW().Skip(1))
        //        {
        //            if (he.IsConstraint)
        //            {
        //                if (he.IsContour)
        //                {
        //                    //       cv = he.Vertex;
        //                }
        //                else
        //                {
        //                    break;
        //                }
        //            }
        //            if (he.Vertex == null)
        //                he.Vertex = bv;
        //        }
        //    }
        //    foreach (var t in triangles)
        //    {
        //        Log($" aca estoy");
        //        Log($" {t.Id} Area:{Mathex.TriangleArea(t.Halfedge.Circulation().Select(h => Points[h.TagPoint]).ToArray())}");
        //        t.Halfedge.Circulation();
        //        assert(t.Halfedge.Circulation().Count() == 3, $" {t} no tiene tres lados");
        //        assert(Mathex.TriangleArea(t.Halfedge.Circulation().Select(h => Points[h.TagPoint]).ToArray()) > 0, $" {t} tiene area negativa");
        //    }

        //    Triangles.AddRange(triangles);
        //}

        //public void Triangulate(IMesh mesh, List<ContourVertex> contourVertices)
        //{
        //    var Dcel = Converter.Instance.ToDCEL((TriangleNet.Mesh)mesh);

        //    IEnumerable<(int a, int b)[]> EnumerateTriangules(IMesh mesh)
        //    {
        //        Otri tri = default;
        //        Otri neighbor = default;
        //        Osub sub = default;

        //        TriangleNet.Geometry.Vertex p1, p2;

        //        int ii = 0;
        //        foreach (var t in mesh.Triangles)
        //        {
        //            ii++;
        //            tri.tri = t;
        //            tri.orient = 0;

        //            (int a, int b)[] segments = new (int a, int b)[3];
        //            for (int i = 0; i < 3; i++)
        //            {
        //                tri.Sym(ref neighbor);

        //                int nid = neighbor.tri.id;

        //                if (true/* (nid != -1)*/)
        //                {
        //                    var backColor = Console.ForegroundColor;
        //                    if ((nid == -1))
        //                    {
        //                        Console.ForegroundColor = ConsoleColor.Blue;
        //                    }
        //                    else if ((tri.tri.id >= nid))
        //                    {
        //                        Console.ForegroundColor = ConsoleColor.Yellow;
        //                    }
        //                    p1 = tri.Org();
        //                    p2 = tri.Dest();

        //                    tri.Pivot(ref sub);

        //                    segments[i] = (p1.ID, p2.ID);

        //                    Console.ForegroundColor = backColor;
        //                }

        //                tri.orient++;
        //            }

        //            yield return segments;
        //        }
        //    }

        //    Triangles = new List<BasicTriangle>(mesh.Triangles.Count);
        //    Vertices = new List<BasicVertex>(contourVertices.Count);

        //    foreach (var v in contourVertices)
        //    {
        //        BasicVertex bv = new BasicVertex(Vertices.Count, v.Point);
        //        bv.degree = v.AngleType == EAngle.Reflex ? 1 : 0;
        //        assert(bv.Id == v.ID);
        //        if (bv.degree == 1)
        //        {
        //            NumExtraBevelingVertices++;
        //        }
        //        Vertices.Add(bv);
        //    }

        //    var tt = mesh.Triangles.Where(t => t.GetHashCode() >= 0);

        //    void assign_halfedge(BasicTriangle t, int a, int b)
        //    {
        //        bool oposite = false;
        //        if (a > b)
        //        {
        //            oposite = true;

        //            (a, b) = (b, a);
        //        }

        //        TrianglesOnEdge value;
        //        if (!halfedges.TryGetValue((a, b), out value))
        //        {
        //            value = new TrianglesOnEdge();
        //            halfedges.Add((a, b), value);
        //        }

        //        if (!oposite)
        //        {
        //            assert(value.T == null);
        //            value.T = t;
        //        }
        //        else
        //        {
        //            assert(value.N == null);
        //            value.N = t;
        //        }
        //    }

        //    void assign_neighbors(BasicTriangle t, int idx, int a, int b, BasicTriangle?[] neighbors)
        //    {
        //        assert(a < Vertices.Count);

        //        var cva = contourVertices[a];
        //        var cvb = contourVertices[b];

        //        if (cva.NextInLAV.ID != b)
        //        {
        //            bool oposite = false;

        //            if (a > b)
        //            {
        //                oposite = true;

        //                (a, b) = (b, a);
        //            }

        //            TrianglesOnEdge value;
        //            if (!halfedges.TryGetValue((a, b), out value))
        //            {
        //                throw new Exception("Esto no deberia pasar");
        //            }

        //            if (oposite)
        //            {
        //                assert(value.T != null);
        //                assert(value.T != t);
        //                assert(value.N == t);
        //                neighbors[idx] = value.T;
        //            }
        //            else
        //            {
        //                assert(value.N != null);
        //                assert(value.N != t);
        //                assert(value.T == t);
        //                neighbors[idx] = value.N;
        //            }
        //        }
        //    }

        //    foreach (var tri in EnumerateTriangules(mesh))
        //    {
        //        var t = new BasicTriangle(Triangles.Count);
        //        assert(t.Id == Triangles.Count);
        //        Triangles.Add(t);
        //        var a = tri[0].a;
        //        var b = tri[1].a;
        //        var c = tri[2].a;

        //        //var wavefrontedges = new WavefrontEdge?[3];
        //        t.Vertices[0] = Vertices[a];
        //        t.Vertices[1] = Vertices[b];
        //        t.Vertices[2] = Vertices[c];

        //        assign_halfedge(t, a, b);
        //        assign_halfedge(t, b, c);
        //        assign_halfedge(t, c, a);

        //        //t.set_wavefronts(wavefrontedges[0], wavefrontedges[1], wavefrontedges[2]);
        //        //foreach (var we in wavefrontedges)
        //        //{
        //        //    if (we != null)
        //        //    {
        //        //        we.set_initial_incident_triangle(t);
        //        //    }
        //        //}
        //    }

        //    foreach (var t in Triangles)
        //    {
        //        var neighbors = t.Neighbors.Array;
        //        assign_neighbors(t, 0, t.Vertices[1].Id, t.Vertices[2].Id, neighbors);
        //        assign_neighbors(t, 1, t.Vertices[2].Id, t.Vertices[0].Id, neighbors);
        //        assign_neighbors(t, 2, t.Vertices[0].Id, t.Vertices[1].Id, neighbors);
        //    }
        //}

        internal int get_num_extra_beveling_vertices() => NumExtraBevelingVertices;

        internal int get_total_degree()
        {
            return halfedges.Count * 2;
        }

        //internal BasicEdge get_edge(int vidxu, int vidxv)
        //{
        //    throw new NotImplementedException();
        //}
        internal bool has_edge(int a, int b)
        {
            if (a > b)
            {
                (a, b) = (b, a);
            }
            return halfedges.ContainsKey((a, b));
        }

        internal void JoinOpposite(BasicHalfedge h1, BasicHalfedge h2)
        {
            var oh1 = h1.Opposite;
            var oh2 = h2.Opposite;
            //    if (h1.Opposite.Face is BasicContour && h2.Opposite.Face is BasicContour)
            //    {
            //        Join two contour halfedges
            //        h1.LinkOpposite(h2);

            //    }
            //    else
            //    {
            //        Join two internal halfedges
            //       assert(h1.Opposite.Face is BasicPolygon && h2.Opposite.Face is BasicPolygon, "Joining non-contour halfedges should be between polygons.");
            //h1.LinkOpposite(h2);
            //        h2.LinkOpposite(h1);
            //    }
        }

        private class TrianglesOnEdge
        {
            public int a;

            public int b;

            public BasicTriangle? N;

            public BasicTriangle? T;

            public TrianglesOnEdge()
            {
            }
        }

        //internal List<BasicEdge> edges()
        //{
        //    throw new NotImplementedException();
        //}
    }
}