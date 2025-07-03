using TriangleNet;
using TriangleNet.Geometry;
using TriangleNet.Meshing;
using TriangleNet.Topology;


namespace SurfNet
{
    using static DebugLog;
    using static Mathex;
    public partial class KineticTriangulation
    {





        internal int HeIDs = 0;
       

    

        public bool IsIntialRefineTriangleDoit { get; set; }=false;



        public KineticTriangulation EnterContour(IEnumerable<SurfNet.Point2> aPoly, bool aCheckValidity = true)
        {
            if (aCheckValidity)
            {
                List<SurfNet.Point2> nPoly = new List<SurfNet.Point2>();
                foreach (var p in aPoly)
                {
                    if (nPoly.Count == 0 || !Mathex.are_near(nPoly.Last(), p, 0.05))
                    {
                        nPoly.Add(p);
                    }
                }
                if (nPoly.Count > 0 && Mathex.are_near(nPoly.First(), nPoly.Last(), 0.05))
                {
                    nPoly.RemoveAt(nPoly.Count - 1);
                }
                if (nPoly.Count < 3)
                    throw new Exception("Degenerate contour (less than 3 non-degenerate vertices).");
                aPoly = nPoly;
            };
            EnterValidContour(aPoly);

            return this;
        }

        private void EnterValidContour(IEnumerable<SurfNet.Point2> points)
        {
            var polygon = new Polygon(Polygons.Count + 1, ContourVertices.Count, points.Count(), Polygons.Count != 0);

            ContourVertices.Capacity = ContourVertices.Count + polygon.Count;
            ContourVertices.AddRange(points.Select((p, i) => new ContourVertex(i + polygon.StartIndex, p)));

            DebugLog.Log("Inserting Connected Component of the Boundary....");

            //Halfedge lFirstCCWBorder = Halfedge.NULL;
            //Halfedge lPrevCCWBorder = Halfedge.NULL;
            //Halfedge lNextCWBorder = Halfedge.NULL;
            ContourVertex lFirstVertex = ContourVertex.NULL;
            ContourVertex lPrevVertex = ContourVertex.NULL;

            // InputPointIterator lPrev = aBegin ;

            int c = 0;

            var last = ContourVertices.Last();
            foreach (var curr in ContourVertices.Skip(polygon.StartIndex))
            {
                this.TopLeft = Mathex.Min(this.TopLeft, curr.Point);

                this.BottomRight = Mathex.Max(this.BottomRight, curr.Point);


                var lVertex = curr;

                if (c == 0)
                {
                    lFirstVertex = lVertex;

                }
                else
                {
                    lVertex.PrevInLAV = lPrevVertex;
                    lPrevVertex.NextInLAV = lVertex;

                    // SetVertexTriedge(lPrevVertex, new Triedge(lPrevCCWBorder, lCCWBorder));



                }
                ++c;

                // lPrev = lCurr ;

                lPrevVertex = lVertex;

            }
            this.Polygons.Add(polygon);

            lFirstVertex.PrevInLAV = lPrevVertex;
            lPrevVertex.NextInLAV = lFirstVertex;


            foreach (var curr in ContourVertices.Skip(polygon.StartIndex))
            {
                curr.AngleType= GetAngleType( curr.PrevInLAV.Point,curr.Point,curr.NextInLAV.Point);
            }
            



            //  SetVertexTriedge(lPrevVertex, new Triedge(lPrevCCWBorder, lFirstCCWBorder));



        }


        public Mesh Triangulate()
        {
            var polygons = new TriangleNet.Geometry.Polygon();

            foreach (var polygon in Polygons.Select(e => e))
            {
                TriangleNet.Geometry.Vertex toVector(SurfNet.Point2 v, int i) => new TriangleNet.Geometry.Vertex(v.X, v.Y, i);
                Contour pts = new Contour(polygon.Vertices(this.ContourVertices)
                                                .Select((p, i) => toVector(p.Point, i + polygon.StartIndex)), polygon.ID);
                bool hole = polygon.Hole == true;
                polygons.Add(pts, hole);
            }

            var options = new ConstraintOptions() { Convex = false };

            var quality = new QualityOptions() { };

            // Generate mesh using the polygons Triangulate extension method.

            Mesh mesh = (Mesh)polygons.Triangulate(options, quality);

            return mesh;
        }

        private class TrianglesOnEdge
        {
            public TrianglesOnEdge()
            {

            }
            public int a;
            public int b;

            public KineticTriangle? T;
            public KineticTriangle? N;

        }
        
        public void BuildKineticTriangles(IMesh mesh)
        {



            IEnumerable<(int a, int b)[]> EnumerateTriangules(IMesh mesh)
            {
                Otri tri = default;
                Otri neighbor = default;
                Osub sub = default;

                TriangleNet.Geometry.Vertex p1, p2;

                int ii = 0;
                foreach (var t in mesh.Triangles)
                {
                    ii++;
                    tri.tri = t;
                    tri.orient = 0;
                  
                    (int a, int b)[] segments = new (int a, int b)[3];
                    for (int i = 0; i < 3; i++)
                    {
                        tri.Sym(ref neighbor);

                        int nid = neighbor.tri.id;

                        if (true/* (nid != -1)*/)
                        {
                            var backColor = Console.ForegroundColor;
                            if ((nid == -1))
                            {
                                Console.ForegroundColor = ConsoleColor.Blue;
                            }
                            else if ((tri.tri.id >= nid))
                            {
                                Console.ForegroundColor = ConsoleColor.Yellow;
                            }
                            p1 = tri.Org();
                            p2 = tri.Dest();

                            tri.Pivot(ref sub);


                            segments[i] = (p1.ID, p2.ID);



                            Console.ForegroundColor = backColor;
                        }

                        tri.orient++;
                    }
                    Console.WriteLine();
                    yield return segments;
                }
            }



       

            this.triangles = new List<KineticTriangle>(mesh.Triangles.Count);


            Dictionary<(int a, int b), TrianglesOnEdge> kineticHalfedges = new Dictionary<(int a, int b), TrianglesOnEdge>();

            var tt = mesh.Triangles.Where(t => t.GetHashCode() >= 0);



            bool assign_wavefrontedge(KineticTriangle t, int idx, int a, int b, WavefrontEdge?[] wavefrontEdges)
            {
                assert(a < vertices.Count);

                var wavefront_edge = this.wavefront_edges[a];
                assert(wavefront_edge != null);

                assert(wavefront_edge.initial_vertex(0).Id == a);
                if (wavefront_edge.initial_vertex(1).Id == b)
                {
                    wavefrontEdges[TriangulationUtils.cw(idx)] = wavefront_edge;

                    return true;
                }
                else
                {
                    bool oposite = false;
                    if (a > b)
                    {
                        oposite = true;

                        (a, b) = (b, a);
                    }

                    TrianglesOnEdge value;
                    if (!kineticHalfedges.TryGetValue((a, b), out value))
                    {
                        value = new TrianglesOnEdge();
                        kineticHalfedges.Add((a, b), value);
                    }

                    if (!oposite)
                    {
                        assert(value.T == null);
                        value.T = t;
                    }
                    else
                    {
                        assert(value.N == null);
                        value.N = t;
                    }
                    return false;
                }




            }

            void assign_neighbors(KineticTriangle t, int idx, int a, int b, KineticTriangle?[] neighbors)
            {
                assert(a < vertices.Count);

                if (t.wavefronts[idx] == null)
                {
                    bool oposite = false;

                    if (a > b)
                    {
                        oposite = true;

                        (a, b) = (b, a);
                    }

                    TrianglesOnEdge value;
                    if (!kineticHalfedges.TryGetValue((a, b), out value))
                    {
                        throw new Exception("Esto no deberia pasar");
                    }

                    if (oposite)
                    {

                        assert(value.T != null);
                        assert(value.T != t);
                        assert(value.N == t);
                        neighbors[idx] = value.T;
                    }
                    else
                    {
                        assert(value.N != null);
                        assert(value.N != t);
                        assert(value.T == t);
                        neighbors[idx] = value.N;
                    }

                }




            }


            foreach (var tri in EnumerateTriangules(mesh))
            {
                var t = new KineticTriangle(triangles.Count, 0, null);
                var a = tri[0].a;
                var b = tri[1].a;
                var c = tri[2].a;

                var wavefrontedges = new WavefrontEdge?[3];
                t.set_vertex(0, vertices[a]);
                t.set_vertex(1, vertices[b]);
                t.set_vertex(2, vertices[c]);

                assign_wavefrontedge(t, 0, a, b, wavefrontedges);
                assign_wavefrontedge(t, 1, b, c, wavefrontedges);
                assign_wavefrontedge(t, 2, c, a, wavefrontedges);

                t.set_wavefronts(wavefrontedges[0], wavefrontedges[1], wavefrontedges[2]);
                foreach (var we in wavefrontedges)
                {
                    if (we != null)
                    {
                        we.set_initial_incident_triangle(t);
                    }
                }

                triangles.Add(t);

            }

            foreach (var t in triangles)
            {

                var neighbors = new KineticTriangle?[3];
                assign_neighbors(t, 0, t.vertex(1).Id, t.vertex(2).Id, neighbors);
                assign_neighbors(t, 1, t.vertex(2).Id, t.vertex(0).Id, neighbors);
                assign_neighbors(t, 2, t.vertex(0).Id, t.vertex(1).Id, neighbors);
                t.set_neighbors(neighbors);
            }

        }

        public List<ContourVertex> ContourVertices { get; set; } = new List<ContourVertex>();

        public List<Polygon> Polygons { get; set; } = new List<Polygon> { };
        public Point2 TopLeft { get; private set; }
        public Point2 BottomRight { get; private set; }


      public void Trangulate(BasicInput input)
        {
            input.Triangulate(Triangulate(), ContourVertices);
        }
        public void Initialize()
        {
            foreach (var cv in ContourVertices)
            {
            
                var wavefront_edge = new WavefrontEdge(cv.Point, cv.NextInLAV.Point, 1, null, null);
                //assert(wavefront_edge.id == cv.ID);
                wavefront_edges.Add(wavefront_edge);
            }
            foreach (var cv in ContourVertices)
            {
                var left = wavefront_edges[cv.PrevInLAV.ID];
                var right = wavefront_edges[cv.ID];


                vertices.make_initial_vertex(cv.Point, left, right, true);
                assert(vertices.Last().Id == cv.ID);

            }
            foreach (var cv in ContourVertices)
            {
               

                var wavefront_edge = this.wavefront_edges[cv.ID];

                wavefront_edge.set_wavefrontedge_vertex(0, this.vertices[cv.ID]);
                wavefront_edge.set_wavefrontedge_vertex(1, this.vertices[cv.NextInLAV.ID]);
                wavefront_edge.set_initial_vertices();
            }
            var mesh = Triangulate();
            BuildKineticTriangles(mesh);
           // BasicInput.Triangulate(mesh,ContourVertices );

            //foreach(var t in triangles)
            //{
            //    var bt = BasicInput.Triangles[t.Id];
            //    assert(bt != null);
            //    for(int i = 0; i < 3; i++)
            //    {
            //        assert((bt.Vertices[i] == null && t.vertex(i) == null) || bt.Vertices[i]?.Id == t.vertex(i)?.Id);
            //        assert((bt.Neighbors[i] ==null && t.neighbor(i)==null)||  bt.Neighbors[i]?.Id == t.neighbor(i)?.Id );

            //    }
            //}


            tidx_in_check_refinement = new bool[triangles.Count];
            assert_valid(0, CORE_ZERO);
            initialized = true;
        }



        //public static DcelMesh ToDCEL(Mesh mesh)
        //{
        //    var dcel = new DcelMesh();

        //    var vertices = new HVertex[mesh.vertices.Count];
        //    var faces = new DcelFace[mesh.triangles.Count];

        //    dcel.HalfEdges.Capacity = 2 * mesh.NumberOfEdges;

        //    mesh.Renumber();

        //    HVertex vertex;

        //    foreach (var v in mesh.vertices.Values)
        //    {
        //        vertex = new HVertex(v.x, v.y);
        //        vertex.ID = v.id;
        //        vertex.Label = v.label;

        //        vertices[v.id] = vertex;
        //    }

        //    // Maps a triangle to its 3 edges (used to set next pointers).
        //    var map = new List<DcelHalfEdge>[mesh.triangles.Count];

        //    DcelFace face;

        //    foreach (var t in mesh.triangles)
        //    {
        //        face = new DcelFace(null);
        //        face.id = t.id;

        //        faces[t.id] = face;

        //        map[t.id] = new List<DcelHalfEdge>(3);
        //    }

        //    Otri tri = default(Otri), neighbor = default(Otri);
        //    TVertex org, dest;

        //    int id, nid, count = mesh.triangles.Count;

        //    DcelHalfEdge edge, twin, next;

        //    var edges = dcel.HalfEdges;

        //    // Count half-edges (edge ids).
        //    int k = 0;

        //    // Maps a vertex to its leaving boundary edge.
        //    var boundary = new Dictionary<int,DcelHalfEdge>();

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
        //                face = faces[id];

        //                // Get the endpoints of the current triangle edge.
        //                org = tri.Org();
        //                dest = tri.Dest();

        //                // Create half-edges.

        //                edge = new DcelHalfEdge(vertices[org.id], face);


        //                twin = new DcelHalfEdge(vertices[dest.id], nid < 0 ? DcelFace.Empty : faces[nid]);

        //                map[id].Add(edge);

        //                if (nid >= 0)
        //                {
        //                    map[nid].Add(twin);
        //                }
        //                else
        //                {
        //                    boundary.Add(dest.id, twin);
        //                }

        //                // Set leaving edges.
        //                edge.origin.Leaving = edge;
        //                twin.origin.Leaving = twin;

        //                // Set twin edges.
        //                edge.twin = twin;
        //                twin.twin = edge;

        //                edge.id = k++;
        //                twin.id = k++;

        //                edges.Add(edge);
        //                edges.Add(twin);
        //            }
        //        }
        //    }

        //    // Set next pointers for each triangle face.
        //    foreach (var t in map)
        //    {
        //        edge = t[0];
        //        next = t[1];

        //        if (edge.twin.origin.ID == next.origin.ID)
        //        {
        //            edge.next = next;
        //            next.next = t[2];
        //            t[2].next = edge;
        //        }
        //        else
        //        {
        //            edge.next = t[2];
        //            next.next = edge;
        //            t[2].next = next;
        //        }
        //    }

        //    // Resolve boundary edges.
        //    foreach (var e in boundary.Values)
        //    {
        //        e.next = boundary[e.twin.origin.ID];
        //    }

        //    dcel.Vertices.AddRange(vertices);
        //    dcel.Faces.AddRange(faces);

        //    return dcel;
        //}


        //public List<WavefrontVertex> WavefrontVertices = new List<WavefrontVertex>();
        //public List<KineticTriangle> KineticTriangles = new List<KineticTriangle>();
        //public List<DcelHalfEdge> Halfedges = new List<DcelHalfEdge>();

        //public void Initialize(TriangleNet.Mesh mesh)
        //{

        //    var polysIdSeek = Polygons.Select(p => p.StartIndex).ToArray();
        //    //~





        //    WavefrontVertices.Capacity = ContourVertices.Count();
        //    WavefrontVertices.AddRange(ContourVertices.Select((cv, idx) => new WavefrontVertex(idx, new Point2(cv.Point.X, cv.Point.Y))));


        //    foreach (var kv in WavefrontVertices)
        //    {

        //        kv.NextInLAV = WavefrontVertices[ContourVertices[kv.ID].NextInLAV.ID];
        //        kv.PrevInLAV = WavefrontVertices[ContourVertices[kv.ID].PrevInLAV.ID];
        //    }




        //    this.KineticTriangles.Capacity = mesh.triangles.Count;
        //    this.Halfedges.Capacity = 2 * mesh.NumberOfEdges;

        //    for (int i = 0; i < this.KineticTriangles.Capacity; i++)
        //    {
        //        this.KineticTriangles.Add(null);
        //    }
        //    for (int i = 0; i < this.Halfedges.Capacity; i++)
        //    {
        //        this.Halfedges.Add(null);
        //    }


        //    List<ContourVertex> contourVertices = new List<ContourVertex>();

        //    var vertices = this.WavefrontVertices;
        //    var faces = this.KineticTriangles;
        //    var edges = this.Halfedges;

        //    mesh.Renumber();

        //    WavefrontVertex vertex;

        //    foreach (var v in mesh.vertices.Values)
        //    {
        //        vertex = new WavefrontVertex(v.ID, new Point2(v.X, v.Y));
        //    }






        //}





        private void Create_supporting_lines()
        {
            ////DBG_FUNC_BEGIN(//DBG_KT_SETUP);
            ///* total number of wavefront edges during the entire propragation.
            // * this is an upper bound.
            // *
            // * num_t is a really rough upper bound on the number of split events.
            // */
            ////int num_t = num_initial_triangles + input.get_num_extra_beveling_vertices();
            ////var num_wavefront_edges = input.edges().size() * 2 + input.get_num_extra_beveling_vertices() + num_t * 2;
            ////wavefront_edges.Capacity = (num_wavefront_edges);

            //foreach (var fit in  triangles)
            //{
            // //   if (!fit.info().matches_component(restrict_component_)) continue;
            //    WavefrontEdge3 w = new WavefrontEdge3(  null, null, null );
            //    KineticTriangle t = fit;


            //    foreach ( var (e,i) in t.Edges.Select( (e,i)=>(e,i)))
            //    {
            //        if ( e.Twin==null)
            //        {
            //            // TODO Validate
            //            var u = t.Edges[TriangulationUtils.ccw(i)].Origin.DcelPoint;
            //            var v = t.Edges[TriangulationUtils.cw(i)].Origin.DcelPoint;

            //            //BasicEdge edge = input.get_edge(vidxu, vidxv);

            //            var wfe = new WavefrontEdge(u, v, 0.0/* edge.weight*/, t, /*T*/ skeleton_dcel_face);
            //            wavefront_edges.Add(wfe);
            //            w[i] = wfe;

            //            //int vidxu = get_basic_vertex_idx_from_triangle_vertex_indices(input, triangle_original_vertex_indices, triangle_idx, ccw(i));
            //            //int vidxv = get_basic_vertex_idx_from_triangle_vertex_indices(input, triangle_original_vertex_indices, triangle_idx, cw(i));
            //            //BasicVertex u = input.vertices()[vidxu];
            //            //BasicVertex v = input.vertices()[vidxv];



            //            //if (!input.has_edge(vidxu, vidxv))
            //            //{
            //            //    //LOG(ERROR) << "Cannot find edge for triangulation constraint (" << vidxu << ", " << vidxv << ") in input.  Input is not a PSLG.";
            //            //    //LOG(ERROR) << " v" << vidxu << ": " << CGAL_point(input.vertices()[vidxu].p);
            //            //    //LOG(ERROR) << " v" << vidxv << ": " << CGAL_point(input.vertices()[vidxv].p);
            //            //    //LOG(ERROR) << " Probably one of these points is on the interior of another edge.";
            //            //    throw new Exception();
            //            //}
            //            //BasicEdge edge = input.get_edge(vidxu, vidxv);

            //            ////DBG(//DBG_KT_SETUP) << "at " << t << "/" << i << ": vertices are " << vidxu << " and " << vidxv;

            //            // WavefrontEdge buddy_wavefront_edge = null;

            //            //var ct_n_hdl = fit.neighbors[i];
            //            //if (fit.neighbors[i]==null)
            //            //{
            //            //    /* neighbor not in map */
            //            //    //assert(restrict_component_ >= 0);
            //            //    //assert(!ct_n_hdl.info().matches_component(restrict_component_));
            //            //}
            //            //else
            //            //{
            //            //    KineticTriangle  n = fit.neighbors[i];
            //            //   // assert(ct_n_hdl.info().matches_component(restrict_component_));
            //            //    int idx_in_n = ct_n_hdl.index(fit);
            //            //    buddy_wavefront_edge = n.wavefronts[idx_in_n]; /* may still be null, and that's OK */
            //            //}

            //            //SkeletonDCELFace skeleton_dcel_face = skeleton.setup_new_input_edge(buddy_wavefront_edge /* may be null and that's OK */);


            //        }
            //        else
            //        {
            //           t.wavefronts[i] = null;
            //        }
            //    }

            //    t.set_wavefronts(w[0], w[1], w[2]);
            //}
            ////DBG_FUNC_END(//DBG_KT_SETUP);
        }





    }

}
