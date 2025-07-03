
using System.Security.Cryptography.X509Certificates;
using TriangleNet.Meshing;
using TriangleNet.Topology;

namespace SurfNet
{
    using static DebugLog;
    public class BasicInput
    {

        private class TrianglesOnEdge
        {
            public TrianglesOnEdge()
            {

            }
            public int a;
            public int b;

            public BasicTriangle? T;
            public BasicTriangle? N;

        }
        Dictionary<(int a, int b), TrianglesOnEdge> halfedges = new Dictionary<(int a, int b), TrianglesOnEdge>();
        public int NumExtraBevelingVertices { get; private set; }
        public List<BasicTriangle> Triangles { get; private set; }
        public List<BasicVertex> Vertices { get; private set; }

        //internal List<BasicEdge> edges()
        //{
        //    throw new NotImplementedException();
        //}

        //internal BasicEdge get_edge(int vidxu, int vidxv)
        //{
        //    throw new NotImplementedException();
        //}


        internal int get_num_extra_beveling_vertices() => NumExtraBevelingVertices;


        internal bool has_edge(int a, int b)
        {
            if (a > b)
            {

                (a, b) = (b, a);
            }
            return halfedges.ContainsKey((a, b));

        }

        internal List<BasicVertex> vertices()
        {
            throw new NotImplementedException();
        }


        public void Triangulate(IMesh mesh, List<ContourVertex> contourVertices)
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




            Triangles = new List<BasicTriangle>(mesh.Triangles.Count);
            Vertices = new List<BasicVertex>(contourVertices.Count);

            foreach (var v in contourVertices)
            {
                BasicVertex bv = new BasicVertex(Vertices.Count, v.Point, v.AngleType == EAngle.Reflex ? 1 : 0);
                assert(bv.Id == v.ID);
                if (bv.degree == 1)
                {
                    NumExtraBevelingVertices++;
                }
                Vertices.Add(bv);
            }






            var tt = mesh.Triangles.Where(t => t.GetHashCode() >= 0);



            void assign_halfedge(BasicTriangle t, int a, int b)
            {
                bool oposite = false;
                if (a > b)
                {
                    oposite = true;

                    (a, b) = (b, a);
                }

                TrianglesOnEdge value;
                if (!halfedges.TryGetValue((a, b), out value))
                {
                    value = new TrianglesOnEdge();
                    halfedges.Add((a, b), value);
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






            }

            void assign_neighbors(BasicTriangle t, int idx, int a, int b, BasicTriangle?[] neighbors)
            {
                assert(a < Vertices.Count);

                var cva = contourVertices[a];
                var cvb = contourVertices[b];

                if (cva.NextInLAV.ID != b)
                {
                    bool oposite = false;

                    if (a > b)
                    {
                        oposite = true;

                        (a, b) = (b, a);
                    }

                    TrianglesOnEdge value;
                    if (!halfedges.TryGetValue((a, b), out value))
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
                var t = new BasicTriangle(Triangles.Count);
                assert(t.Id == Triangles.Count);
                Triangles.Add(t);
                var a = tri[0].a;
                var b = tri[1].a;
                var c = tri[2].a;

                //var wavefrontedges = new WavefrontEdge?[3];
                t.Vertices[0] = Vertices[a];
                t.Vertices[1] = Vertices[b];
                t.Vertices[2] = Vertices[c];

                assign_halfedge(t, a, b);
                assign_halfedge(t, b, c);
                assign_halfedge(t, c, a);

                //t.set_wavefronts(wavefrontedges[0], wavefrontedges[1], wavefrontedges[2]);
                //foreach (var we in wavefrontedges)
                //{
                //    if (we != null)
                //    {
                //        we.set_initial_incident_triangle(t);
                //    }
                //}



            }

            foreach (var t in Triangles)
            {

                var neighbors = t.Neighbors.Array;
                assign_neighbors(t, 0, t.Vertices[1].Id, t.Vertices[2].Id, neighbors);
                assign_neighbors(t, 1, t.Vertices[2].Id, t.Vertices[0].Id, neighbors);
                assign_neighbors(t, 2, t.Vertices[0].Id, t.Vertices[1].Id, neighbors);

            }

        }

        internal int get_total_degree()
        {
            return halfedges.Count * 2;
        }

    }
   
}