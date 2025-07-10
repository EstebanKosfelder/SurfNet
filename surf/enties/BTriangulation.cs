namespace SurfNet
{
    public class BTriangulation
    {


        //public static DcelMesh ToDCEL(Mesh mesh)
        //{
        //    var dcel = new DcelMesh();

        //    var vertices = new DcelVertex[mesh.vertices.Count];
        //    var faces = new Face[mesh.triangles.Count];

        //    dcel.HalfEdges.Capacity = 2 * mesh.NumberOfEdges;

        //    mesh.Renumber();

        //    HVertex vertex;

        //    foreach (var v in mesh.vertices.Values)
        //    {
        //        vertex = new HVertex(v.X, v.Y);
        //        vertex.Id = v.Id;
        //        vertex.label = v.label;

        //        vertices[v.Id] = vertex;
        //    }

        //    // Maps a triangle to its 3 edges (used to set next pointers).
        //    var map = new List<HalfEdge>[mesh.triangles.Count];

        //    Face face;

        //    foreach (var t in mesh.triangles)
        //    {
        //        face = new Face(null);
        //        face.Id = t.Id;

        //        faces[t.Id] = face;

        //        map[t.Id] = new List<HalfEdge>(3);
        //    }

        //    Otri tri = default(Otri), neighbor = default(Otri);
        //    TVertex org, dest;

        //    int Id, nid, count = mesh.triangles.Count;

        //    HalfEdge edge, twin, next;

        //    var edges = dcel.HalfEdges;

        //    // Count half-edges (edge ids).
        //    int k = 0;

        //    // Maps a vertex to its leaving boundary edge.
        //    var boundary = new Dictionary<int, HalfEdge>();

        //    foreach (var t in mesh.triangles)
        //    {
        //        Id = t.Id;

        //        tri.tri = t;

        //        for (int i = 0; i < 3; i++)
        //        {
        //            tri.orient = i;
        //            tri.Sym(ref neighbor);

        //            nid = neighbor.tri.Id;

        //            if (Id < nid || nid < 0)
        //            {
        //                face = faces[Id];

        //                // Get the endpoints of the current triangle edge.
        //                org = tri.Org();
        //                dest = tri.Dest();

        //                // Create half-edges.
        //                edge = new DcelHalfEdge(vertices[org.Id], face);
        //                twin = new DcelHalfEdge(vertices[dest.Id], nid < 0 ? Face.Empty : faces[nid]);

        //                map[Id].Add(edge);

        //                if (nid >= 0)
        //                {
        //                    map[nid].Add(twin);
        //                }
        //                else
        //                {
        //                    boundary.Add(dest.Id, twin);
        //                }

        //                // Set leaving edges.
        //                edge.origin.Leaving = edge;
        //                twin.origin.Leaving = twin;

        //                // Set twin edges.
        //                edge.twin = twin;
        //                twin.twin = edge;

        //                edge.Id = k++;
        //                twin.Id = k++;

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

        //        if (edge.twin.origin.Id == next.origin.Id)
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
        //        e.next = boundary[e.twin.origin.Id];
        //    }

        //    dcel.Vertices.AddRange(vertices);
        //    dcel.Faces.AddRange(faces);

        //    return dcel;
        //}
    }
}
