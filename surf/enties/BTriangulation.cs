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
        //        vertex.id = v.id;
        //        vertex.label = v.label;

        //        vertices[v.id] = vertex;
        //    }

        //    // Maps a triangle to its 3 edges (used to set next pointers).
        //    var map = new List<HalfEdge>[mesh.triangles.Count];

        //    Face face;

        //    foreach (var t in mesh.triangles)
        //    {
        //        face = new Face(null);
        //        face.id = t.id;

        //        faces[t.id] = face;

        //        map[t.id] = new List<HalfEdge>(3);
        //    }

        //    Otri tri = default(Otri), neighbor = default(Otri);
        //    TVertex org, dest;

        //    int id, nid, count = mesh.triangles.Count;

        //    HalfEdge edge, twin, next;

        //    var edges = dcel.HalfEdges;

        //    // Count half-edges (edge ids).
        //    int k = 0;

        //    // Maps a vertex to its leaving boundary edge.
        //    var boundary = new Dictionary<int, HalfEdge>();

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
        //                twin = new DcelHalfEdge(vertices[dest.id], nid < 0 ? Face.Empty : faces[nid]);

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

        //        if (edge.twin.origin.id == next.origin.id)
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
        //        e.next = boundary[e.twin.origin.id];
        //    }

        //    dcel.Vertices.AddRange(vertices);
        //    dcel.Faces.AddRange(faces);

        //    return dcel;
        //}
    }
}
