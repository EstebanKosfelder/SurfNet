namespace SurfNet
{
    using static DebugLog;
    using static Mathex;
    using TriangleOriginalVertexIndexList = List<int>;
    public partial class KineticTriangulation
    {


        public int ccw(int i) => i.ccw();
        public int cw(int i) => i.cw();

        public int mod3(int i) => i.mod3();


        internal partial int get_num_initial_triangles(DcelMesh mesh)
        {
            int num_initial_triangles;
            //if (restrict_component_ < 0) {
            //    num_initial_triangles = ct.tds().number_of_faces();
            //} else {
            num_initial_triangles = 0;
            foreach (var fit in mesh.Faces)
            {
                num_initial_triangles += 1;// (fit.info().component == restrict_component_);
            }
            //}
            return num_initial_triangles;
        }

        ///** initialize the TDS, the Triangulation Data Structure
        // *
        // * After this function, we'll have Kinetic Triangles with neighborhood
        // * references but no wavefront vertices yet.
        // *
        // * Given a basic constrained triangulation (CT) of the input (point and edge
        // * set), start building up our kinetic triangulation.
        // *
        // * We begin with allocating sufficient memory to hold all wavefront vertices
        // * and triangles during the entire straight-skeleton wavefront propagation.
        // *
        // * Next, we make kinetic triangles, one for each of the basic constrained
        // * triangulation.  These don't have vertex references initially as the wavefront
        // * vertices don't exist yet; those will come later.  The only exception is that
        // * references to the infinite vertex are set already.  We do however, create a map
        // * that associates a CT face with a kinetic triangle list index.
        // *
        // * Then, we set up neighborhood references between adjacent triangles that are
        // * not seperated by a constraint.
        // */

        //private partial void initialize_tds(DcelMesh mesh)
        //{
        //    //DBG_FUNC_BEGIN(//DBG_KT_SETUP);

        //    // build kinetic triangulation
        //    int num_t = mesh.Faces.Count;
        //    // at each input vertex we start one kinetic vertex per sector, and the number of sectors equals the vertex degree.
        //    // Additionally, we'll split vertices for beveling (including degree-1 vertices)

        //    int num_events = num_t;
        //    var num_v = mesh.HalfEdges.Count;

        //    /** Allocate sufficient space for our arrays */
        //    this.vertices = new WavefrontVertexList();
        //    this.vertices.Capacity = (num_v + num_t);
        //    triangles.Capacity = (num_t);


        //    ///   vertices.Add(WavefrontVertex.make_infinite_vertex());
        //    WavefrontVertex infinite = null; /// vertices.back();

        //    /** Create kinetic triangles.  Without any vertices initially.
        //     * The only exception is that we already point to the correct infinite vertex.
        //     * Also set up the CT-face to KT index map.
        //     */
        //    //DBG(//DBG_KT_SETUP) << "#faces: " << num_t;

        //    foreach (var fit in mesh.Faces)
        //    {
        //        Debug.Assert(triangles.Count == fit.Id);
        //        var t = new KineticTriangle(triangles.Count,/* fit.ID*/ 0 );
        //        triangles.Add(t);

        //        //int i = 0;
        //        //foreach (var he in fit.EnumerateEdges())
        //        //{
        //        //    triangle_original_vertex_indices.Add(he.Origin.ID);

        //        //}
        //        //assert(triangles.size() * 3 == triangle_original_vertex_indices.size());
        //        //face_to_triangle_idx.emplace((fit, triangles.size() - 1));
        //        ////DBG(//DBG_KT_SETUP) << "Added " << &triangles.back()
        //        ////  << "; input vidx: "
        //        ////  << fit.vertex(0).info().original_vertex_idx << ","
        //        ////  << fit.vertex(1).info().original_vertex_idx << ","
        //        ////  << fit.vertex(2).info().original_vertex_idx;
        //    }
        //    /** Set up neighborhood references between triangles that are not
        //     * seperated by a constraint.
        //     */
        //    foreach (var fit in mesh.Faces)
        //    {
        //        var fedges = fit.EnumerateEdges().ToArray();
        //        Debug.Assert(fedges.Length == 3);

        //        // if (!fit.info().matches_component(restrict_component_)) continue;

        //        KineticTriangle?[] n = fedges.Select(e => (e.Twin.IsContour) ? null : triangles[e.Twin.Face.Id]).ToArray();

        //        var t = triangles[fit.Id];
        //        //DBG(//DBG_KT_SETUP) << "set neighbors for " << t << " to "
        //        //<< n[0] << ", "
        //        //<< n[1] << ", "
        //        //<< n[2] << ".";
        //        t.set_neighbors(n);
        //    }

        //    //   assert(triangles.size() == num_initial_triangles);
        //    //DBG_FUNC_END(//DBG_KT_SETUP);
        //}

        //internal partial void create_supporting_lines(DcelMesh mesh)
        //{
        //    // DBG_FUNC_BEGIN(DBG_KT_SETUP);
        //    /* total number of wavefront edges during the entire propragation.
        //     * this is an upper bound.
        //     *
        //     * num_t is a really rough upper bound on the number of split events.
        //     */
        //    // int num_t = num_initial_triangles + input.get_num_extra_beveling_vertices();
        //    //  var num_wavefront_edges = input.edges().size() * 2 +  num_t * 2;
        //    this.wavefront_edges = new WavefrontEdgeList();
        //    //  this.wavefront_edges.Capacity = num_wavefront_edges;

        //    foreach (var fit in mesh.Faces)
        //    {
        //        //if (!fit.info().matches_component(restrict_component_)) continue;
        //        WavefrontEdge3 w = new WavefrontEdge3();

        //        KineticTriangle t = triangles[fit.Id];
        //        var edges = fit.EnumerateEdges().ToArray();

        //        Debug.Assert(edges.Length == 3);
        //        for (int i = 0; i < 3; ++i)
        //        {
        //            if (edges[i].is_constrained())
        //            {
        //                /*
        //                int vidxu = get_basic_vertex_idx_from_triangle_vertex_indices(input, triangle_original_vertex_indices, triangle_idx, ccw(i));
        //                int vidxv = get_basic_vertex_idx_from_triangle_vertex_indices(input, triangle_original_vertex_indices, triangle_idx, cw(i));
        //                BasicVertex u = input.vertices()[vidxu];
        //                BasicVertex v = input.vertices()[vidxv];
        //                */

        //                var vidxu = edges[ccw(i)].origin.ID;
        //                var vidxv = edges[cw(i)].origin.ID;

        //                var u = edges[ccw(i)].origin;
        //                var v = edges[cw(i)].origin;

        //                //if (!input.has_edge(vidxu, vidxv))
        //                //{
        //                //    //LOG(ERROR) << "Cannot find edge for triangulation constraint (" << vidxu << ", " << vidxv << ") in input.  Input is not a PSLG.";
        //                //    //LOG(ERROR) << " v" << vidxu << ": " << CGAL_point(input.vertices()[vidxu].p);
        //                //    //LOG(ERROR) << " v" << vidxv << ": " << CGAL_point(input.vertices()[vidxv].p);
        //                //    //LOG(ERROR) << " Probably one of these points is on the interior of another edge.";
        //                //    //exit(1);
        //                //}
        //                var edge = mesh.HalfEdges.FirstOrDefault(e => e.Origin.ID == vidxu && e.Twin.Origin.ID == vidxv);

        //                //   DBG(DBG_KT_SETUP) << "at " << t << "/" << i << ": vertices are " << vidxu << " and " << vidxv;

        //                WavefrontEdge buddy_wavefront_edge = null;

        //                var ct_n_hdl = edges[i].Twin.Face;
        //                if (ct_n_hdl == null || ct_n_hdl.Id < 0)
        //                {
        //                    /* neighbor not in map */
        //                    assert(ct_n_hdl?.Id < 0);
        //                }
        //                else
        //                {
        //                    KineticTriangle n = triangles[ct_n_hdl.Id];
        //                    //   assert(ct_n_hdl.info().matches_component(restrict_component_));

        //                    int idx_in_n = Array.IndexOf(ct_n_hdl.EnumerateEdges().Select(e => e.Twin.Face).ToArray(), fit);
        //                    Debug.Assert(idx_in_n > 0 && idx_in_n < 3);
        //                    buddy_wavefront_edge = n.wavefronts[idx_in_n]; /* may still be null, and that's OK */
        //                }

        //                SkeletonDCELFace skeleton_dcel_face = skeleton.setup_new_input_edge(buddy_wavefront_edge /* may be NULL and that's OK */);

        //                w[i] = new WavefrontEdge(u.DcelPoint, v.DcelPoint, /*edge.weight*/0, t, skeleton_dcel_face);
        //                wavefront_edges.Add(w[i]);
        //            }
        //            else
        //            {
        //                w[i] = null;
        //            }
        //        }

        //        t.set_wavefronts(w[0], w[1], w[2]);
        //    }
        //    // DBG_FUNC_END(DBG_KT_SETUP);
        //}

        /** Create wavefront edges for all constraints and references them from kinetic triangles.
         *
         * Also create straight skeleton (DCEL) faces and halfedge pairs for them all.
         */

        /** Return a BasicVertex index from a triangle index and vertex position within the triangle.
         *
         * This only works for triangles that were part of the initial triangulation
         * and where we set up the mapping in triangle_original_vertex_indices.
         *
         * So it gets invalidated partially as soon as we start splitting triangles.
         *
         * This function primary checks bounds and some basic validity.
         */

        private static partial int get_basic_vertex_idx_from_triangle_vertex_indices(
          BasicInput input,
          TriangleOriginalVertexIndexList triangle_original_vertex_indices,
          int t_idx,
          int i
        )
        {
            assert(i < 3);

            /* Make sure we do not overrun the original triangulation size.
             */
            int idx = t_idx * 3 + i;
            assert(idx < triangle_original_vertex_indices.size());
            assert(triangle_original_vertex_indices[idx] >= 0);
            int vidx = triangle_original_vertex_indices[idx];
            assert(vidx < input.vertices().size());
            return vidx;
        }

        /** Return a BasicVertex from a triangle index and vertex position within the triangle.
         *
         * See get_basic_vertex_idx_from_triangle_vertex_indices() for details.
         */

        private static partial BasicVertex get_basic_vertex_from_triangle_vertex_indices(
          BasicInput input,
          TriangleOriginalVertexIndexList triangle_original_vertex_indices,
          int t_idx,
          int i
        )
        {
            BasicVertex bv = input.vertices()[get_basic_vertex_idx_from_triangle_vertex_indices(input, triangle_original_vertex_indices, t_idx, i)];
            return bv;
        }

        /** mark the original vertex as unusable.
         *
         * once we have created a kinetic vertex for a basic (input) vertex,
         * we no longer want to touch the basic vertex from input.vertices().
         *
         * This function marks the entry in the triangle-vertex to original-vertex map
         * as invalid, preventing further accesses.
         */

        [Obsolete]
        private static partial void invalidate_basic_vertex_idx_in_triangle_vertex_indices(
          BasicInput input,
          TriangleOriginalVertexIndexList triangle_original_vertex_indices,
          int t_idx,
          int i
        )
        {
            throw new NotImplementedException();
            int idx = t_idx * 3 + i;
            assert(idx < triangle_original_vertex_indices.size());
            assert(triangle_original_vertex_indices[idx] >= 0);
            triangle_original_vertex_indices[idx] = -1;
        }

        /** Create wavefront vertices, except those that need beveling (such as degree-1 vertices)
         */

        private partial void create_kinetic_vertices(DcelMesh mesh)
        {
            ////DBG_FUNC_BEGIN(//DBG_KT_SETUP);

            ///* We iterate over all triangles, and in each triangle over each
            // * vertex v.  If the triangle is the start of a clockwise triangle
            // * fan about v, that is, if the correct edge incident to v is contraint,
            // * we set up the kinetic vertex.
            // */
            //foreach (var t_it in triangles)
            //{
            //    /* Invariant: all previous triangles have vertices set up on
            //     * the "head" side of each constraint they have.
            //     *
            //     * They not necessarily have vertices set up on the "tail" side of
            //     * each constraint or at vertices with no incident (in this triangle)
            //     * constraint.
            //     */
            //    int t_idx = t_it.id;
            //    //assert(t_idx == t_it - triangles.begin());
            //    //DBG(//DBG_KT_SETUP) << "setting up kinetic vertices for " << &*t_it;
            //    for (int i = 0; i < 3; ++i)
            //    {
            //        //DBG(//DBG_KT_SETUP) << "vertex at idx " << i;
            //        if (!t_it.is_constrained(ccw(i))) continue;

            //        /* the edge cw of i (opposite i's ccw vertex) is constrained. */
            //        AroundVertexIterator faces_it = incident_faces_iterator(t_it, i);
            //        assert(!faces_it.next_triangle_ccw());
            //        AroundVertexIterator most_cw_triangle = faces_it.most_cw();

            //        WavefrontEdge l = t_it.wavefront(ccw(i));
            //        WavefrontEdge r = most_cw_triangle.t().wavefront(cw(most_cw_triangle.v_in_t_idx()));

            //        BasicVertex bv = get_basic_vertex_from_triangle_vertex_indices(input, triangle_original_vertex_indices, t_idx, i);

            //        if (bv.degree == 1)
            //        {
            //            continue; // Need to split/bevel in any case.
            //        }

            //        assert(bv.degree > 1);

            //        // XXX do not make a vertex if this is reflex and the reflex_beveling_add is greater than 1.
            //        WavefrontVertex v = vertices.make_initial_vertex(bv.p, l, r);

            //        // Should not have created a vertex that needed beveling.
            //        assert(bv.reflex_beveling_add == 0 || v.is_convex_or_straight());

            //        for (AroundVertexIterator it = incident_faces_iterator(t_it, i); it != incident_faces_end(); it = it.walk_dir_cw)
            //        {
            //            it.set_vertex(it.v_in_t_idx(), v);
            //            invalidate_basic_vertex_idx_in_triangle_vertex_indices(input, triangle_original_vertex_indices, (*it).id, it.v_in_t_idx());
            //            //DBG(//DBG_KT_SETUP) << "  setting vertex to " << v << " in " << &*it << "(" << it.v_in_t_idx() << ")";
            //        };

            //        //DBG(//DBG_KT_SETUP) << "  setting vertex 1 to " << v << " in " << *l;
            //        //DBG(//DBG_KT_SETUP) << "  setting vertex 0 to " << v << " in " << *r;
            //        assert(l.vertex(1) == v);
            //        assert(r.vertex(0) == v);
            //    }
            //}
            //DBG_FUNC_END(//DBG_KT_SETUP);
        }

        /** Split vertex
         *
         * For creating bevels (such as when dealing with degree-1 vertices),
         * we need to split triangulation vertices.
         *
         * This function splits a vertex v, as given by a triangle t and vertex index.
         * The triangle is duplicated, and the new triangle is ccw of t at v.
         *
         * Returns a pointer to the new KineticTriangle.  The new edges is between
         * vertices 0 and 1.
         */

        private partial KineticTriangle split_vertex(
          KineticTriangle t,
          int i,
          WavefrontEdge new_e,
          WavefrontVertex new_v
        )
        {
            //DBG_FUNC_BEGIN(//DBG_KT_SETUP);
            assert(new_e != null);
            assert(new_v != null);
            assert(t != null);

            /* Split the triangle */
            triangles.Add(new KineticTriangle(triangles.size(), t.component));
            assert(triangles.back().Id == triangles.size() - 1);
            KineticTriangle new_t = triangles.back();
            //DBG(//DBG_KT_SETUP) << " New triangle " << &new_t;

            new_t.wavefronts[0] = null;
            new_t.wavefronts[1] = new_e;
            new_t.wavefronts[2] = t.wavefronts[cw(i)];
            new_t.neighbors[0] = t;
            new_t.neighbors[1] = null;
            new_t.neighbors[2] = t.neighbors[cw(i)];

            assert(new_t.vertices[1] == null);
            assert(new_t.vertices[2] == null);
            new_t.set_vertex(0, new_v);
            if (t.vertices[ccw(i)] != null) { new_t.set_vertex(1, t.vertices[ccw(i)]); }
            if (t.vertices[i] != null) { new_t.set_vertex(2, t.vertices[i]); }
            //DBG(//DBG_KT_SETUP) << "  setting vertex to " << new_v << " in " << &new_t << "(" << 0 << ")";

            t.wavefronts[cw(i)] = null;
            t.neighbors[cw(i)] = new_t;

            assert((new_t.wavefronts[2] == null) != (new_t.neighbors[2] == null));
            if (new_t.wavefronts[2] != null)
            {
                new_t.wavefronts[2].set_incident_triangle(new_t);
                //DBG(//DBG_KT_SETUP) << "  setting incident triangle for wavefronts[2], " << *new_t.wavefronts[2];
            }
            else
            {
                int pos = new_t.neighbors[2].index(t);
                new_t.neighbors[2].neighbors[pos] = new_t;
                //DBG(//DBG_KT_SETUP) << "  setting neighbor triangle for neighbors[2], " << new_t.neighbors[2];
            }

            if (new_e.incident_triangle() == null)
            {
                new_e.set_incident_triangle(new_t);
                //DBG(//DBG_KT_SETUP) << "  setting incident triangle for current_edge, " << *new_e;
            }

            //DBG(//DBG_KT_SETUP) << " Old triangle details: " << t;
            //DBG(//DBG_KT_SETUP) << " New triangle details: " << &new_t;

            //DBG_FUNC_END(//DBG_KT_SETUP);
            return new_t;
        }

        /** return (position, is_infinite) of vertex i of triangle t.
         *
         * This uses the kinetic triangulation's vertex location
         * if already set, and falls back to using the vertex
         * position from the original, underlying, constrained triangulation.
         *
         * The latter only works for triangles that were created initially,
         * not any that are the result of splits (but those should all have
         * kinetic vertices yet).  NO, WRONG, XXX NOT GUARANTEED.
         */

        private partial (Point2, bool) get_vertex_pos(BasicInput input,
          TriangleOriginalVertexIndexList triangle_original_vertex_indices,
          KineticTriangle t,
          int i
        )
        {
            ////DBG_FUNC_BEGIN(//DBG_KT_SETUP);
            //DBG_INDENT_INC();

            WavefrontVertex kv = t.vertex(i);

            Point2 pos = Point2.NaN;
            bool is_inf = false;

            if (kv != null)
            {
                //DBG(//DBG_KT_SETUP) << "  v is " << kv;
                if (kv.is_infinite)
                { /* we already have a vertex set at t's cw */
                    is_inf = true;
                }
                else
                {
                    pos = kv.pos_zero;
                }
            }
            else
            {
                //DBG(//DBG_KT_SETUP) << "  kinetic vertex is not yet set up, using vertex from input triangulation.";
                /* The infinite vertex is associated in kinetic triangles already
                 * in initialize_tds at the very start of setting things up.
                 * So if we are here, it's not the infinite vertex.
                 *
                 * However, it may be a vertex we'll have to bevel later, so it might
                 * still not be set in the kinetic triangulation.  In that case,
                 * get its position from the underlying original triangulation and input.
                 */
                assert(t.Id < triangle_original_vertex_indices.size() / 3); // XXX -- check if this holds
                BasicVertex bv = get_basic_vertex_from_triangle_vertex_indices(input, triangle_original_vertex_indices, t.Id, i);
                pos = bv.p;
            }

            //DBG_INDENT_DEC();
            ////DBG_FUNC_END(//DBG_KT_SETUP);
            return (pos, is_inf);
        }

        /** Create bevels at one vertex.
         *
         * .) find the neighboring constraints,
         * .) create beveling wavefront supporting lines
         * .) create the kinetic vertices between adjacent wavefronts, and
         * .) split one or more triangles and set the vertices accordingly.
         */

        private partial void create_bevels_at_vertex(DcelMesh mesh, KineticTriangle t, int i)

        // TODO [AroundVertexIterator y punteros] Reimplementar

        {
            throw new NotImplementedException();

            //   //DBG_FUNC_BEGIN(//DBG_KT_SETUP);
            //   /* We should not have any triangles with undefined vertices after
            //    * we went past the original size.  In particular, since the
            //    * triangle_vertex_handle_idx array is not defined for the new ones,
            //    * and in fact may have become invalid for some of the existing entries
            //    * too.
            //    */
            //   //DBG(//DBG_KT_SETUP) << t << "; vertex " << i;
            //   WavefrontEdge l, r;

            //   AroundVertexIterator most_ccw_triangle = incident_faces_iterator(t, i).most_ccw();
            //   AroundVertexIterator most_cw_triangle = incident_faces_iterator(t, i).most_cw();
            //   l = most_ccw_triangle.t().wavefront(ccw(most_ccw_triangle.v_in_t_idx()));
            //   r = most_cw_triangle.t().wavefront(cw(most_cw_triangle.v_in_t_idx()));

            //   //DBG(//DBG_KT_SETUP) << " incident edges: " << *l << ", " << *r;

            //   /* Build the list of wavefront edges */
            //   /*************************************/
            //   List<WavefrontEdge> edges = new List<WavefrontEdge>();
            //   edges.Add(l);
            //   BasicVertex bv = get_basic_vertex_from_triangle_vertex_indices(input, triangle_original_vertex_indices, t.id, i);

            //   if (bv.degree == 1)
            //   {
            //       if (bv.reflex_beveling_add >= 2)
            //       {
            //           throw new NotImplementedException("Beveling degree-one vertex with 2 or more extra vertices not implemented yet");
            //       }
            //       {
            //           assert(l.l().l == r.l().l.opposite());

            //           if (l.l().weight != r.l().weight)
            //           {
            //               throw new NotImplementedException("Unclear what to do when the incidents weights do not match in beveling");

            //           }

            //           WavefrontEdge w = new WavefrontEdge(new WavefrontSupportingLine(l.l().l.opposite().perpendicular(bv.p), l.l().weight));
            //           wavefront_edges.Add(w);
            //           //DBG(//DBG_KT_SETUP) << " New wavefront edge: " << *w;
            //       }
            //   }
            //   else
            //   {
            //       throw new NotImplementedException("Beveling reflex vertices not implemented yet");
            //   }
            //   edges.Add(r);

            //   List<WavefrontVertex> new_wavefront_vertices = new List<WavefrontVertex>();
            //   assert(edges.size() > 2);

            //   //DBG(//DBG_KT_SETUP) << " Edges:";
            //   foreach (var e in edges)
            //   {
            //       //DBG(//DBG_KT_SETUP) << "  - " << *e;
            //   }

            //   /* Build the list of wavefront vertices */
            //   /***************************************/
            ////   var current_edge = edges.First();
            //   var previous_edge = edges.First();
            //   WavefrontVertex prev_vertex = null;
            //   //DBG(//DBG_KT_SETUP) << " New wavefront vertices:";
            //   foreach (var current_edge in edges.Skip(1)) /*; current_edge != edges.end(); previous_edge = current_edge++*/ {
            //       WavefrontVertex v = vertices.make_initial_vertex(bv.p, previous_edge, current_edge, true);
            //       if (prev_vertex != null)
            //       {
            //           v.link_tail_to_tail(prev_vertex);
            //       }
            //       prev_vertex = v;

            //       new_wavefront_vertices.Add(v);
            //       previous_edge = current_edge;
            //       //DBG(//DBG_KT_SETUP) << "  - " << v;
            //   }

            //   /* And update the triangulation, splitting vertices */
            //   /****************************************************/
            //   var new_v = new_wavefront_vertices.FirstOrDefault();
            //   var currentEdge = edges[1]; //.First();
            //   //++current_edge;
            //   var t_it = most_ccw_triangle;
            //   while (t_it != incident_faces_end())
            //   {
            //       int t_idx = t_it.t().id;
            //       assert(t_it.t() == triangles[t_idx]);

            //       /* if we split in the iteration before, this is a duplicate (as we set it in split_vertex()), but who cares :) */
            //       t_it.t().set_vertex(t_it.v_in_t_idx(), new_v);
            //       //DBG(//DBG_KT_SETUP) << " Setting vertex to " << (*new_v) << " in " << t_it.t() << "(" << t_it.v_in_t_idx() << ")";

            //       /* If the normal of the right edge points into the current triangle,
            //        * we move to the next wavefront/vertex.  We split the triangle by
            //        * duplicating the triangulation vertex here, assigning the previous
            //        * and new wavefront vertex to each of the triangulation vertices from
            //        * the split.
            //        */
            //       bool split_this;
            //       bool need_flip = false;
            //       assert(currentEdge == (new_v).incident_wavefront_edge(1));

            //       if ((new_v) == (new_wavefront_vertices.back()))
            //       {
            //           //DBG(//DBG_KT_SETUP) << "  At last vertex, nothing to split anymore";
            //           split_this = false;
            //       }
            //       else if (t_it.next_triangle_cw() == null)
            //       {
            //           //DBG(//DBG_KT_SETUP) << "  At last triangle, we need to split here.";
            //           split_this = true;
            //       }
            //       else
            //       {
            //           //DBG(//DBG_KT_SETUP) << " Checking if we need to split in " << t_it.t();
            //           Point_2 pos_t_ccw;
            //           bool infinite_ccw_vertex;
            //           (pos_t_ccw, infinite_ccw_vertex) = get_vertex_pos(input, triangle_original_vertex_indices, t_it.t(), ccw(t_it.v_in_t_idx()));

            //           if (!infinite_ccw_vertex)
            //           {
            //               //DBG(//DBG_KT_SETUP) << "  All nice and finite";
            //               Point_2 pos_plus_normal = (bv.p +(Point_2) (currentEdge).l().normal_direction.perpendicular(RIGHT_TURN));

            //               //DBG(//DBG_KT_SETUP) << "   p     " << bv.p;
            //               //DBG(//DBG_KT_SETUP) << "   p+    " << pos_plus_normal;
            //               //DBG(//DBG_KT_SETUP) << "   ptccw " << pos_t_ccw;
            //               var orientation = Tools.orientation(bv.p, pos_plus_normal, pos_t_ccw);
            //               //DBG(//DBG_KT_SETUP) << "   orientation " << orientation;
            //               split_this = (orientation == LEFT_TURN);
            //           }
            //           else
            //           {
            //               //DBG(//DBG_KT_SETUP) << "  Unbounded triangle: ccw vertex is the infinite vertex";
            //               //WavefrontVertex const * const cw_v = t_it.t().vertex(cw(t_it.v_in_t_idx()));
            //               ////DBG(//DBG_KT_SETUP) << "   cw_vertex point              " << cw_v;

            //               var next_t = t_it;
            //               ++next_t;
            //               Point_2 next_pos_t_ccw;
            //               bool next_infinite_ccw_vertex;
            //               (next_pos_t_ccw, next_infinite_ccw_vertex) = get_vertex_pos(input, triangle_original_vertex_indices, next_t.t(), ccw(next_t.v_in_t_idx()));
            //               assert(!next_infinite_ccw_vertex);

            //               //DBG(//DBG_KT_SETUP) << "   current_edge supporting line " << CGAL_line((*current_edge).l().l);
            //               //DBG(//DBG_KT_SETUP) << "   next_pos_t_ccw " << next_pos_t_ccw;

            //               /* We can err on the side of having triangulation spokes on the CH boundary,
            //                * or we can have collinear points on the CH boundary.
            //                * Whichever we pick, we will then have to figure out what the right way
            //                * to flip is in one case or another.
            //                */
            //               // split_this = (*current_edge).l().l.has_on_negative_side(next_pos_t_ccw);
            //               split_this = !(current_edge).l().l.has_on_positive_side(next_pos_t_ccw);
            //               if (split_this)
            //               {
            //                   Point_2 pos_t_cw;
            //                   bool infinite_cw_vertex;
            //                   (pos_t_cw, infinite_cw_vertex) = get_vertex_pos(input, triangle_original_vertex_indices, t_it.t(), cw(t_it.v_in_t_idx()));
            //                   assert(!infinite_cw_vertex);

            //                   //DBG(//DBG_KT_SETUP) << "   has on pos " << (*current_edge).l().l.has_on_positive_side(pos_t_cw);
            //                   //DBG(//DBG_KT_SETUP) << "   has on neg " << (*current_edge).l().l.has_on_negative_side(pos_t_cw);

            //                   /* The triangle created by our splitting the edge from v to ccw(v) will be incident
            //                    * to the infinite vertex, but the old vertex will not be on the convex hull,
            //                    * so we need to split things.
            //                    */
            //                   need_flip = (current_edge).l().l.has_on_positive_side(pos_t_cw);
            //                   //DBG(//DBG_KT_SETUP) << "   need flip " << need_flip;
            //               }
            //           }
            //       }

            //       if (split_this)
            //       {
            //           //DBG(//DBG_KT_SETUP) << " Splitting triangle " << t_it.t();
            //           assert(t_it.t().vertex(t_it.v_in_t_idx()) == new_v);
            //           ++new_v;
            //           //DBG(//DBG_KT_SETUP) << "  next vertex is " << (*new_v);

            //           KineticTriangle new_t = split_vertex(t_it.t(), t_it.v_in_t_idx(), current_edge, *new_v);
            //           assert((current_edge).incident_triangle() == new_t || &(current_edge) == &(edges.back()));

            //           if (need_flip)
            //           {
            //               //DBG(//DBG_KT_SETUP) << " need flip.";
            //               //DBG(//DBG_KT_SETUP) << "  t  " << t_it.t();
            //               //DBG(//DBG_KT_SETUP) << "  tn " << new_t;
            //               LOG(WARNING) << __FILE__ << ":" << __LINE__ << " " << "untested code path.";
            //               new_t.do_raw_flip_inner(0);
            //               //DBG(//DBG_KT_SETUP) << "  t  " << t_it.t();
            //               //DBG(//DBG_KT_SETUP) << "  tn " << new_t;
            //           };

            //           t_it = AroundVertexIterator(new_t, 0);
            //           ++current_edge;
            //       }
            //       else
            //       {
            //           //DBG(//DBG_KT_SETUP) << " Moving on, not splitting this triangle " << t_it.t();
            //           ++t_it;
            //       }
            //   }
            //   assert((new_v) == (new_wavefront_vertices.back()));
            //   //DBG_FUNC_END(//DBG_KT_SETUP);
        }

        /** Create bevels
         *
         * Split vertices for bevels at degree-1 vertices and where we want other reflex vertices to bevel
         *
         * For each vertex in the triangulation that is not yet initialized, call create_bevels_at_vertex()
         */

        private partial void create_bevels(DcelMesh mesh)
        {
            //DBG_FUNC_BEGIN(//DBG_KT_SETUP);

            int initial_triangles_size = triangles.size();

            /* verify all existing wavefronts have an incident triangle that
             * reference back to the wavefront.
             */
            foreach (var wf in wavefront_edges)
            {
                assert(wf.incident_triangle() != null);
                assert(wf.incident_triangle().has_wavefront(wf));
            }

            /* set up bevels */
            var ft_id = 0;
            foreach (var t_it in triangles)
            {
                assert(t_it.Id == ft_id);
                ft_id++;
                for (int i = 0; i < 3; ++i)
                {
                    if (t_it.vertex(i) != null) continue;
                    assert(t_it.Id < initial_triangles_size);

                    //DBG(//DBG_KT_SETUP) << "Doing bevels at t " << &*t_it << "; vertex " << i;
                    create_bevels_at_vertex(mesh, t_it, i);
                }
            }

            /* assert all kinetic triangulation vertices have an assigned
             * wavefront vertex.
             */
            foreach (var t_it in triangles)
            {
                int t_idx = t_it.Id;
                assert(t_it == triangles[t_it.Id]);
                for (int i = 0; i < 3; ++i)
                {
                    assert(t_it.vertex(i) != null);
                }
            }

            /* verify all existing wavefronts have an incident triangle that
             * reference back to the wavefront.
             */
            foreach (var wf in wavefront_edges)
            {
                //DBG(//DBG_KT_SETUP) << " wf: " << wf;
                assert(wf.incident_triangle() != null);
                assert(wf.incident_triangle().has_wavefront(wf));
            }

            //DBG_FUNC_END(//DBG_KT_SETUP);
        }

        private partial void store_initial_wavefront_vertices()
        {
            foreach (var e in wavefront_edges)
            {
                assert(e.is_initial);
                e.set_initial_vertices();
            }
        }

        

       

       

        



        public partial void initialize(
            DcelMesh mesh,
            WavefrontEdgeList p_wavefront_edges,
            int p_restrict_component)
        {
            //DBG_FUNC_BEGIN(//DBG_KT_SETUP);

            //assert(!initialized);
            //wavefront_edges = p_wavefront_edges;
            //restrict_component_ = p_restrict_component;

            ///* set up basic triangulation */
            //BasicTriangulation ct = new BasicTriangulation();
            //ct.initialize(input);

            ////DBG(//DBG_KT_SETUP) << "Input has " << ct.max_component() + 1 << " component(s): 0 .. " << ct.max_component();
            ////DBG(//DBG_KT_SETUP) << "Requested straight skeleton of component " << restrict_component_;

            //if (ct.max_component() < restrict_component_)
            //{
            //    throw new Exception($"Requested straight skeleton of component {restrict_component_} but max component index is {ct.max_component()}.");
            //};

            //FaceToTriangleIdxMap face_to_triangle_idx = new FaceToTriangleIdxMap();
            //TriangleOriginalVertexIndexList triangle_original_vertex_indices = new TriangleOriginalVertexIndexList();

            //// initialze the basic triangulation data structure, setting up neighborhoods
            //int num_initial_triangles = get_num_initial_triangles(ct);

            //DBG(//DBG_KT_SETUP) << "  Have " << num_initial_triangles << " initial kinetic triangles.";
            //   initialize_tds(mesh);

            //    create_supporting_lines(mesh);
            create_kinetic_vertices(mesh);
            /* until here, triangle_original_vertex_indices is consistent.
             * create_bevels is the first that may flip things around.
             */
            create_bevels(mesh);
            store_initial_wavefront_vertices();

            //    tidx_in_check_refinement.Capacity = triangles.size();
            refine_triangulation_initial();

            assert_valid(-1, CORE_ZERO);
            initialized = true;

            //DBG_FUNC_END(//DBG_KT_SETUP);
        }

        public partial void assert_valid()
        {
            foreach (var t in triangles)
            {
                if (t.is_dead()) continue;
                t.assert_valid();
            }
            assert(wavefront_edges != null);
            foreach (var e in wavefront_edges)
            {
                if (e.is_dead()) continue;
                assert(e.incident_triangle() != null);
                assert(e.incident_triangle().has_wavefront(e));
                assert(!e.incident_triangle().is_dead());
                // Remaining combinatorics is checked by triangle's assert_valid
            }
        }


       

        

        /** actually performs the flip, assert the triangulation is consistent before. */

        private partial void do_raw_flip(KineticTriangle t, int edge_idx, double time, bool allow_collinear)
        {
            assert(t != null);
            assert(!t.is_constrained(edge_idx));
            KineticTriangle n = t.neighbor(edge_idx);
            int nidx = n.index(t);

            WavefrontVertex v = t.vertex(edge_idx);
            WavefrontVertex v1 = t.vertex(ccw(edge_idx));
            WavefrontVertex v2 = t.vertex(cw(edge_idx));
            WavefrontVertex o = n.vertex(nidx);
            assert(v1 == n.vertex(cw(nidx)));
            assert(v2 == n.vertex(ccw(nidx)));

            /* We may also call this for two unbounded triangles.  However, right now we
             * only call this in one very specific way, so v1 is always the infinite and
             * v1 the finite one, and v, v1, and o are collinear on the CH boundary right
             * now.
             */
            assert(!v.is_infinite);
            assert(!v2.is_infinite);

            if (v1.is_infinite)
            {
                Point2 pos_v = v.p_at(time);
                Point2 pos_v2 = v2.p_at(time);
                Point2 pos_o = o.p_at(time);
                assert(orientation(pos_v, pos_o, pos_v2) == OrientationEnum.COLLINEAR);
            }
            else
            {
                Point2 pos_v = v.p_at(time);
                Point2 pos_v1 = v1.p_at(time);
                Point2 pos_v2 = v2.p_at(time);
                assert(orientation(pos_v1, pos_v2, pos_v) != OrientationEnum.RIGHT_TURN); // v may be on line(v1,v2)

                if (o.is_infinite)
                {
                    /* Flipping to an unbounded triangle. */
                    // nothing to do here.
                }
                else
                {
                    Point2 pos_o = o.p_at(time);

                    //DBG(//DBG_KT_EVENT2) << " o(v,o,v1):  " << SurfNet.orientation(pos_v, pos_o, pos_v1);
                    //DBG(//DBG_KT_EVENT2) << " o(v,o,v2):  " << SurfNet.orientation(pos_v, pos_o, pos_v2);
                    //DBG(//DBG_KT_EVENT2) << " o(v1,v2,o):  " << SurfNet.orientation(pos_v1, pos_v2, pos_o);

                    if (allow_collinear || true)
                    {
                        assert(orientation(pos_v, pos_o, pos_v1) != OrientationEnum.LEFT_TURN);
                        assert(orientation(pos_v, pos_o, pos_v2) != OrientationEnum.RIGHT_TURN);
                    }
                    else
                    {
                        assert(orientation(pos_v, pos_o, pos_v1) == OrientationEnum.RIGHT_TURN);
                        assert(orientation(pos_v, pos_o, pos_v2) == OrientationEnum.LEFT_TURN);
                    }
                    assert(orientation(pos_v1, pos_v2, pos_o) != OrientationEnum.LEFT_TURN); // The target triangle may be collinear even before.
                }
            }

            // not strictly necessary for flipping purpuses, but we probably
            // should not end up here if this doesn't hold:
            assert(!t.is_constrained(cw(edge_idx)) || !t.is_constrained(ccw(edge_idx)) || allow_collinear);
            t.do_raw_flip(edge_idx);
        }

        /** perform a flip, marking t and its neighbor as modified. */

        private partial void do_flip(KineticTriangle t, int idx, double time, bool allow_collinear)
        {
            assert(t != null);
            assert(!t.is_constrained(idx));
            KineticTriangle n = t.neighbor(idx);

          
            do_raw_flip(t, idx, time, allow_collinear);
          
            modified(t, true);
            modified(n);
        }

        /** process a flip evnt, checking all the assertions hold.  Calls do_flip(). */

        private partial void do_flip_event(double time, KineticTriangle t, int edge_idx)
        {
            //DBG_FUNC_BEGIN(//DBG_KT_EVENT);
            //DBG(//DBG_KT_EVENT) << &t << "; flip edge " << edge_idx;

            assert(0 <= edge_idx && edge_idx < 3);

            assert(t.wavefront(edge_idx) == null);
            assert(!t.is_constrained(edge_idx));

            WavefrontVertex[] v = new[]{ t.vertices[    edge_idx ],
                            t.vertices[ccw(edge_idx)],
                            t.vertices[cw (edge_idx)] };
            Point2[] p = new[] { v[0].p_at(time),
                         v[1].p_at(time),
                         v[2].p_at(time) };
            double[] squared_lengths = new[]{ squared_distance(p[1],p[2]),
                                  squared_distance(p[2],p[0]),
                                  squared_distance(p[0],p[1]) };
            assert(squared_lengths[0] > squared_lengths[1]);
            assert(squared_lengths[0] > squared_lengths[2]);
            assert(squared_lengths[1] != CORE_ZERO);
            assert(squared_lengths[2] != CORE_ZERO);

            do_flip(t, edge_idx, time);
            //DBG_FUNC_END(//DBG_KT_EVENT);
        }




        private partial void process_check_refinement_queue(double time)
        {
            //DBG_FUNC_BEGIN(//DBG_KT_EVENT);

            while (check_refinement.Count > 0)
            {
                KineticTriangle t = check_refinement_pop();
                refine_triangulation(t, time);
            }

            //DBG_FUNC_END(//DBG_KT_EVENT);
        }

        private partial KineticTriangle check_refinement_pop()
        {
            KineticTriangle t = check_refinement.First();
            check_refinement.RemoveAt(0);
            assert(tidx_in_check_refinement[t.Id]);
            tidx_in_check_refinement[t.Id] = false;
            return t;
        }

        private partial void put_on_check_refinement(KineticTriangle t, bool front)
        {
            assert(t != null);
            assert(tidx_in_check_refinement.Length > t.Id);

            if (tidx_in_check_refinement[t.Id]) return;

            tidx_in_check_refinement[t.Id] = true;
            if (front)
            {
                check_refinement.Insert(0, t);
            }
            else
            {
                check_refinement.Add(t);
            }
        }

        /** deal with spoke collapse
         *
         * The spoke at edge_idx collapsed, the incident vertices have been stopped
         * already.
         *
         * See if this triangle collapses entirely, or if its 3rd vertex (v) is elsewhere.
         *
         * If the triangle collapses entirely, we may have zero, one, or two constraints.
         * If it does not, we still may have zero or one constraint.
         * (Or two, if there is an infinitely fast vertex as v.)
         */

        private partial void do_spoke_collapse_part2(KineticTriangle t, int edge_idx, double time)
        {
            //DBG_FUNC_BEGIN(//DBG_KT_EVENT);
            //DBG(//DBG_KT_EVENT) << "t  " << &t;

            int num_constraints = t.is_constrained(0).ToInt() + t.is_constrained(1).ToInt() + t.is_constrained(2).ToInt();

            WavefrontVertex v = t.vertices[edge_idx];
            WavefrontVertex va = t.vertices[ccw(edge_idx)];
            WavefrontVertex vb = t.vertices[cw(edge_idx)];
            assert(va.has_stopped());
            assert(vb.has_stopped());
            assert(va.pos_stop() == vb.pos_stop());
            assert(!v.has_stopped());
            assert(t.neighbors[edge_idx] == null);

            var posa = va.pos_stop();

            if (!v.is_infinite && v.p_at(time).AreNear(posa))
            {
                t.set_dying();
                v.stop(time);
                /* triangle collapses completely */
                //DBG(//DBG_KT_EVENT) << "triangle collapses completely as spoke collapses; spoke is " << edge_idx;
                //LOG(WARNING) << __FILE__ << ":" << __LINE__ << " " << "untested code path.";

                int num_neighbors = 0;
                for (int i = 1; i <= 2; ++i)
                {
                    int edge = mod3(edge_idx + i);
                    //DBG(//DBG_KT_EVENT) << " at edge " << edge;
                    if (t.neighbors[edge] != null)
                    {
                        //DBG(//DBG_KT_EVENT) << "  we have a neighbor";
                        num_neighbors++;

                        KineticTriangle n = t.neighbors[edge];
                        int idx_in_n = n.index(t);

                        t.neighbors[edge] = null;
                        n.neighbors[idx_in_n] = null;

                        do_spoke_collapse_part2(n, idx_in_n, time);
                    }
                    else
                    {
                        //DBG(//DBG_KT_EVENT) << "  we have no neighbor";
                        //LOG(WARNING) << __FILE__ << ":" << __LINE__ << " " << "untested code path: DECL linking.";
                        t.vertices[ccw(edge)].set_next_vertex(1, t.vertices[cw(edge)], false);
                        if (t.wavefront(edge) != null)
                        {
                            t.wavefront(edge).set_dead();
                        }
                    }
                }
                queue.needs_dropping(t);
            }
            else
            {
                bool call_constraint_collapse = true;
                if (v.infinite_speed != InfiniteSpeedType.NONE)
                {
                    throw new NotImplementedException();
                    
                    int i_1 = mod3(edge_idx + 1);
                    int i_2 = mod3(edge_idx + 2);
                    if (t.neighbors[i_1] == null && t.neighbors[i_2] == null)
                    {
                        //DBG(//DBG_KT_EVENT) << "triangle is fully constraint and has an infinitely fast vertex.";

                        t.set_dying();
                        v.stop(time, posa);
                        for (int i = 1; i <= 2; ++i)
                        {
                            int edge = mod3(edge_idx + i);
                            //    LOG(WARNING) << __FILE__ << ":" << __LINE__ << " " << "untested code path: DECL linking.";
                            t.vertices[ccw(edge)].set_next_vertex(1, t.vertices[cw(edge)], false);
                            if (t.wavefront(edge) != null)
                            {
                                t.wavefront(edge).set_dead();
                            }
                        }
                        call_constraint_collapse = false;
                        queue.needs_dropping(t);
                    }
                }

                if (call_constraint_collapse)
                {
                    /* triangle does not collapse completely */
                    //DBG(//DBG_KT_EVENT) << "triangle does not collapse completely";
                    do_constraint_collapse_part2(t, edge_idx, time);
                }
            }

            //DBG_FUNC_END(//DBG_KT_EVENT);
        }

       
     

       

       

       

       

        // TODO
      

     
      

       
        

     

        public partial void update_event_timing_stats(double now)
        {
            if (now != this.last_event_time)
            {
                this.max_events_per_time = Math.Max(max_events_per_time, this.events_per_current_event_time);
                this.avg_events_per_time_sum += this.events_per_current_event_time;
                ++this.avg_events_per_time_ctr;

                this.last_event_time = now;
                this.events_per_current_event_time = 1;
            }
            else
            {
                ++this.events_per_current_event_time;
            }
        }

       
    }
}