namespace SurfNet
{
    using static DebugLog;
    using static Mathex;
   
    public partial class KineticTriangulation
    {
        TriangleOriginalVertexIndexList triangle_original_vertex_indices = new TriangleOriginalVertexIndexList();

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
            assert(vidx < input.Vertices.Count);
            return vidx;
        }

        /** Return a BasicVertex from a triangle index and vertex position within the triangle.
         *
         * See get_basic_vertex_idx_from_triangle_vertex_indices() for Debug.
         */

        private static partial BasicVertex get_basic_vertex_from_triangle_vertex_indices(
          BasicInput input,
          TriangleOriginalVertexIndexList triangle_original_vertex_indices,
          int t_idx,
          int i
        )
        {
            BasicVertex bv = input.Vertices[get_basic_vertex_idx_from_triangle_vertex_indices(input, triangle_original_vertex_indices, t_idx, i)];
            return bv;
        }

       

        /** Create wavefront vertices, except those that need beveling (such as degree-1 vertices)
         */

        
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

            //DBG(//DBG_KT_SETUP) << " Old triangle Debug: " << t;
            //DBG(//DBG_KT_SETUP) << " New triangle Debug: " << &new_t;

            //DBG_FUNC_END(//DBG_KT_SETUP);
            return new_t;
        }

       
       

       

       

        



        //public partial void initialize(
        //    DcelMesh mesh,
        //    WavefrontEdgeList p_wavefront_edges,
        //    int p_restrict_component)
        //{
        //    //DBG_FUNC_BEGIN(//DBG_KT_SETUP);

        //    //assert(!initialized);
        //    //wavefront_edges = p_wavefront_edges;
        //    //restrict_component_ = p_restrict_component;

        //    ///* set up basic triangulation */
        //    //BasicTriangulation ct = new BasicTriangulation();
        //    //ct.initialize(input);

        //    ////DBG(//DBG_KT_SETUP) << "Input has " << ct.max_component() + 1 << " component(s): 0 .. " << ct.max_component();
        //    ////DBG(//DBG_KT_SETUP) << "Requested straight skeleton of component " << restrict_component_;

        //    //if (ct.max_component() < restrict_component_)
        //    //{
        //    //    throw new Exception($"Requested straight skeleton of component {restrict_component_} but max component index is {ct.max_component()}.");
        //    //};

        //    //FaceToTriangleIdxMap face_to_triangle_idx = new FaceToTriangleIdxMap();
        //    //TriangleOriginalVertexIndexList triangle_original_vertex_indices = new TriangleOriginalVertexIndexList();

        //    //// initialze the basic triangulation data structure, setting up neighborhoods
        //    //int num_initial_triangles = get_num_initial_triangles(ct);

        //    //DBG(//DBG_KT_SETUP) << "  Have " << num_initial_triangles << " initial kinetic triangles.";
        //    //   initialize_tds(mesh);

        //    //    create_supporting_lines(mesh);
        //    create_kinetic_vertices(mesh);
        //    /* until here, triangle_original_vertex_indices is consistent.
        //     * create_bevels is the first that may flip things around.
        //     */
        //    create_bevels(mesh);
        //    store_initial_wavefront_vertices();

        //    //    tidx_in_check_refinement.Capacity = triangles.size();
        //    refine_triangulation_initial();

        //    assert_valid(-1, CORE_ZERO);
        //    initialized = true;

        //    //DBG_FUNC_END(//DBG_KT_SETUP);
        //}

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