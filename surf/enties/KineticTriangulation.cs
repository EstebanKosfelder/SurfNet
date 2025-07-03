
using System.Collections.Generic;
using System.Text.Json;
using TriangleNet.Geometry;
using TriangleNet.Topology.DCEL;

namespace SurfNet
{
    using static DebugLog;
    using static Mathex;
    using TriangleList = List<KineticTriangle>;
 


    public partial class KineticTriangulation
    {

    }

    public partial class KineticTriangulation
    {
        /* Iterate around the vertex identified by v_in_t_idx in t,
         * either clockwise (++) or counterclockwise (--).
         *
         * The neighborhood relation must be consistent across boundaries we traverse.
         */

        private class AroundVertexIterator //:
                                           //public std.iterator<std.bidirectional_iterator_tag
                                           //public std.iterator<std.forward_iterator_tag,
                                           //                     std.pair<KineticTriangle *,int>>
        {
            //  public:
            // using iterator_category = std.bidirectional_iterator_tag;
            //using value_type = std.pair<KineticTriangle *,int>;
            //using value_type = void;
            //using difference_type = void;
            //using pointer = void;
            //using reference = void;
            private KineticTriangle t_;

            private int v_in_t_idx_;




            private KineticTriangle next_triangle(int[] direction)
            {
                assert(v_in_t_idx_ < 3);
                return t_.neighbor(direction[v_in_t_idx_]);
            }

            private AroundVertexIterator walk_dir(int[] direction)
            {
                assert(v_in_t_idx_ < 3);
                KineticTriangle next = next_triangle(direction);
                int v_in_n_idx;
                if (next != null)
                {
                    int t_in_next_idx = next.index(t_);
                    assert(t_in_next_idx < 3);
                    v_in_n_idx = direction[t_in_next_idx];
                }
                else
                {
                    v_in_n_idx = 0;
                }
                t_ = next;
                v_in_t_idx_ = v_in_n_idx;
                return this;
            }

            public bool IsEnd => t_ == null && v_in_t_idx_ == 0;

            private AroundVertexIterator most_dir(int[] direction)
            {
                var  res = this;
                var t = res.t_;
                while (res.next_triangle(direction) != null)
                {
                    res.walk_dir(direction);
                    assert(res.t_ != t);
                   
                }
                return res;
            }

            public AroundVertexIterator(KineticTriangle t = null, int v_in_t_idx = 0)

            {
                t_ = t;
                v_in_t_idx_ = v_in_t_idx;
                assert(v_in_t_idx < 3);
            }

           public  AroundVertexIterator Copy()
            {
                return new AroundVertexIterator(t_, v_in_t_idx_);
            
            }

            public KineticTriangle next_triangle_cw()
            { return next_triangle(TriangulationUtils._cw); }

            public KineticTriangle next_triangle_ccw()
            { return next_triangle(TriangulationUtils._ccw); }

            public AroundVertexIterator most_cw()
            { return most_dir(TriangulationUtils._cw); }

            public AroundVertexIterator most_ccw()
            { return most_dir(TriangulationUtils._ccw); }

            public AroundVertexIterator next() { return walk_dir(TriangulationUtils._cw); }
            public AroundVertexIterator prev() { return walk_dir(TriangulationUtils._ccw); }

            public AroundVertexIterator walk_dir_cw => walk_dir(TriangulationUtils._cw);
            public AroundVertexIterator walk_dir_ccw => walk_dir(TriangulationUtils._ccw);
            //public static bool operator ==(AroundVertexIterator me, AroundVertexIterator other)
            //{
            //    return me.t_ == other.t_ && me.v_in_t_idx_ == other.v_in_t_idx_;
            //}

            //public static bool operator !=(AroundVertexIterator me, AroundVertexIterator other)
            //{
            //    return !(me == other);
            //}

            //public static implicit operator KineticTriangle(AroundVertexIterator me)
            //{ return me.t_; }

            public KineticTriangle t()
            { return t_; }

            public int v_in_t_idx()
            { return v_in_t_idx_; }

           
        };

        //friend std.ostream & operator <<(std.ostream& os, const KineticTriangulation.AroundVertexIterator it);

        private bool initialized = false;
        private bool finalized = false;
        public WavefrontEdgeList wavefront_edges = new WavefrontEdgeList();
        public WavefrontVertexList vertices = new WavefrontVertexList();
        public TriangleList triangles = new TriangleList();
        private int restrict_component_ = -1;

        public int restrict_component()
        { return restrict_component_; }

        /** The straight skeleton as a DCEL data structure.
         *
         * This is built in two parts.  Initially, we only create halfedges for all
         * input segments, and corresponding faces.  We do not created faces or
         * anything for beveled vertices.
         *
         * Once the propagation is done, we walk through all wavefront vertices (which
         * correspond to arcs in the straight skeleton and created the corresponding
         * halfedges and, where requred (beveling), also new faces.
         * Here we also create vertices.
         */
        private SkeletonDCEL skeleton=new SkeletonDCEL();

        public EventQueue queue;

        private partial (Point2, bool) get_vertex_pos(
              BasicInput input,
              TriangleOriginalVertexIndexList triangle_original_vertex_indices,
              KineticTriangle t,
              int i);

        private partial KineticTriangle split_vertex(
          KineticTriangle t,
          int i,
          WavefrontEdge new_e,
          WavefrontVertex new_v);

    

        private static partial int get_basic_vertex_idx_from_triangle_vertex_indices(
            BasicInput input,
            TriangleOriginalVertexIndexList triangle_original_vertex_indices,
            int t_idx,
            int i);

        private static partial BasicVertex get_basic_vertex_from_triangle_vertex_indices(
            BasicInput input,
            TriangleOriginalVertexIndexList triangle_original_vertex_indices,
            int t_idx,
            int i);

        private static partial void invalidate_basic_vertex_idx_in_triangle_vertex_indices(
            BasicInput input,
            TriangleOriginalVertexIndexList triangle_original_vertex_indices,
            int t_idx,
            int i);


        internal partial int get_num_initial_triangles(DcelMesh mesh);

        //private partial void initialize_tds(DcelMesh mesh);

        //internal partial void create_supporting_lines(DcelMesh mesh);

        

        private partial void store_initial_wavefront_vertices();



        ///<summary>
        /// build the remaining skeleton dcel
        ///
        /// For each "initial" vertex, build the structure of its right-hand face, going
        /// around it in a ccw direction.Its corresponding dcel face already exists if
        /// the face is traced out by an input edge.For edges that are the result of
        /// beveling, we need to create the dcel face (and the cbb) first.
        ///
        /// At each kinetic vertex in this walk, if it's not degenerate (i.e. if it has
        /// length > 0), create a halfedge pair unless we already have it.Then link up
        /// the correct halfedge with its predecessor in the list, and set the
        /// appropriate dcel halfedge in the kinetic vertex.
        ///
        /// may only be called once, at the end.
        ///</summary>

        public  void create_remaining_skeleton_dcel()
        {
            //   throw new NotImplementedException();
            //DBG_FUNC_BEGIN(//DBG_SKEL);
            assert(!finalized);

            /* set up halfedges */

            foreach (var v_it in vertices.Take(this.ContourVertices.Count))
            {
                WavefrontVertex v = v_it;
                bool new_face;
                SkeletonDCELFace f;
                SkeletonDCELCbb c;
                SkeletonDCELHalfedge @base;
                //DBG(//DBG_SKEL) << " v: " << v;

                if (v.is_infinite) continue;
                assert(v.is_initial);
                if (v.is_degenerate()) continue;

                assert(!v.is_degenerate());
                assert(v.incident_wavefront_edge(0) != null);
                if (null != (f = v.incident_wavefront_edge(0).skeleton_face))
                {
                    new_face = false;

                    assert(f.outer_ccbs_begin() != f.outer_ccbs_end());
                    @base = f.outer_ccbs_begin();
                    c = @base.outer_ccb();
                }
                else
                {
                    /* the result of beveling */
                    new_face = true;

                    f = skeleton.new_face();
                    c = skeleton.new_outer_ccb();
                    c.set_face(f);
                    @base = null;
                };
                //DBG(//DBG_SKEL) << "new_face: " << new_face;
                SkeletonDCELHalfedge last_he = create_remaining_skeleton_dcel_one_face(v, @base, c, f);
                if (new_face)
                {
                    assert(@base == null);
                    f.add_outer_ccb(c, last_he);
                }
                //DBG(//DBG_SKEL) << "skeleton_face: " << *f;
                if (@base == null)
                {
                    @base = f.outer_ccbs_begin();
                }
                //DBG(//DBG_SKEL) << "base         : " << *base;
                for (SkeletonDCELHalfedge he = @base.next(); he != @base; he = he.next())
                {
                    //DBG(//DBG_SKEL) << " he        : " << *he;
                }
            }
            foreach (var v_it in vertices.Skip(this.ContourVertices.Count))
            {
                assert(!v_it.is_initial);
            };

            if (restrict_component_ >= 0)
            {
                link_dcel_halfedges_on_ignored_side();
            }

            foreach (var v in vertices)
            {
                assert(v.is_infinite.ToInt() + v.is_degenerate().ToInt() == 1 || v.skeleton_dcel_halfedge(0) != null);
                assert(v.is_infinite.ToInt() + v.is_degenerate().ToInt() == 1 || v.skeleton_dcel_halfedge(1) != null);
            }

            /* set up vertices */
            bool have_infinite_vertex;
            if (restrict_component_ >= 0)
            {
                have_infinite_vertex = false;
                foreach (var v in vertices)
                {
                    if (v.is_infinite)
                    {
                        continue;
                    }
                    if (!v.has_stopped())
                    {
                        have_infinite_vertex = true;
                        break;
                    }
                }
                //DBG(//DBG_SKEL) << " have an infinite vertex: " << have_infinite_vertex;
                if (!have_infinite_vertex)
                {
                    skeleton.set_num_v_skew(0);
                };
            }
            else
            {
                have_infinite_vertex = true;
            }
            skeleton.set_number_of_points_and_curves();

            bool did_infinite_vertex = false;
            foreach (var v in vertices)
            {
                if (v.is_infinite || v.is_degenerate())
                {
                    continue;
                };
                //DBG(//DBG_SKEL) << "  Setting dcel vertex for " << v.details();
                SkeletonDCELHalfedge he;
                if (v.is_initial)
                {
                    he = v.skeleton_dcel_halfedge(1);
                    //DBG(//DBG_SKEL) << "  he (initial): " << *he;
                    if (he.vertex() == null)
                    {
                        assert(v.time_start == CORE_ZERO);
                        set_dcel_vertex(he, v.pos_start, CORE_ZERO);
                    }
                }

                he = v.skeleton_dcel_halfedge(0);
                //DBG(//DBG_SKEL) << "  he: " << *he;
                if (he.vertex() != null)
                {
                    if (v.has_stopped())
                    {
                        set_dcel_vertex(he, v.pos_stop(), v.time_stop());
                    }
                    else
                    {
                        assert(!did_infinite_vertex);
                        // DEBUG_STMT(did_infinite_vertex = true);
                        set_dcel_vertex(he, null, CORE_ZERO);
                    }
                }
            }
            assert(did_infinite_vertex ^ !have_infinite_vertex);

            /* set up segments/rays for arcs */
            foreach (var v in vertices)
            {
                if (v.is_infinite || v.is_degenerate())
                {
                    continue;
                };
                //DBG(//DBG_SKEL) << "  at v: " << &v;

                SkeletonDCELHalfedge he0 = v.skeleton_dcel_halfedge(0);
                SkeletonDCELHalfedge he1 = v.skeleton_dcel_halfedge(1);

                //DBG(//DBG_SKEL) << "   he0: " << *he0;
                //DBG(//DBG_SKEL) << "   he1: " << *he1;

                X_monotone_curve p;
                if (v.has_stopped())
                {
                    assert(!he0.vertex().has_null_point());
                    assert(!he1.vertex().has_null_point());
                    p = skeleton.new_segment(new Segment_3(he0.vertex().Point(), he1.vertex().Point()));
                }
                else
                {
                    assert(he0.vertex().has_null_point());
                    assert(!he1.vertex().has_null_point());
                    Vector_3 vec = new Vector_3(v.velocity.x(), v.velocity.y(), 0);
                    p = skeleton.new_ray(new Ray3(he1.vertex().point(), vec));
                }
                he0.set_curve(p);
                he1.set_curve(p);
            }
            /* and input segments */
            for (var he = skeleton.halfedges_begin(); he != skeleton.halfedges_end(); he = he.next())
            {
                if (he.is_emanating_input())
                {
                    if (he.has_null_curve())
                    {
                        X_monotone_curve p = skeleton.new_segment(new Segment_3(he.vertex().point(), he.opposite().vertex().point()));
                        he.set_curve(p);
                        he.opposite().set_curve(p);
                    }
                }
                else
                {
                    assert(!he.has_null_curve());
                }
            }
            skeleton.assert_sane();

            finalized = true;

        }

        private void link_dcel_halfedges_on_ignored_side()
        {
            
            //DBG_FUNC_BEGIN(//DBG_SKEL);

            foreach (var he in skeleton.halfedges)
            {
                if (!he.is_emanating_input() && he.opposite().is_emanating_input())
                {
                    /* This is the "ignored" side of an input edge. */
                    if (he.next() != null)
                    { // already processed
                        continue;
                    }

                    SkeletonDCELFace f;
                    SkeletonDCELCbb c;
                    f = skeleton.new_face();
                    c = skeleton.new_outer_ccb();
                    c.set_face(f);
                    var cur_he = he;

                    do
                    {
                        /* find next */
                        var next = cur_he.opposite();
                        while ( next.prev()!=null)
                        {
                            next = next.prev().opposite();
                            assert(next!=null);
                        }
                        cur_he.set_next(next);
                        cur_he = next;
                    } while (cur_he != he);

                    f.add_outer_ccb(c, he);
                }
            }
           
        }

        /** Create a dcel vertex and set it at all halfedges incident to that vertex.
       */

        private void set_dcel_vertex(SkeletonDCELHalfedge start, Point2? p, double time)
        {

            SkeletonDCELVertex new_v = skeleton.new_vertex();
            new_v.set_halfedge(start);
            if (p != null)
            {
                Point_3 pp = skeleton.new_point(new Point_3(p.Value.X, p.Value.Y, time));
                new_v.set_point(pp);
            }

            int degree = 0;
            SkeletonDCELHalfedge he = start;
            do
            {
                he.set_vertex(new_v);
                //DBG(//DBG_SKEL) << "  Setting vertex for " << *he;
                assert(he.next() != null);
                if (he.next() == null)
                {
                    /* If we miss events, then we can end up in an inconsistent state.
                     * See for instance https://github.com/cgalab/surfer2/issues/2
                     * This guards against the segfault and causes an abort instead.  Not
                     * convinced that's the way to go, but still.
                     */
                    throw new Exception("We were asked to create a DCEL vertex but the half-edge has a null next pointer.  This should not happen.");

                }
                he = he.next().opposite();
                degree++;
            } while (he != start);
        }

        /** Set up DCEL half-edges around one straight skeleton face
       *
       * The kinetic vertex start is the "right" vertex of base, which traces out f.  We are walking
       * along the boundary of f counter-clockwise.
       */

        private SkeletonDCELHalfedge create_remaining_skeleton_dcel_one_face(WavefrontVertex start, SkeletonDCELHalfedge heBase, SkeletonDCELCbb c, SkeletonDCELFace f)
        {


            //DBG_FUNC_BEGIN(//DBG_SKEL);
            assert(start.is_initial);
            assert(c!=null);
            assert(f!=null);

            if (heBase!=null)
            {
                //DBG(//DBG_SKEL) << "starting at : " << start << "; base: " << *base;
            }
            else
            {
                //DBG(//DBG_SKEL) << "starting at : " << start << "; base: -";
            }
            SkeletonDCELHalfedge prev_he = heBase;
            WavefrontVertex cur_v = start;
            WavefrontVertex prev_v;
            int side_of_f = 0;
            bool looped_around = false;
            while (!looped_around )
            {
                //DBG(//DBG_SKEL) << " doing " << cur_v << "; " << cur_v.details();
                assert(!cur_v.is_degenerate());

                SkeletonDCELHalfedge new_he;
                /* Make a new halfedge or use the one we already have for this vertex. */
                assert(cur_v.skeleton_dcel_halfedge(side_of_f)==null);
                if (cur_v.skeleton_dcel_halfedge(1 - side_of_f)!=null)
                {
                    new_he = cur_v.skeleton_dcel_halfedge(1 - side_of_f).opposite();
                    assert(new_he.is_on_outer_ccb());
                    assert(new_he.outer_ccb() == null);
                }
                else
                {
                    new_he = skeleton.new_edge();
                }
                new_he.set_outer_ccb(c);
                //DBG(//DBG_SKEL) << "  he for " << cur_v << ": " << *new_he;

                if (prev_he!=null)
                {
                    assert(prev_he.next() == null);
                    new_he.set_prev(prev_he);
                }
                prev_he = new_he;
                cur_v.set_skeleton_dcel_halfedge(side_of_f, new_he);

                //DBG(//DBG_SKEL) << "  looking for next. ";
                bool foundvertex = false;
                do
                {
                    prev_v = cur_v;
                    /* if side_of_f is 1 (i.e. if the face is to the right of us, we are
                     * going forward on this vertex, else we are going backwards. */
                    if (side_of_f!=0)
                    {
                        cur_v = cur_v.prev_vertex(side_of_f);
                    }
                    else
                    {
                        cur_v = cur_v.next_vertex(side_of_f);
                    }
                    //DBG(//DBG_SKEL) << "  considering " << cur_v;

                    if (cur_v == null && side_of_f == 0)
                    { /* prev wavefront vertex escapes to infinity */
                        //DBG(//DBG_SKEL) << "  went to infinity (and beyond?).  jumping over using wavefront edge";
                        assert(!prev_v.has_stopped());
                        WavefrontEdge e = prev_v.incident_wavefront_edge(0);
                        assert(!e.is_dead());
                        assert(e.vertex(1) == prev_v);
                        cur_v = e.vertex(0);
                        side_of_f = 1;
                        foundvertex = true;
                        
                    }
                    else if (cur_v == null && side_of_f == 1)
                    {
                        //DBG(//DBG_SKEL) << "  done now, arrived at base.";
                        assert(prev_v.is_initial);
                       looped_around=true;
                    }
                    else if (cur_v == start)
                    {
                        //DBG(//DBG_SKEL) << "  done now, looped around";
                        assert(cur_v.prev_vertex(0) == prev_v);
                        side_of_f = 0;
                        assert(prev_v.is_initial);
                        assert(cur_v.is_beveling || prev_v.is_beveling);
                        looped_around = true;
                    }
                    else
                    {
                        assert(cur_v!=null);

                        /* Figure out which side f is one */
                        if (cur_v.prev_vertex(0) == prev_v)
                        {
                            side_of_f = 0;
                        }
                        else
                        {
                            assert(cur_v.next_vertex(1) == prev_v);
                            side_of_f = 1;
                        }
                    }
                } while ( cur_v.is_degenerate() && !foundvertex && !looped_around);
                if (foundvertex)
                {
                    //DBG(//DBG_SKEL) << "  found " << cur_v;
                }
                foundvertex = false;
            };

       
            if (cur_v!=null)
            {
                //DBG(//DBG_SKEL) << "back at start, so this is the result of beveling";
                /* beveling vertex, this is the start vertex, again, hopefully. */
                assert(cur_v == start);
                assert(side_of_f == 0);
                assert(cur_v.skeleton_dcel_halfedge(side_of_f)!=null);
                cur_v.skeleton_dcel_halfedge(side_of_f).set_prev(prev_he);
            }
            else
            {
                assert(heBase.prev() == null);
                heBase.set_prev(prev_he);
            }
            //DBG_FUNC_END(//DBG_SKEL);
            return prev_he;
        }

        private AroundVertexIterator incident_faces_iterator(KineticTriangle t, int v_in_t)
        {
            return new AroundVertexIterator(t, v_in_t);
        }

      

        public KineticTriangulation() { }

       

        public partial void assert_valid();


        public void assert_valid(int current_component, double time)
        {
#if !DEBUG_EXPENSIVE_PREDICATES
            //DBG_FUNC_BEGIN(//DBG_KT_EVENT2);
            assert_valid();
            foreach (var t in triangles)
            {
                if (restrict_component_ >= 0)
                {
                    /* We only do SK in one component */
                    if (t.component != restrict_component_) continue;
                }
                else
                {
                    /* We do SK in all components, but the time in the others will be earlier/later so triangles will not be right */
                    if (t.component != current_component) continue;
                }
                if (t.is_dead()) continue;
                if (t.unbounded())
                {
                    // recall that unbounded triangles witness that the vertex ccw of the infinite
                    // vertex is on the boundary of the CH.
                    //DBG(//DBG_KT_EVENT2) << "t" << &t;
                    int idx = t.infinite_vertex_idx();
                    KineticTriangle n = t.neighbor(cw(idx));
                    int nidx = n.infinite_vertex_idx();

                    WavefrontVertex u = t.vertex(cw(idx));
                    WavefrontVertex v = t.vertex(ccw(idx));
                    WavefrontVertex V = n.vertex(cw(nidx));
                    WavefrontVertex w = n.vertex(ccw(nidx));
                    assert(v == V);
                    assert(orientation(u.p_at(time),
                                             v.p_at(time),
                                             w.p_at(time)) != OrientationEnum.RIGHT_TURN);
                }
                else
                {
                    //DBG(//DBG_KT_EVENT2) << "t" << &t;
                    Point2 a = (t.vertex(0).p_at(time));
                    Point2 b = (t.vertex(1).p_at(time));
                    Point2 c = (t.vertex(2).p_at(time));
                    double det = compute_determinant(
                        a.X, a.Y,
                        b.X, b.Y,
                        c.X, c.Y);
                    if (det.IsZero()) det = 0;
                    if ((orientation(a, b, c) == OrientationEnum.RIGHT_TURN) != (det < 0))
                    {
                        assert($"CGAL is confused about orientation of triangle {t}: determinant is {det} but CGAL thinks orientation is {orientation(a, b, c)}");

                    }
                    assert(orientation(a, b, c) != OrientationEnum.RIGHT_TURN, $"kt:{t.Id} time:{time} o:{(int)orientation(a, b, c)}");
                }
            }
            foreach (var v in vertices)
            {
                v.assert_valid();
            }
            //DBG_FUNC_END(//DBG_KT_EVENT2);
#else
            assert_valid();
            foreach (var v in vertices)
            {
                v.assert_valid();
            }
#endif
        }

        

        public void set_queue(EventQueue q)
        { queue = q; }

        //        //const FixedVector<WavefrontVertex>& get_vertices() { return vertices; };
        //        public const VertexList.const_iterator vertices_begin() { return vertices.begin(); }

        //    public const VertexList.const_iterator vertices_end() { return vertices.end(); }

        //public const TriangleList.const_iterator triangles_begin() { return triangles.begin(); }
        //public const TriangleList.const_iterator triangles_end() { return triangles.end(); }
        public int triangles_size() { return triangles.size(); }

        /* profiling/statistics counters */
        public int[] event_type_counter = new int[(int)(CollapseType.NEVER)];

        public int max_triangles_per_edge_event = 0;
        public int avg_triangles_per_edge_event_sum = 0;
        public int avg_triangles_per_edge_event_ctr = 0;

        public int max_triangles_per_split_event = 0;
        public int avg_triangles_per_split_event_sum = 0;
        public int avg_triangles_per_split_event_ctr = 0;

        public double last_event_time = 0;
        public int events_per_current_event_time = 0;
        public int max_events_per_time = 0;
        public int avg_events_per_time_sum = 0;
        public int avg_events_per_time_ctr = 0;
        public partial void update_event_timing_stats(double now);

        private partial void do_raw_flip(KineticTriangle t, int edge_idx, double time, bool allow_collinear);
        private partial void do_flip(KineticTriangle t, int edge_idx, double time, bool allow_collinear = false);
        private partial void do_flip_event(double time, KineticTriangle t, int edge_idx);

        private partial void do_spoke_collapse_part2(KineticTriangle t, int edge_idx, double time);

        /** Handle the 2nd part of a collapse of one constraint from a constraint evnt.
        *
        * After the old kinetic vertices are stopped, this creates a new kv, and updates
        * all incident triangles.
        *
        * We may also call this when a spoke collapses, in which case t.wavefront(edge_idx)
        * is null, but it will also not have a neighbor there anymore.
        */
        private void do_constraint_collapse_part2(KineticTriangle t, int edge_idx, double time)
        {



            // throw new NotImplementedException();
            //DBG_FUNC_BEGIN(//DBG_KT_EVENT);

            WavefrontVertex va = t.vertices[ccw(edge_idx)];
            WavefrontVertex vb = t.vertices[cw(edge_idx)];
            assert(va.has_stopped());
            assert(vb.has_stopped());
            assert(va.pos_stop().AreNear(vb.pos_stop()));

            var pos = va.pos_stop();

            t.set_dying();

            move_constraints_to_neighbor(t, edge_idx);

            if (t.wavefront(edge_idx) != null)
            {
                /* constraint collapsed */
                assert(t.wavefront(edge_idx) == va.incident_wavefront_edge(1));
                assert(t.wavefront(edge_idx) == vb.incident_wavefront_edge(0));
            }
            else
            {
                /* spoke collapsed */
                assert(t.neighbor(edge_idx) == null);
            }
            WavefrontEdge ea = va.incident_wavefront_edge(0);
            WavefrontEdge eb = vb.incident_wavefront_edge(1);

            KineticTriangle? na = t.neighbors[cw(edge_idx)];
            KineticTriangle? nb = t.neighbors[ccw(edge_idx)];

            WavefrontVertex v = vertices.make_vertex(pos, time, ea, eb);
            //DBG(//DBG_KT) << " va: " << va;
            //DBG(//DBG_KT) << " vb: " << vb;
            va.set_next_vertex(0, v);
            vb.set_next_vertex(1, v);

            {
                int affected_triangles = 0;

                //DBG(//DBG_KT_EVENT) << "updating vertex in affected triangles";
              
                AroundVertexIterator i = incident_faces_iterator(t, ccw(edge_idx));
                bool first = true;
                KineticTriangle ti;
                //DBG(//DBG_KT_EVENT) << " ccw:";
                for (i.prev(); !i.IsEnd; i.prev())
                {
                    ti = i.t();

                    assert(!first || na == ti);
                    first = false;
                    ti.set_vertex(i.v_in_t_idx(), v);
                  
                    modified(ti);
                    ++affected_triangles;
                };

                i = incident_faces_iterator(t, cw(edge_idx));
                first = true;
                //DBG(//DBG_KT_EVENT) << " cw:";
                for (i.next(); !i.IsEnd; i.next())
                {
                    ti = i.t();
                    assert(!first || nb == ti);
                    first = false;
                    ti.set_vertex(i.v_in_t_idx(), v);
                  
                    modified(ti);
                    ++affected_triangles;
                }

                max_triangles_per_edge_event = Math.Max(max_triangles_per_edge_event, affected_triangles);
                avg_triangles_per_edge_event_sum += affected_triangles;
                ++avg_triangles_per_edge_event_ctr;
            };

            //DBG(//DBG_KT_EVENT) << " cw:";
            if (na != null)
            {
                var idx = na.index(t);
                assert(na.index(v) == ccw(idx));
                na.neighbors[idx] = nb;
            }
            if (nb != null)
            {
                var idx = nb.index(t);
                assert(nb.index(v) == cw(idx));
                nb.neighbors[idx] = na;
            }

            if (t.wavefront(edge_idx) != null)
            {
                /* if it was a constraint collapse rather than a spoke collapse */
                t.wavefront(edge_idx).set_dead();
            }
            queue.needs_dropping(t);


        }

        private  void handle_constraint_event(CollapseEvent evnt)
        {

            assert(evnt.type() == CollapseType.CONSTRAINT_COLLAPSE);
            KineticTriangle t = triangles[evnt.t.Id];
            double time = evnt.time();
            int edge_idx = evnt.relevant_edge();


          

            assert(t.wavefront(edge_idx) != null);
            assert(t.wavefront(edge_idx).get_collapse(t.component, time, edge_idx) == evnt);

            assert(t.is_constrained(edge_idx));

            WavefrontVertex va = t.vertices[ccw(edge_idx)];
            WavefrontVertex vb = t.vertices[cw(edge_idx)];
            
            va.stop(time);
            vb.stop(time);
            assert(va.p_at(time).AreNear(vb.p_at(time)) ,$" va:{va.p_at(time)} vb:{vb.p_at(time)}" );


            // update prev/next for the DCEL that is the wavefront vertices
            va.set_next_vertex(1, vb, false);

            do_constraint_collapse_part2(t, edge_idx, time);
          
            assert_valid(t.component, time);

            
        }


        private  void handle_spoke_collapse_event(CollapseEvent evnt)
        {
            //DBG_FUNC_BEGIN(//DBG_KT_EVENT);
            //DBG(//DBG_KT_EVENT) << evnt;

            assert(evnt.type() == CollapseType.SPOKE_COLLAPSE);
            KineticTriangle t = triangles[evnt.t.Id];
            double time = evnt.time();
            int edge_idx = evnt.relevant_edge();

            assert(!t.is_constrained(edge_idx));

            WavefrontVertex va = t.vertices[ccw(edge_idx)];
            WavefrontVertex vb = t.vertices[cw(edge_idx)];
            va.stop(time);
            vb.stop(time);
            assert(va.pos_stop() == vb.pos_stop());

            KineticTriangle n = t.neighbors[edge_idx];
            assert(n != null);
            int idx_in_n = n.index(t);

            t.neighbors[edge_idx] = null;
            n.neighbors[idx_in_n] = null;

            do_spoke_collapse_part2(t, edge_idx, time);
            do_spoke_collapse_part2(n, idx_in_n, time);

            // update prev/next for the DCEL that is the wavefront vertices
            /* actually, nothing to do here, the two do_spoke_collapse_part2 calls did everything. */
            //LOG(WARNING) << __FILE__ << ":" << __LINE__ << " " << "untested code path: DECL linking.";

            assert_valid(t.component, time);

            //DBG_FUNC_END(//DBG_KT_EVENT);
        }

        private void handle_triangle_collapse_event(CollapseEvent evnt)
        {
            //DBG_FUNC_BEGIN(//DBG_KT_EVENT);
            //DBG(//DBG_KT_EVENT) << evnt;

            assert(evnt.type() == CollapseType.TRIANGLE_COLLAPSE);
            KineticTriangle t = triangles[evnt.t.Id];
            double time = evnt.time();

            int num_constraints = t.is_constrained(0).ToInt() + t.is_constrained(1).ToInt() + t.is_constrained(2).ToInt();
            //DBG(//DBG_KT_EVENT) << "have " << num_constraints << " constraints";

          

            for (int i = 0; i < 3; ++i)
            {
                t.vertices[i].stop(time);
                //DBG(//DBG_KT_EVENT) << "v[" << i << "]: " << t.vertices[i];
            }
            for (int i = 1; i < 3; ++i)
            {
                assert(t.vertices[0].p_at(time).AreNear(t.vertices[i].p_at(time)));
                assert(t.vertices[0].pos_stop().AreNear(t.vertices[i].pos_stop()));
            }

            t.set_dying();
            for (int i = 0; i < 3; ++i)
            {
                if (t.is_constrained(i))
                {
                    assert(t.wavefront(i) == t.vertex(ccw(i)).incident_wavefront_edge(1));
                    assert(t.wavefront(i) == t.vertex(cw(i)).incident_wavefront_edge(0));

                    t.wavefront(i).set_dead();

                    // update prev/next for the DCEL that is the wavefront vertices
                    t.vertices[cw(i)].set_next_vertex(0, t.vertices[ccw(i)], false);
                }
                else
                {
                    // from the other triangle's point of view, a spoke collapsed.  deal with that.
                    KineticTriangle n = t.neighbors[i];
                    assert(n != null);
                    int idx_in_n = n.index(t);

                    t.neighbors[i] = null;
                    n.neighbors[idx_in_n] = null;

                    do_spoke_collapse_part2(n, idx_in_n, time);
                }
            }
            queue.needs_dropping(t);

            assert_valid(t.component, time);

            //DBG_FUNC_END(//DBG_KT_EVENT);
        }

        private void handle_split_event(CollapseEvent evnt)
        {

            Log($"Handle Split Event ");
            LogIndent();
            assert(evnt.type() == CollapseType.SPLIT_OR_FLIP_REFINE);
           
            KineticTriangle  t = triangles[evnt.t.Id];
            
            double time= evnt.time();
            int edge_idx = evnt.relevant_edge();

          

            //DBG(//DBG_KT_EVENT2) << " t:  " << &t;
            WavefrontVertex v = t.vertices[edge_idx];

            WavefrontEdge e = t.wavefront(edge_idx);
            WavefrontEdge eb = v.incident_wavefront_edge(0);
            WavefrontEdge ea = v.incident_wavefront_edge(1);
            assert(ea.vertex(0) == v);
            assert(eb.vertex(1) == v);

            KineticTriangle na = t.neighbors[cw(edge_idx)];
            KineticTriangle nb = t.neighbors[ccw(edge_idx)];

            { 
               // So, we have a wavefront edge e and we have an opposite reflex vertex v.
               // Assert that the vertex can actually hit the edge, i.e., check if the
               // wavefront edges e0 and e1 incident at v point away from e when
               // starting at v.
               //
               // The edges we store with v are directed such that e0 points towards v
               // and e1 away from it.
               //
               // The orientation needed for an evnt thus is that e and e0 form a right
               // turn, and e and e1 form a left turn.
               //
               // If either are collinear, things can still happen.
              

                Vector2  e0 =(t.vertices[edge_idx].incident_wavefront_edge(0).l().l.to_vector());
                Vector2  e1 =(t.vertices[edge_idx].incident_wavefront_edge(1).l().l.to_vector());
                var o0 = (orientation(e.l().l.to_vector(), e0));
                var o1 = (orientation(e.l().l.to_vector(), e1));

                assert(o0 != OrientationEnum.LEFT_TURN);
                assert(o1 != OrientationEnum.RIGHT_TURN);
            }
            v.stop(time);
            Point2 pos = v.pos_stop();
            // Stop the vertex
          

            // Split edge e into parts,
            var new_edges = e.split(wavefront_edges);
            WavefrontEdge nea = new_edges.Left;
            WavefrontEdge neb = new_edges.Right;

            // Create new wavefront vertices with the new edges
            WavefrontVertex nva = vertices.make_vertex(pos, time, nea, ea, true);
            WavefrontVertex nvb = vertices.make_vertex(pos, time, eb, neb, true);

            // And set these new vertices on the edges.
            nea.vertex(0).set_incident_wavefront_edge(1, nea);
            nea.set_wavefrontedge_vertex(1, nva);
            neb.set_wavefrontedge_vertex(0, nvb);
            neb.vertex(1).set_incident_wavefront_edge(0, neb);

            // And update prev/next for the DCEL that is the wavefront vertices
            v.set_next_vertex(0, nvb);
            v.set_next_vertex(1, nva);
            nva.link_tail_to_tail(nvb);

            t.set_dying();
            {
                int affected_triangles = 0;

                AroundVertexIterator i = incident_faces_iterator(t, edge_idx);
                KineticTriangle lasta = null;
                //DBG(//DBG_KT_EVENT2) << " split: updating vertex on a side:";
                for (i.next(); !i.IsEnd; i.next())
                {
                    var ti = i.t();
                   
                    //DBG(//DBG_KT_EVENT2) << " split:   updating vertex on a side in " << &*i;
                    ti.set_vertex(i.v_in_t_idx(), nva);
                    modified(ti);
                    lasta = ti;
                    ++affected_triangles;
                };
                assert(lasta != null);
                assert(lasta.wavefront(cw(lasta.index(nva)))!=null);
                assert(lasta.wavefront(cw(lasta.index(nva))).vertex(0) == nva);

                i = incident_faces_iterator(t, edge_idx);
                KineticTriangle? lastb = null;
                //DBG(//DBG_KT_EVENT2) << " split: updating vertex on b side: ";
                for (i.prev(); !i.IsEnd; i.prev())
                {
                    //DBG(//DBG_KT_EVENT2) << " split:   updating vertex on b side in " << &*i;
                    var ti = i.t();
                   
                    ti.set_vertex(i.v_in_t_idx(), nvb);
                    modified(ti);
                    lastb = ti;
                    ++affected_triangles;
                }

                max_triangles_per_split_event = Math.Max(max_triangles_per_split_event, affected_triangles);
                avg_triangles_per_split_event_sum += affected_triangles;
                ++avg_triangles_per_split_event_ctr;

                assert(lastb != null);
                assert(lastb.wavefront(ccw(lastb.index(nvb))) != null);
                assert(lastb.wavefront(ccw(lastb.index(nvb))).vertex(1) == nvb);

                assert(na.index(nva) == cw(na.index(t)));
                na.set_wavefront(na.index(t), nea);

                assert(nb.index(nvb) == ccw(nb.index(t)));
                nb.set_wavefront(nb.index(t), neb);

                //DBG(//DBG_KT_EVENT2) << " nea:" << *nea;
                //DBG(//DBG_KT_EVENT2) << " neb:" << *neb;
                //DBG(//DBG_KT_EVENT2) << " ea: " << *ea;
                //DBG(//DBG_KT_EVENT2) << " eb: " << *eb;
                //DBG(//DBG_KT_EVENT2) << " na: " << na;
                //DBG(//DBG_KT_EVENT2) << " nb: " << nb;

                na.assert_valid();
                nb.assert_valid();
                lasta.assert_valid();
                lastb.assert_valid();
              
            }

            queue.needs_dropping(t);
            assert_valid(t.component, time);

            LogIndent();
        }



        /// <summary>
        /// handle a split evnt, or refine this as a flip evnt.
        ///
        /// A triangle with exactly one constraint, e, has the opposite vertex v moving
        /// onto the supporting line of e.
        ///
        /// This can be a split evnt, in which case we handle it right here,
        /// or it can be a flip evnt, in which case we pump it down the line
        /// so we can first deal with other, real split events or higher-priority flip
        /// events.
        /// </summary>
        private void handle_split_or_flip_refine_event(CollapseEvent evnt)
        {
            //DBG_FUNC_BEGIN(//DBG_KT_EVENT);
            //DBG(//DBG_KT_EVENT) << evnt;

            assert(evnt.type() == CollapseType.SPLIT_OR_FLIP_REFINE);
            KineticTriangle t = triangles[evnt.t.Id];
            double time = (evnt.time());
            int edge_idx = evnt.relevant_edge();

            //DBG(//DBG_KT_EVENT2) << " t:  " << &t;
            assert(t.wavefront(edge_idx) != null);
            assert(t.is_constrained(edge_idx));
            assert(!t.is_constrained(cw(edge_idx)));
            assert(!t.is_constrained(ccw(edge_idx)));

            WavefrontVertex v = t.vertices[edge_idx];
            WavefrontVertex va = t.vertices[ccw(edge_idx)];
            WavefrontVertex vb = t.vertices[cw(edge_idx)];
            assert(!v.has_stopped());

            var posa = va.p_at(time);
            var posb = vb.p_at(time);
            var pos = v.p_at(time);
            Segment2 s = new Segment2(posa, posb);
            assert(s.supporting_line().has_on(pos));

            /* 3 cases:
             * (1) the vertex v is outside of the segment s, on the supporting line,
             * (2) it's at one of the endpoints of s,
             * (1) or v is on the interior of the edge.
             */

            double sq_length_constraint = s.squared_length();
            double sq_length_v_to_va = squared_distance(pos, posa);
            double sq_length_v_to_vb = squared_distance(pos, posb);
            double longest_spoke;

            // //DBG(//DBG_KT_EVENT) << " v.pos: " << CGAL_point(pos);
            // //DBG(//DBG_KT_EVENT) << " va.pos: " << CGAL_point(posa);
            // //DBG(//DBG_KT_EVENT) << " vb.pos: " << CGAL_point(posb);
            // //DBG(//DBG_KT_EVENT) << " sqlength s: " << SurfNet.to_double(sq_length_constraint);
            // //DBG(//DBG_KT_EVENT) << " sqlength v-va: " << SurfNet.to_double(sq_length_v_to_va);
            // //DBG(//DBG_KT_EVENT) << " sqlength v-vb: " << SurfNet.to_double(sq_length_v_to_vb);
            if (!(s.CollinearHasOn(pos)|| pos.AreNear(posa)|| pos.AreNear(posa)))
            {
                // case 1
                // VLOG(2) << "A potential split evnt is actually a flip evnt.  Maybe refinement should have prevented that?";
                //DBG(//DBG_KT_EVENT) << "Re-classifying as flip evnt as v is not on the constrained segment.";

                assert(KineticTriangle.edge_is_faster_than_vertex(t.vertex(edge_idx), t.wavefront(edge_idx).l()) != Mathex.NEGATIVE);

                /** there are basically two types of flip events involving constrained triangles.
                 * One is where the vertex is coming towards the supporting line of the constraint
                 * edge e and is passing left or right of of e.  The other is where (the
                 * supporting line of an edge e overtakes a vertex.
                 *
                 * They are generally handled identically, however the assertions are slightly different
                 * ones as these cases differ in which of v's incident edges is relevant.
                 */
                Vector2 e = (t.wavefront(edge_idx).l().l.to_vector());
                Vector2 e0 = (t.vertex(edge_idx).incident_wavefront_edge(0).l().l.to_vector());
                Vector2 e1 = (t.vertex(edge_idx).incident_wavefront_edge(1).l().l.to_vector());
                var o0 = (orientation(e, e0));
                var o1 = (orientation(e, e1));
                int flip_edge;

                /* Figure out which side of the constraint edge we're on. */
                if (sq_length_v_to_va > sq_length_v_to_vb)
                {
                    //DBG(//DBG_KT_EVENT) << "(v, va) is the longest spoke, so vb moves over that.";
                    longest_spoke = sq_length_v_to_va;

                    assert(sq_length_v_to_va > sq_length_constraint); // Check that v to va is the longest spoke //
                    assert(t.index(vb) == cw(edge_idx));
                    assert(vb.is_reflex_or_straight(),$"vb:{vb.Id}");

                    /* If we come from behind, we don't really care about the first of these things in the discunjtion,
                     * if we face it head on, we don't care about the second.  hmm. */
                    assert(o1 != OrientationEnum.RIGHT_TURN || (v.is_reflex_or_straight() && o0 != OrientationEnum.RIGHT_TURN));

                    flip_edge = cw(edge_idx);
                }
                else /* (pos, posa) < (pos, posb) */
                {
                    //DBG(//DBG_KT_EVENT) << "(v, vb) is the longest spoke, so va moves over that.";
                    longest_spoke = sq_length_v_to_vb;
                    assert(sq_length_v_to_va != sq_length_v_to_vb); /* they really shouldn't be able to be equal. */

                    assert(sq_length_v_to_vb > sq_length_constraint); /* Check that v to vb is the longest spoke */
                    assert(t.index(va) == ccw(edge_idx));
                    assert(va.is_reflex_or_straight());

                    /* If we come from behind, we don't really care about the first of these things in the discunjtion,
                     * if we face it head on, we don't care about the second.  hmm. */
                    assert(o0 != OrientationEnum.LEFT_TURN || (v.is_reflex_or_straight() && o1 != OrientationEnum.LEFT_TURN));

                    flip_edge = ccw(edge_idx);
                };
                CollapseSpec c = t.refine_collapse_spec(new CollapseSpec(t.component, CollapseType.VERTEX_MOVES_OVER_SPOKE, time, flip_edge, longest_spoke));
                //DBG(//DBG_KT_EVENT) << " Refining to " << c;
                queue.needs_update(t, true);
            }
            else if (pos.AreNear(posa))
            {
                // case 2
                CollapseSpec c = t.refine_collapse_spec(new CollapseSpec(t.component, CollapseType.SPOKE_COLLAPSE, time, cw(edge_idx)));
                //DBG(//DBG_KT_EVENT) << " v is incident to va Refining to " << c;
                queue.needs_update(t, true);
            }
            else if (pos.AreNear(posb))
            {
                // case 2
                CollapseSpec c = t.refine_collapse_spec(new CollapseSpec(t.component, CollapseType.SPOKE_COLLAPSE, time, ccw(edge_idx)));
                //DBG(//DBG_KT_EVENT) << " v is incident to vb Refining to " << c;
                queue.needs_update(t, true);
            }
            else
            {
                //DBG(//DBG_KT_EVENT) << "We have a real split evnt.";
                handle_split_event(evnt);
            }

           
        }

        private void handle_vertex_moves_over_spoke_event(CollapseEvent evnt)
        {
            

            assert(evnt.type() == CollapseType.VERTEX_MOVES_OVER_SPOKE);
            KineticTriangle t = triangles[evnt.t.Id];
         
            do_flip_event(evnt.time(), t, evnt.relevant_edge());
           

            
        }
        private  void handle_ccw_vertex_leaves_ch_event(CollapseEvent evnt)
        {
            throw new NotImplementedException();
            ////DBG_FUNC_BEGIN(//DBG_KT_EVENT);
            ////DBG(//DBG_KT_EVENT) << evnt;

            //assert(evnt.type() == CollapseType.CCW_VERTEX_LEAVES_CH);
            //KineticTriangle t = triangles[evnt.t.Id];
            //int idx = evnt.relevant_edge(); /* finite edge idx == infinite vertex idx */

            //assert((int)idx == t.infinite_vertex_idx());
            //do_flip(t, cw(idx), evnt.time());

            ////DBG_FUNC_END(//DBG_KT_EVENT);
        }

        private  void handle_face_with_infintely_fast_weighted_vertex(CollapseEvent evnt)
        {
            throw new NotImplementedException();
            ////DBG_FUNC_BEGIN(//DBG_KT_EVENT);
            ////DBG(//DBG_KT_EVENT) << evnt;

            //assert(evnt.type() == CollapseType.FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED);
            //KineticTriangle tref = triangles[evnt.t.id];
            //KineticTriangle t = tref;
            //double time(evnt.time());
            //int edge_idx = evnt.relevant_edge();

            ////DBG(//DBG_KT_EVENT2) << "t:  " << t;

            //assert((t.vertex(0).infinite_speed == InfiniteSpeedType.WEIGHTED) ||
            //       (t.vertex(1).infinite_speed == InfiniteSpeedType.WEIGHTED) ||
            //       (t.vertex(2).infinite_speed == InfiniteSpeedType.WEIGHTED));

            ///* Find vertex with the fastest incident edge,
            // * This will then gobble up one incident slower edge. */
            //assert(t.wavefront(edge_idx));
            //assert(t.wavefronts[edge_idx].vertex(0) == t.vertices[ccw(edge_idx)]);
            //assert(t.wavefronts[edge_idx].vertex(1) == t.vertices[cw (edge_idx)]);

            //assert((t.vertices[ccw(edge_idx)].infinite_speed == InfiniteSpeedType.WEIGHTED) ||
            //       (t.vertices[cw (edge_idx)].infinite_speed == InfiniteSpeedType.WEIGHTED));

            //WavefrontEdge winning_edge = t.wavefront(edge_idx);

            //WavefrontVertex v_fast;
            //KineticTriangle most_cw_triangle;
            //int idx_fast_in_most_cw_triangle;
            //int winning_edge_idx_in_v;
            ////DBG(//DBG_KT_EVENT2) << " edge_idx:  " << edge_idx;
            ////DBG(//DBG_KT_EVENT2) << " t.vertices[    edge_idx ]:  " << t.vertices[    edge_idx ];
            ////DBG(//DBG_KT_EVENT2) << " t.vertices[ccw(edge_idx)]:  " << t.vertices[ccw(edge_idx)];
            ////DBG(//DBG_KT_EVENT2) << " t.vertices[cw (edge_idx)]:  " << t.vertices[cw (edge_idx)];

            ///* If both vertices are of type InfiniteSpeedType.WEIGHTED, pick one. */
            //if (t.vertices[ccw(edge_idx)].infinite_speed == InfiniteSpeedType.WEIGHTED) {
            //  /* The left vertex of the edge is the one in question. */
            //  v_fast = t.vertices[ccw(edge_idx)];
            //  winning_edge_idx_in_v = 1;
            //  most_cw_triangle = t;
            //  idx_fast_in_most_cw_triangle = ccw(edge_idx);
            //} else {
            //  /* The right vertex of the edge is the one in question,
            //   * find the triangle that is incident to the other edge */
            //  v_fast = t.vertices[cw(edge_idx)];
            //  winning_edge_idx_in_v = 0;
            //  assert(t.wavefronts[edge_idx].vertex(1).wavefronts()[0] == t.wavefront(edge_idx));

            //  most_cw_triangle = t.wavefronts[edge_idx].vertex(1).wavefronts()[1].incident_triangle();
            //  idx_fast_in_most_cw_triangle = most_cw_triangle.index(v_fast);
            //};

            //if (t.vertices[edge_idx].is_infinite) {
            //  //DBG(//DBG_KT_EVENT2) << "Unbounded triangle";
            //  int nidx = ccw(idx_fast_in_most_cw_triangle);
            //  KineticTriangle *n = most_cw_triangle.neighbor( nidx);
            //  assert(n.unbounded());
            //  int idx_in_n = n.index(most_cw_triangle);
            //  assert(n.vertices[cw(idx_in_n)].is_infinite);

            //  //DBG(//DBG_KT_EVENT2) << "- flipping: " << most_cw_triangle;
            //  //DBG(//DBG_KT_EVENT2) << "  towards:  " << n;

            //  do_raw_flip(most_cw_triangle, nidx, time, true);
            //  modified(n);
            //} else {
            //  //DBG(//DBG_KT_EVENT2) << "Bounded triangle";

            //  /* The triangle should not have collapsed yet. */
            //  /*
            //   * Actually it may have
            //  {
            //    //DBG(//DBG_KT_EVENT2) << "most_cw_triangle is " << most_cw_triangle;
            //    const Point_2 pos_v0 = most_cw_triangle.vertex(0).p_at(time);
            //    const Point_2 pos_v1 = most_cw_triangle.vertex(1).p_at(time);
            //    const Point_2 pos_v2 = most_cw_triangle.vertex(2).p_at(time);
            //    assert(SurfNet.orientation(pos_v0, pos_v1, pos_v2) == SurfNet.LEFT_TURN);
            //  }
            //  */

            //  /* Flip away any spoke at the infinitely fast vertex. */
            //  //DBG(//DBG_KT_EVENT2) << "flipping all spokes away from " << v_fast;
            //  AroundVertexIterator flipping_triangle = incident_faces_iterator(most_cw_triangle, idx_fast_in_most_cw_triangle);
            //  int nidx_in_most_cw_triangle = ccw(idx_fast_in_most_cw_triangle);
            //  while (true) {
            //    assert(most_cw_triangle.vertices[idx_fast_in_most_cw_triangle] == v_fast);
            //    if (most_cw_triangle.is_constrained(nidx_in_most_cw_triangle)) {
            //      break;
            //    }

            //    Point_2 pos_v0 = flipping_triangle.t().vertex( ccw(flipping_triangle.v_in_t_idx()) ).p_at(time);
            //    Point_2 pos_v1 = flipping_triangle.t().vertex( cw (flipping_triangle.v_in_t_idx()) ).p_at(time);
            //    if (flipping_triangle.t() != most_cw_triangle) {
            //      /* We already delayed flipping at least once.  Let's see if we can go back */
            //      int nidx = cw(flipping_triangle.v_in_t_idx()); /* previous guy, the one cw */
            //      KineticTriangle n = flipping_triangle.t().neighbor( nidx );
            //      int idx_in_n = n.index(flipping_triangle.t());
            //      assert(n.vertex( ccw (idx_in_n) ) == flipping_triangle.t().vertex( ccw (flipping_triangle.v_in_t_idx()) ) );
            //      Point_2 pos_v2 = n.vertex(idx_in_n).p_at(time);

            //      if (SurfNet.orientation(pos_v2, pos_v0, pos_v1) != SurfNet.RIGHT_TURN) {
            //        //DBG(//DBG_KT_EVENT2) << "- We can go back, and do a flip: " << flipping_triangle.t();
            //        ++flipping_triangle;
            //        continue; /* We will flip in the next iteration. */
            //      } else {
            //        //DBG(//DBG_KT_EVENT2) << "- We cannot go back just yet";
            //      }
            //    };

            //    /* Check if we can flip to the neighbor ccw. */
            //    int nidx = ccw(flipping_triangle.v_in_t_idx()); /* next guy, the one ccw */
            //    KineticTriangle *n = flipping_triangle.t().neighbor( nidx );
            //    int idx_in_n = n.index(flipping_triangle.t());

            //    assert(n.vertex( cw (idx_in_n) ) == flipping_triangle.t().vertex( cw (flipping_triangle.v_in_t_idx()) ) );
            //    const Point_2 pos_v2 = n.vertex(idx_in_n).p_at(time);
            //    if (SurfNet.orientation(pos_v0, pos_v1, pos_v2) == SurfNet.RIGHT_TURN) {
            //      /* No, not right now.  Try in the next ccw triangle. */
            //      //DBG(//DBG_KT_EVENT2) << "- not flipping right now: " << flipping_triangle.t();
            //      --flipping_triangle;
            //    } else {
            //      //DBG(//DBG_KT_EVENT2) << "- flipping: " << flipping_triangle.t();
            //      //DBG(//DBG_KT_EVENT2) << "  towards:  " << n;
            //      do_raw_flip(flipping_triangle.t(), nidx, time, true);
            //      modified(n);
            //    }
            //  }
            //  //DBG(//DBG_KT_EVENT2) << "flipping done; " << most_cw_triangle;
            //}

            //assert(winning_edge == v_fast.wavefronts()[winning_edge_idx_in_v]);
            //const WavefrontEdge * const losing_edge = v_fast.wavefronts()[1-winning_edge_idx_in_v];
            //assert(v_fast == losing_edge.vertex(winning_edge_idx_in_v));
            //WavefrontVertex* o = losing_edge.vertex(1-winning_edge_idx_in_v);

            ////DBG(//DBG_KT_EVENT) << "v_fast " << v_fast;
            ////DBG(//DBG_KT_EVENT) << "o      " << o;
            ////DBG(//DBG_KT_EVENT) << "most_cw_triangle:  " << most_cw_triangle;
            ////DBG(//DBG_KT_EVENT) << "winning edge at v: " << winning_edge_idx_in_v;

            //if (o.infinite_speed != InfiniteSpeedType.NONE) {
            //  assert(o.infinite_speed == InfiniteSpeedType.WEIGHTED);
            //  o.stop(time, o.pos_start);
            //} else {
            //  o.stop(time);
            //}
            //v_fast.stop(time, o.pos_stop());

            //// update prev/next for the DCEL that is the wavefront vertices
            //v_fast.set_next_vertex(1-winning_edge_idx_in_v, o, false);

            //do_constraint_collapse_part2(*most_cw_triangle, most_cw_triangle.index(losing_edge), time);

            ////DBG_FUNC_END(//DBG_KT_EVENT);
        }



        private void handle_face_with_infintely_fast_opposing_vertex(CollapseEvent evnt)
        {
            throw new NotImplementedException();

            //DBG_FUNC_BEGIN(//DBG_KT_EVENT);
            //DBG(//DBG_KT_EVENT) << evnt;

            assert(evnt.type() == CollapseType.FACE_HAS_INFINITELY_FAST_VERTEX_OPPOSING);
            KineticTriangle tref = triangles[evnt.t.Id];
            KineticTriangle t = tref;
            double time = evnt.time();

            //DBG(//DBG_KT_EVENT2) << "t:  " << t;

            int num_constraints = t.is_constrained(0).ToInt() + t.is_constrained(1).ToInt() + t.is_constrained(2).ToInt();
            int num_fast = (t.vertex(0).infinite_speed != InfiniteSpeedType.NONE).ToInt() +
                                (t.vertex(1).infinite_speed != InfiniteSpeedType.NONE).ToInt() +
                                (t.vertex(2).infinite_speed != InfiniteSpeedType.NONE).ToInt();
            assert(num_fast >= 1);
            assert(num_fast < 3);
            if (num_constraints == 3)
            {
                //DBG(//DBG_KT_EVENT2) << "infinitely fast triangle with 3 constraints.";
                t.set_dying();

                Point2 p = Point2.NaN;
                bool first = true;
                for (int i = 0; i < 3; ++i)
                {
                    if (t.vertex(i).infinite_speed != InfiniteSpeedType.NONE) continue;
                    t.vertices[i].stop(time);
                    if (first)
                    {
                        p = t.vertices[i].pos_stop();
                        first = false;
                    }
                    else
                    {
                        assert(p == t.vertices[i].pos_stop());
                    }
                }
                assert(!first);
                for (int i = 0; i < 3; ++i)
                {
                    if (t.vertex(i).infinite_speed == InfiniteSpeedType.NONE) continue;
                    t.vertices[i].stop(time, p);
                }

                for (int i = 0; i < 3; ++i)
                {
                    t.wavefront(i).set_dead();
                }

                // update prev/next for the DCEL that is the wavefront vertices
                for (int i = 0; i < 3; ++i)
                {
                    t.vertices[i].set_next_vertex(0, t.vertices[cw(i)], false);
                }
                //LOG(WARNING) << __FILE__ << ":" << __LINE__ << " " << "untested code path: DECL linking.";

                queue.needs_dropping(t);
            }
            else
            {
                //DBG(//DBG_KT_EVENT2) << "infinitely fast triangle with fewer than 3 constraints.";
                assert(num_fast <= 2);
                int t_fast_idx = t.infinite_speed_opposing_vertex_idx();
                WavefrontVertex v = t.vertices[t_fast_idx];
                //DBG(//DBG_KT_EVENT2) << "infinitely fast vertex: " << t_fast_idx;

                AroundVertexIterator faces_it = incident_faces_iterator(t, t_fast_idx);
                AroundVertexIterator most_cw_triangle = faces_it.most_cw();

                //DBG(//DBG_KT_EVENT2) << "cw most triangle: " << most_cw_triangle;

                /* Flip away any spoke at the infinitely fast vertex. */
                /* ================================================== */
                //DBG(//DBG_KT_EVENT2) << "flipping all spokes away from " << most_cw_triangle.t().vertex( most_cw_triangle.v_in_t_idx() );
                List<KineticTriangle> mark_modified = new List<KineticTriangle>();
                while (true)
                {
                    int nidx = ccw(most_cw_triangle.v_in_t_idx());
                    if (most_cw_triangle.t().is_constrained(nidx))
                    {
                        break;
                    }
                    KineticTriangle n = most_cw_triangle.t().neighbor(nidx);
                    //DBG(//DBG_KT_EVENT2) << "- flipping: " << most_cw_triangle.t();
                    //DBG(//DBG_KT_EVENT2) << "  towards:  " << n;
                    do_raw_flip(most_cw_triangle.t(), nidx, time, true);
                    mark_modified.Add(n);
                }
                //DBG(//DBG_KT_EVENT2) << "flipping done; " << most_cw_triangle;

                t = most_cw_triangle.t();
                //DBG(//DBG_KT_EVENT2) << "most cw triangle is: " << t;
                int vidx_in_tc = most_cw_triangle.v_in_t_idx();

                assert(t.vertex(vidx_in_tc) == v);
                assert(t.is_constrained(cw(vidx_in_tc)));
                assert(t.is_constrained(ccw(vidx_in_tc)));

                /* Figure out which edge to retire */
                /* =============================== */
                WavefrontVertex v_cw = t.vertices[cw(vidx_in_tc)];
                WavefrontVertex v_ccw = t.vertices[ccw(vidx_in_tc)];

                WavefrontEdge l = t.wavefront(ccw(vidx_in_tc));
                WavefrontEdge r = t.wavefront(cw(vidx_in_tc));

                //DBG(//DBG_KT_EVENT) << "l is " << *l;
                //DBG(//DBG_KT_EVENT) << "r is " << *r;

                //DBG(//DBG_KT_EVENT) << "v     is " << v;
                //DBG(//DBG_KT_EVENT) << "v_cw  is " << v_cw;
                //DBG(//DBG_KT_EVENT) << "v_ccw is " << v_ccw;

                assert((v_cw.infinite_speed == InfiniteSpeedType.NONE) ||
                       (v_ccw.infinite_speed == InfiniteSpeedType.NONE));

                int collapse = 0;
                WavefrontVertex o = null;
                bool spoke_collapse = false;
                /* collapse the shorter edge, or the edge to the non-fast vertex (i.e, the one opposite of the fast). */
                if (v_cw.infinite_speed != InfiniteSpeedType.NONE)
                {
                    collapse = cw(vidx_in_tc);
                    o = v_ccw;
                    //DBG(//DBG_KT_EVENT) << "v_cw has infinite speed";
                }
                else if (v_ccw.infinite_speed != InfiniteSpeedType.NONE)
                {
                    collapse = ccw(vidx_in_tc);
                    o = v_cw;
                    //DBG(//DBG_KT_EVENT) << "v_ccw has infinite speed";
                }
                else
                {
                    var pos = v.p_at(time);
                    var poscw = v_cw.p_at(time);
                    var posccw = v_ccw.p_at(time);
                    double sq_length_v_to_vcw = squared_distance(pos, poscw);
                    double sq_length_v_to_vccw = squared_distance(pos, posccw);
                    if (sq_length_v_to_vcw < sq_length_v_to_vccw)
                    {
                        collapse = ccw(vidx_in_tc);
                        o = v_cw;
                        //DBG(//DBG_KT_EVENT) << "sq_length_v_to_vcw < sq_length_v_to_vccw";
                    }
                    else if (sq_length_v_to_vcw > sq_length_v_to_vccw)
                    {
                        collapse = cw(vidx_in_tc);
                        o = v_ccw;
                        //DBG(//DBG_KT_EVENT) << "sq_length_v_to_vcw > sq_length_v_to_vccw";
                    }
                    else
                    {
                        //DBG(//DBG_KT_EVENT) << "sq_length_v_to_vcw == sq_length_v_to_vccw";
                        spoke_collapse = true;
                    }
                }
                if (spoke_collapse)
                {
                    //DBG(//DBG_KT_EVENT) << "both edges incident to the infinitely fast vertex have the same length";
                    t.set_dying();

                    v_cw.stop(time);
                    v_ccw.stop(time);
                    assert(v_cw.pos_stop() == v_ccw.pos_stop());

                    v.stop(time, v_cw.pos_stop());
                    assert(t.wavefront(cw(vidx_in_tc)) != null);
                    assert(t.wavefront(ccw(vidx_in_tc)) != null);
                    t.wavefront(cw(vidx_in_tc)).set_dead();
                    t.wavefront(ccw(vidx_in_tc)).set_dead();

                    KineticTriangle n = t.neighbors[vidx_in_tc];
                    assert(n != null);
                    int idx_in_n = n.index(t);

                    t.neighbors[vidx_in_tc] = null;
                    n.neighbors[idx_in_n] = null;

                    // update prev/next for the DCEL that is the wavefront vertices
                    v.set_next_vertex(0, v_cw, false);
                    v.set_next_vertex(1, v_ccw, false);

                    assert(!t.is_dead());
                    do_spoke_collapse_part2(n, idx_in_n, time);

                    assert(!t.is_dead());
                    queue.needs_dropping(t);
                }
                else
                {
                    //DBG(//DBG_KT_EVENT) << "collapse " << collapse;
                    //DBG(//DBG_KT_EVENT) << "v " << v;
                    //DBG(//DBG_KT_EVENT) << "o " << o;
                    assert((t.index(v) + t.index(o) + collapse) == 3);

                    o.stop(time);
                    v.stop(time, o.pos_stop());

                    // update prev/next for the DCEL that is the wavefront vertices
                    assert(o == v_ccw || o == v_cw);
                    v.set_next_vertex(o == v_ccw ? 1 : 0, o, false);
                    //LOG(WARNING) << __FILE__ << ":" << __LINE__ << " " << "untested code path: DECL linking.";

                    do_constraint_collapse_part2(t, collapse, time);
                }
                foreach (KineticTriangle tm in mark_modified)
                {
                    if (!tm.is_dead())
                    {
                        modified(tm);
                    }
                    else
                    {
                        //DBG(//DBG_KT_EVENT) << "Not marking " << tm << " as modified because it is dead already.";
                        assert(spoke_collapse);
                    }
                }
            }
            assert_valid(t.component, time);

            //DBG_FUNC_END(//DBG_KT_EVENT);
        }


        /// <summary>
        /// Transfer constraints to neighbors
        ///
        /// Face f is going away because it collapsed.  Transfer the constraints
        /// where appropriate.
        ///
        /// f is a triangle with edges a, b, c.  edge a has collapsed.  Not both
        /// of b and c may be constrained.
        ///
        /// If either b or c is constrained, move the constraint to the neighbor
        /// on the other side.
        /// </summary>
        private void move_constraints_to_neighbor(KineticTriangle t, int idx)
        {
            assert(! (t.is_constrained(cw(idx)) && t.is_constrained(ccw(idx))) );

            for (int i = 1; i <= 2; ++i)
            {
                if (t.is_constrained((idx + i).mod3()))
                {
                    var n = t.neighbor((idx + 3 - i).mod3());
                    int nidx = n.index(t);
                    assert(! n.is_constrained(nidx) );

                    n.move_constraint_from(nidx, t, (idx + i).mod3());
                }
            }
        }



        ///<summary>
        /// Note that triangle t has been modified.
        /// Put it in out check_refinement queue (at front or back), and
        /// put it in the needs_update queue of the evnt queue.
        ///</summary>
        private void modified(KineticTriangle t, bool front=false)
        {
            put_on_check_refinement(t, front);
            // In the initial refinement calls, we will
            // not have a queue yet.
            if (queue != null)
            {
                queue.needs_update(t);
            }

            if (t.unbounded())
            {
                /* t's neighbor also witnesses the common vertex remaining in the wavefront.
                 * let them know things might have changed.
                 */
                int idx = t.infinite_vertex_idx();
                KineticTriangle n = t.neighbor(ccw(idx));
                assert(n != null);
                if (!n.is_dying())
                {
                    n.invalidate_collapse_spec();
                    if (queue != null)
                    {
                        queue.needs_update(n);
                    }
                }
            }
        }

        ///////////////
        // Refinement things
        private List<KineticTriangle> check_refinement = new TriangleList();
        private bool[] tidx_in_check_refinement;

        /** refines the triangulation locally around t
         *
         * refines the triangulation locally around t with at most one flip.
         */
#if !REFINE_TRIANGULATION

        private void refine_triangulation(KineticTriangle t, double time)
        {
            //DBG_FUNC_BEGIN(//DBG_KT_EVENT);
            //DBG(//DBG_KT_EVENT) << t;

            assert(t != null);
            t.assert_valid();

            if (!t.unbounded())
            {
                WavefrontVertex v0 = t.vertex(0);
                WavefrontVertex v1 = t.vertex(1);
                WavefrontVertex v2 = t.vertex(2);

                int num_reflex = (!v0.is_convex_or_straight()).ToInt() +
                                 (!v1.is_convex_or_straight()).ToInt() +
                                 (!v2.is_convex_or_straight()).ToInt();
                switch (num_reflex)
                {
                    case 0:
                        break;

                    case 1:
                        {
                            int reflex_idx = (!v0.is_convex_or_straight()) ? 0 :
                                             (!v1.is_convex_or_straight()) ? 1 :
                                                                              2;
                            assert(t.vertex(cw(reflex_idx)).is_convex_or_straight());
                            assert(t.vertex(ccw(reflex_idx)).is_convex_or_straight());

                            if (t.is_constrained(reflex_idx)) break;
                            KineticTriangle n = t.neighbor(reflex_idx);
                            if (n.unbounded()) break;

                            /* If on of the corners of the quadrilateral is actually straight (it won't be reflex), do not flip */
                            int idx_in_n = n.index(t);
                            //if (t.is_constrained(ccw(reflex_idx)) && n.is_constrained(cw (idx_in_n)) && t.vertex(cw (reflex_idx)).is_reflex_or_straight()) break;
                            //if (t.is_constrained(cw (reflex_idx)) && n.is_constrained(ccw(idx_in_n)) && t.vertex(ccw(reflex_idx)).is_reflex_or_straight()) break;

                            WavefrontVertex v = t.vertices[reflex_idx];
                            WavefrontVertex va = t.vertices[ccw(reflex_idx)];
                            WavefrontVertex vb = t.vertices[cw(reflex_idx)];
                            WavefrontVertex o = n.vertices[idx_in_n];
                            for (int i = 0; i <= 1; ++i)
                            {
                                assert(v.wavefronts()[i] != null);
                                assert(va.wavefronts()[i] != null);
                                assert(vb.wavefronts()[i] != null);
                                assert(o.wavefronts()[i] != null);

                                assert(v.wavefronts()[i].vertex(1 - i) == v);
                                assert(va.wavefronts()[i].vertex(1 - i) == va);
                                assert(vb.wavefronts()[i].vertex(1 - i) == vb);
                                assert(o.wavefronts()[i].vertex(1 - i) == o);
                            }
                            /* If on of the corners of the quadrilateral is actually straight (it won't be reflex), do not flip */
                            if (v.wavefronts()[1].vertex(1) == va && va.is_reflex_or_straight()) break;
                            if (v.wavefronts()[0].vertex(0) == vb && vb.is_reflex_or_straight()) break;
                            /* Either at v, or at the opposite vertex */
                            // //DBG(//DBG_KT_REFINE) << "   o.wavefronts()[0].vertex(0) " << o.wavefronts()[0].vertex(0);
                            // //DBG(//DBG_KT_REFINE) << "   o.wavefronts()[0].vertex(1) " << o.wavefronts()[0].vertex(1);
                            // //DBG(//DBG_KT_REFINE) << "   o.wavefronts()[1].vertex(0) " << o.wavefronts()[1].vertex(0);
                            // //DBG(//DBG_KT_REFINE) << "   o.wavefronts()[1].vertex(1) " << o.wavefronts()[1].vertex(1);
                            if (o.wavefronts()[0].vertex(0) == va && va.is_reflex_or_straight()) break;
                            if (o.wavefronts()[1].vertex(1) == vb && vb.is_reflex_or_straight()) break;

                            //DBG(//DBG_KT_REFINE) << "  Flipping " << t.get_name() << " along " << reflex_idx;
                            //DBG(//DBG_KT_REFINE) << "   t: " << t;
                            //DBG(//DBG_KT_REFINE) << "   n: " << n;
                            //DBG(//DBG_KT_REFINE) << "   v[" <<     reflex_idx  << "]: " << v  << " - " << CGAL_point(v .p_at(time));
                            //DBG(//DBG_KT_REFINE) << "   v[" << ccw(reflex_idx) << "]: " << va << " - " << CGAL_point(va.p_at(time));
                            //DBG(//DBG_KT_REFINE) << "   v[" << cw (reflex_idx) << "]: " << vb << " - " << CGAL_point(vb.p_at(time));
                            //DBG(//DBG_KT_REFINE) << "   n[o"                      "]: " << o  << " - " << CGAL_point(o .p_at(time));

                           
                            do_flip(t, reflex_idx, time);
                           
                            break;
                        };
                    case 2:
                        break;

                    case 3:
                        break;

                        //default:
                        //  assert(false);
                }
            }

            //DBG_FUNC_END(//DBG_KT_EVENT);
        }



#else
        private partial void refine_triangulation(KineticTriangle t, double time)
        {
        }
#endif
       

       
        

        public  void refine_triangulation_initial()
        {
            foreach (var t in triangles)
            {
                // Do _one_ refinement of each triangle.  If it actually does change
                // anything, then refine_triangulation() will have put the triangle
                // onto the check_refinement queue.  If not, no harm done.
                refine_triangulation(t, CORE_ZERO);
            }
            process_check_refinement_queue(CORE_ZERO);
            assert_valid(-1, CORE_ZERO);
            IsIntialRefineTriangleDoit = true;  
        }

        private partial void process_check_refinement_queue(double time);
        private partial void put_on_check_refinement(KineticTriangle t, bool front = false);
        private partial KineticTriangle check_refinement_pop();
        ///////////////

        private static double handle_event_last_time = 0;
        private static int handle_event_count = 0;

        public  void handle_event(CollapseEvent evnt)
        {
            //DBG_FUNC_BEGIN(//DBG_KT | //DBG_KT_EVENT);
            //DBG(//DBG_KT_EVENT);
            //DBG(//DBG_KT_EVENT);
            //DBG(//DBG_KT_EVENT);
            //DBG(//DBG_KT | //DBG_KT_EVENT) << evnt;
            assert(!finalized);

            Log($"Handle:{evnt}");

            double time = evnt.time();

            if (time <= handle_event_last_time)
            {
                ++handle_event_count;
                if (handle_event_count > 10000)
                {
                    throw new Exception("In double loop at line ");
                }
            }
            else
            {
                handle_event_count = 0;
                handle_event_last_time = time;
            };

            assert(triangles[evnt.t.Id] == evnt.t);

            assert(check_refinement.Count == 0);

            assert((CollapseSpec)(evnt) == evnt.t.get_collapse(time));

            ++event_type_counter[(int)(CollapseType.UNDEFINED)];
            ++event_type_counter[(int)(evnt.type())];
            update_event_timing_stats(time);

            switch (evnt.type())
            {
                case CollapseType.TRIANGLE_COLLAPSE:
                    handle_triangle_collapse_event(evnt);

                    break;

                case CollapseType.CONSTRAINT_COLLAPSE:
                    handle_constraint_event(evnt);
                    //throw new NotImplementedException();
                    break;

                case CollapseType.FACE_HAS_INFINITELY_FAST_VERTEX_OPPOSING:

                    handle_face_with_infintely_fast_opposing_vertex(evnt);
                    break;

                case CollapseType.FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED:
                    handle_face_with_infintely_fast_weighted_vertex(evnt);
                    break;

                case CollapseType.SPOKE_COLLAPSE:
                    handle_spoke_collapse_event(evnt);
                    break;

                case CollapseType.SPLIT_OR_FLIP_REFINE:
                    handle_split_or_flip_refine_event(evnt);
                    break;

                case CollapseType.VERTEX_MOVES_OVER_SPOKE:
                    handle_vertex_moves_over_spoke_event(evnt);
                    break;

                case CollapseType.CCW_VERTEX_LEAVES_CH:
                    handle_ccw_vertex_leaves_ch_event(evnt);
                    break;
                /*
                case CollapseType.GENERIC_FLIP_EVENT:
                  assert(false);
                  exit(1);
                  break;
                */
                case CollapseType.INVALID_EVENT:
                case CollapseType.UNDEFINED:
                case CollapseType.NEVER:
                    throw new Exception($"Should not get evnt {evnt} to handle.");
                    break;

                default:
                    throw new Exception($"Unexpected evnt {evnt}");
                    break;
            }
            //DBG(//DBG_KT_REFINE) << evnt << " - done.  doing refinement";
            process_check_refinement_queue(time);

            //DBG_FUNC_END(//DBG_KT | //DBG_KT_EVENT);
        }
        public SkeletonDCEL get_skeleton() { return skeleton; }


        public static List<List<Point2>> CargarListaVector2(string rutaArchivo)
        {
            if (!File.Exists(rutaArchivo))
                throw new FileNotFoundException("El archivo no existe.", rutaArchivo);

            string json = File.ReadAllText(rutaArchivo);
            return JsonSerializer.Deserialize<List<List<Point2>>>(json);
        }
    }


}