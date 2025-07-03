namespace SurfNet
{
    using static DebugLog;
    using static Mathex;

    /** Set to -1 for slower edge wins.  Changing this during the propagation will
 * yield undefined behaviour as the event queue is not correct then.  */

    /* A triangle, part of our kinetic triangulation.
     *
     * We maintain a kinetic triangulation of the area not yet swept by the
     * wavefront.  When sweeping to the outside, this area can be unbounded.
     *
     * - A triangle of this triangulation always has three vertices (one of which
     *   may be the "infinite" vertex.
     * - Triangle edges may be part of the wavefront, in which case we point to the
     *   input edge from which the corresponding wavefront edge emanated.
     * - Also, each triangle has up to three neighbors, exactly one for each
     *   non-constrained edge.  Thus, each side either has a neighbor or a wavefront set.
     *
     * As such, the kinetic triangulation is a triangulation of a set of polygons,
     * one of them unbounded.
     */

    public partial class KineticTriangle
    {
        //  friend class KineticTriangulation;


        public readonly int component;

        private bool is_dead_ = false;
        private bool is_dying_ = false;
        public HalfEdge3 Edges = null;
        public WavefrontVertex3 vertices = new WavefrontVertex3 { };
        public WavefrontEdge3 wavefronts = new WavefrontEdge3 { };
        public KineticTriangle3 neighbors = new KineticTriangle3 { };




        //        public IEnumerable<KineticTriangulation.HalfEdge> Edges { 
        //            get {
        //               if (edge == null) yield break;
        //                var e = edge;

        //#if !NO_PRECONDITION
        //                int i= 0;
        //#endif
        //                do
        //                {
        //#if !NO_PRECONDITION
        //                    if (i > 2) throw new Exception("contiene mas de 3 lados");
        //                    i++;
        //#endif

        //                    yield return e;
        //                    e = e.Next;

        //                } while (e != null && edge != e);
        //                yield break;

        //            } 

        //        }




        private CollapseSpec collapse_spec;
        private bool collapse_spec_valid = false;
#if !SURF_NDEBUG
        private WavefrontVertex3 collapse_spec_computed_with_vertices = new WavefrontVertex3();

#endif





        public static int cw(int i)
        { return TriangulationUtils.cw(i); }

        public static int ccw(int i)
        { return TriangulationUtils.ccw(i); }

        public readonly TriangleNet.Geometry.ITriangle OriginalTriangle;

        public int Id { get; private set; }

        public KineticTriangle(int pd_id, int p_component, TriangleNet.Geometry.ITriangle originalTriangle = null)

        {
            OriginalTriangle = originalTriangle;
            Id = pd_id;
            component = p_component;
            vertices = new WavefrontVertex3();
            //  edge = null;
            wavefronts = new WavefrontEdge3();
            neighbors = new KineticTriangle3();
            collapse_spec = new CollapseSpec(p_component, CollapseType.NEVER);
            collapse_spec_valid = false;

#if !SURF_NDEBUG
            collapse_spec_computed_with_vertices = new WavefrontVertex3();
#endif
        }

        //    KineticTriangle(const KineticTriangle&) = delete;
        //KineticTriangle& operator =(const KineticTriangle&) = delete;
        //KineticTriangle(KineticTriangle&&) = default;
        //KineticTriangle& operator = (KineticTriangle&&) = default;

        internal partial void set_neighbors(params KineticTriangle?[] n);

        internal void set_wavefronts(WavefrontEdge? w0, WavefrontEdge? w1, WavefrontEdge? w2)
        {
            wavefronts[0] = w0;
            wavefronts[1] = w1;
            wavefronts[2] = w2;
            invalidate_collapse_spec();
        }

        /** Set a vertex of the triangle, updating the wavefront edge if this side is constrained. */

        internal partial void set_vertex(int i, WavefrontVertex v);

        internal partial bool has_neighbor(KineticTriangle needle);

        internal partial bool has_vertex(WavefrontVertex needle);

        internal partial bool has_wavefront(WavefrontEdge needle);

        internal partial int index(KineticTriangle needle);

        internal partial int index(WavefrontVertex needle);

        internal partial int index(WavefrontEdge needle);

        internal partial KineticTriangle neighbor(int i);

        /** Check if this side has a constraint.
         *
         * This usually implies it has no neighbor, except during construction and while
         * manipulating the triangulation.  Outside of methods, it should always hold
         * that we have exactly either neighbor[i] or constraint[i].
         */

        internal partial bool is_constrained(int i);

        internal partial WavefrontEdge wavefront(int i);

        internal WavefrontVertex vertex(int i) { assert(i < 3); return vertices[i]; }

#if !SURF_NDEBUG

        internal partial void assert_valid();

#else
    void assert_valid() {};
#endif

        //friend std::ostream& operator<<(std::ostream& os, const KineticTriangle  const kt);
        //std::string get_name() const {
        //  return "kt" + std::to_string(id);
        //}

        private void assert_is_id(int q)
        { assert(Id == q); }

        public CollapseSpec get_collapse(double time_now)
        {
            assert(!is_dead_);
            assert(!is_dying_);
            if (!collapse_spec_valid)
            {
                collapse_spec = compute_collapse(time_now);
                collapse_spec_valid = true;
#if !SURF_NDEBUG
                set_to_cur_wf_vertices(collapse_spec_computed_with_vertices);
#endif
            };
            return get_cached_collapse();
        }

        internal CollapseSpec refine_collapse_spec(CollapseSpec c)
        {

            assert(this.collapse_spec_valid);
#if !SURF_NDEBUG
            assert_cur_wf_vertices(collapse_spec_computed_with_vertices);
#endif
            assert(this.collapse_spec.allows_refinement_to(c));

            this.collapse_spec = c;
            return collapse_spec;
        }

        private CollapseSpec get_cached_collapse()
        {
            assert(!is_dying_);
            assert(collapse_spec_valid);
#if !SURF_NDEBUG
            assert_cur_wf_vertices(collapse_spec_computed_with_vertices);
#endif
            return collapse_spec;
        }

        internal void invalidate_collapse_spec()
        {
            assert(!is_dead_);
            collapse_spec_valid = false;
#if !SURF_NDEBUG
            invalidate_cur_wf_vertices(collapse_spec_computed_with_vertices);
#endif
            /*
            _is_squeezed = false;
            determinant_valid = false;
            invalidate_cur_wf_vertices(determinantComputedWithWavefrontVertices);
            */
        }

#if !SURF_NDEBUG

        private partial void set_to_cur_wf_vertices(WavefrontVertex3 arr);

        private partial void invalidate_cur_wf_vertices(WavefrontVertex3 arr);

        private partial void assert_cur_wf_vertices(WavefrontVertex3 arr);

#endif

        private partial CollapseSpec compute_collapse(double time_now);

        private partial CollapseSpec compute_collapse_bounded(double time_now);

        private partial CollapseSpec compute_collapse_unbounded(double time_now);

        private partial CollapseSpec compute_collapse_bounded_constrained_0(double time_now);

        private partial CollapseSpec compute_collapse_bounded_constrained_1(double time_now);

        private partial CollapseSpec compute_collapse_bounded_constrained_2(double time_now);

        private partial CollapseSpec compute_collapse_bounded_constrained_3(double time_now);

        private static partial Polynomial1D compute_determinant_from_vertices(WavefrontVertex v0, WavefrontVertex v1, WavefrontVertex v2);

        private static partial CollapseSpec event_that_will_not_happen(int component, double time_now, Polynomial1D determinant);

        private static partial bool get_generic_collapse_time(double time_now, Polynomial1D det, out double collapse_time);

        private static partial bool accept_collapse_bounded_constrained_1(double collapse_time, Polynomial1D determinant, bool collapse_time_is_edge_collapse);

        private partial CollapseSpec get_generic_collapse(double time_now, Polynomial1D determinant);

        // CollapseSpec compute_constraint_collapse(const NT& time_now, const Polynomial_1& determinant) const;
        // CollapseSpec compute_split_event(const NT& time_now) const;
        // CollapseSpec determine_split_or_flip_bounded_constrained_1(const NT& collapse_time, int c_idx) const;
        private partial CollapseSpec compute_split_or_flip_event_bounded_constrained_1(double time_now, int c_idx, Polynomial1D det);

        private partial CollapseSpec compute_flip_event(double time_now, Polynomial1D determinant);

        /* called from KineticTriangulation only */

        internal partial void move_constraint_from(int idx, KineticTriangle src, int src_idx);

        internal void set_dying(){ 
            is_dying_ = true;
            
            LogRich($"kt{Id} is dying v:[{ string.Join(", ",vertices.Select(v=>v.has_stopped()?$"[s]{v.Id}[/s]":$"[b]{v.Id}[/b]"))}]");
        }

        private partial void set_neighbor(int idx, KineticTriangle n);

        /** set wavefront and update wavefront's incident triangle.
         *
         * - wavefront needs to be not null,
         * - wavefront's vertices already need to be right,
         * - there needs to be a neighbor on the side that is about to become a
         *   wavefront, and
         * - that neighbor needs to be dying.
         * - wavefront's (old) incident_triangle needs to be neighbor
         *
         * - sets neighbor to null
         * - sets wavefront
         *   sets wavefront's incident_triangle
         */

        internal void set_wavefront(int idx, WavefrontEdge e)
        { // {{{
            //CGAL_precondition(idx < 3);
            assert(e != null);

            assert(wavefronts[idx] == null);
            assert(neighbors[idx] != null);
            assert(neighbors[idx].is_dying());
            assert(e.incident_triangle() == neighbors[idx]);

            wavefronts[idx] = e;
            neighbors[idx] = null;
            e.set_incident_triangle(this);
        } // }}}

        



        public bool is_dead() => is_dead_;

        /* called by EventQ */

        public bool is_dying() => is_dying_;

        /** Mark this triangle as dead.  May only be called once. */

        public partial void set_dead();

        public bool is_collapse_spec_valid() => collapse_spec_valid;

        public partial bool unbounded();

        public partial int infinite_vertex_idx();

        public partial InfiniteSpeedType has_vertex_infinite_speed();

        public partial int infinite_speed_opposing_vertex_idx();

        public enum VertexOnSupportingLineType
        { ONCE, NEVER, ALWAYS };

        //    friend std::ostream& operator<<(std::ostream& os, const KineticTriangle::VertexOnSupportingLineType a);

        public (double, VertexOnSupportingLineType) get_time_vertex_on_supporting_line(WavefrontVertex v, WavefrontSupportingLine e)
        { 
            double collapse_time;
            VertexOnSupportingLineType vertex_on_line_type;

            Log($" ");

            /* Let n be some normal to e,
             * let P be some point on e (at time zero),
             * let s be the speed (well, velocity) vector of the vertex v.
             * let Q be the location of v at time zero, and
             * let w be the weight (speed) of e.
             */
            Vector2 n = (e.normal_direction);
            Point2 P = (e.l.point());
            Vector2 s = (v.velocity);
            Point2 Q = (v.pos_zero);
            double w = (e.weight);
            /* Then PQ.n is the length of the projection of PQ onto n, times the length of n.
             * Likewise, s.n is the length of the projection of s onto n, times the length of n.
             * Per time unit, v and e will approach each other by (w - s.n/|n|).
             *
             * So v will hit (the supporting line of) e at time t := PQ.n/|n| / (w - s.n/|n|) ==
             * == PQ.n / (w |n| - s.n)
             */
            Vector2 PQ = new Vector2(P, Q);
            double scaled_distance = PQ.Dot(n);
            double scaled_edge_speed = w * sqrt(n.squared_length());
            double scaled_vertex_speed = s.Dot(n);
            double scaled_speed_approach = scaled_edge_speed - scaled_vertex_speed;

            //DBG(//DBG_TRIANGLE_TIMING2) << " -- n: " << CGAL_vector(n);
            //DBG(//DBG_TRIANGLE_TIMING2) << " -- P: " << CGAL_point(P);
            //DBG(//DBG_TRIANGLE_TIMING2) << " -- s: " << CGAL_vector(s);
            //DBG(//DBG_TRIANGLE_TIMING2) << " -- Q: " << CGAL_point(Q);
            //DBG(//DBG_TRIANGLE_TIMING2) << " -- w: " << to_double(w);
            //DBG(//DBG_TRIANGLE_TIMING2) << " -- PQ: " << CGAL_vector(PQ);
            //DBG(//DBG_TRIANGLE_TIMING2) << " -- num (∝ distance      ): " << to_double(scaled_distance);
            //DBG(//DBG_TRIANGLE_TIMING2) << " -- den (∝ approach speed): " << to_double(scaled_speed_approach);

            if (scaled_speed_approach == 0)
            {
                collapse_time = CORE_ZERO;
                if (scaled_distance == 0)
                {
                    vertex_on_line_type = VertexOnSupportingLineType.ALWAYS;
                }
                else
                {
                    vertex_on_line_type = VertexOnSupportingLineType.NEVER;
                }
            }
            else
            {
                collapse_time = scaled_distance / scaled_speed_approach;
                vertex_on_line_type = VertexOnSupportingLineType.ONCE;
            }

            var areaT0 = area(0.0);
            var areaTC = area(collapse_time);
            areaT0 = !areaT0.IsZero() ? areaT0 : 0.0;
            areaTC = !areaTC.IsZero() ? areaTC : 0.0;
            assert(sign(areaT0) == sign(areaTC) || areaTC.IsZero(), $" kt:{Id} a0:{areaT0} at: {areaTC} time:{collapse_time}");
            //DBG(//DBG_TRIANGLE_TIMING2) << "returning " << to_double(collapse_time) << " with VertexOnSupportingLineType " << vertex_on_line_type;
            //DBG_FUNC_END(//DBG_TRIANGLE_TIMING2);
            return (collapse_time, vertex_on_line_type);
        } 
        public double area(double time = 0.0)
        {
            var p0 = vertices[0].p_at(time);
            var p1 = vertices[1].p_at(time);
            var p2 = vertices[2].p_at(time);

            return (p1.X - p0.X) * (p2.Y - p0.Y) - (p1.Y - p0.Y) * (p2.X - p0.X);
        }
        public static partial int edge_is_faster_than_vertex(WavefrontVertex v, WavefrontSupportingLine e);
    }

    //std::ostream& operator<<(std::ostream& os, const KineticTriangle::VertexOnSupportingLineType a);
}