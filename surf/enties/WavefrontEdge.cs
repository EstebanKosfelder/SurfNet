namespace SurfNet
{
    using static DebugLog;

    public class WavefrontEdgeList : List<WavefrontEdge> { };

    public partial class WavefrontEdge
    {

        public override string ToString()
        {
            return $" we:{Id:##0} {vertices[0]?.Id:##0}->{vertices[1]?.Id:##0} {supporting_line.l.ToVector().Normal().Debug()} ";
        }

        private static int wavefront_edge_ctr;
        public readonly int Id;


        private bool is_dead_ = false;        /** stopped propagating */

        /** The left and right wavefront vertex right
                                   *  now.  This changes over the course of the
                                   *  propagation period.
                                   */

        ///<summary>
        /// The left and right wavefront vertex right now.
        /// This changes over the course of the propagation period.
        ///</summary>
        private WavefrontVertex[] vertices = new WavefrontVertex[2];

        ///<summary>
        ///  The supporting line backing this wavefront vertex.
        ///</summary>
        private WavefrontSupportingLine supporting_line;

        ///<summary>
        /// The triangle incident right now at this wavefront edge.
        ///</summary>
        private KineticTriangle? incident_triangle_ = null;


        /// <summary>
        /// Is this wavefront edge one that was
        /// created initially as part of the input or
        /// from beveling (true), or is it the result
        /// of a wavefront edge having been split during the
        /// propagation period (false).
        /// </summary>
        public readonly bool is_initial;

        /// <summary>
        /// Is this wavefront edge the result of beveling,
        /// and thus degenerate at time zero?
        /// </summary>
        public readonly bool is_beveling;


        /// <summary>
        /// The pointers to the left and right straight
        /// skeleton arcs (==kinetic wavefront vertex).
        /// Only for is_initial wavefront edges, so
        /// we can then find the faces nicely. 
        /// </summary>
        private WavefrontVertex[] initial_vertices = new WavefrontVertex[2];


        ///<summary>
        /// The straight skeleton face that this edge traces out (at least partially.
        /// With splitting, multiple edges are needed to trace out a single face.
        ///</summary>
        public readonly SkeletonDCELFace skeleton_face;

        private EdgeCollapseSpec collapse_spec;
        private bool collapse_spec_valid = false;
#if  !SURF_NDEBUG
        private WavefrontVertex[] collapse_spec_computed_with_vertices = new WavefrontVertex[2];
#endif

        /// <summary>
        /// Used when setting up initial wavefront edges for all constraints 
        /// </summary>
        public WavefrontEdge(Point2 u, Point2 v, double weight, KineticTriangle incident_triangle, SkeletonDCELFace p_skeleton_face)

        {
#if !SURF_NDEBUG
            Id = (wavefront_edge_ctr++);
#endif
            vertices = new WavefrontVertex[] { null, null };
            supporting_line = new WavefrontSupportingLine(u, v, weight);
            incident_triangle_ = (incident_triangle);
            is_initial = (true);
            is_beveling = (false);
            initial_vertices = new WavefrontVertex[] { null, null };
            skeleton_face = p_skeleton_face;
#if !SURF_NDEBUG
            collapse_spec_computed_with_vertices = new WavefrontVertex[] { null, null };
#endif
            //     assert((skeleton_face!=null) ^ is_beveling);
        }

        /// <summary>
        /// Used when setting up bevels 
        /// </summary>
        public WavefrontEdge(WavefrontSupportingLine p_supporting_line)
        {
#if !SURF_NDEBUG
            Id = (wavefront_edge_ctr++);
#endif
            vertices = new WavefrontVertex[] { null, null };
            supporting_line = p_supporting_line;
            incident_triangle_ = null;
            is_initial = true;
            is_beveling = true;
            initial_vertices = new WavefrontVertex[] { null, null };
            skeleton_face = null;
# if !SURF_NDEBUG
            collapse_spec_computed_with_vertices = new WavefrontVertex[] { null, null };
#endif

          //  assert(skeleton_face != null ^ is_beveling);
        }

        private WavefrontEdge(WavefrontVertex va,
                          WavefrontVertex vb,
                          WavefrontSupportingLine p_supporting_line,
                          KineticTriangle incident_triangle,
                          SkeletonDCELFace p_skeleton_face,
                          bool p_is_beveling)
        {
# if !SURF_NDEBUG
            Id = (wavefront_edge_ctr++);
#endif
            vertices = new WavefrontVertex[] { va, vb };
            supporting_line = (p_supporting_line);
            incident_triangle_ = (incident_triangle);
            is_initial = (false);
            is_beveling = (p_is_beveling);
            initial_vertices = new WavefrontVertex[] { null, null };
            skeleton_face = (p_skeleton_face);
#if !SURF_NDEBUG
            collapse_spec_computed_with_vertices = new WavefrontVertex[] { null, null };
#endif

          //  assert((skeleton_face != null) ^ is_beveling);
        }

        public void set_dead()
        {
            assert(!is_dead_);
            is_dead_ = true;
        }

        public WavefrontSupportingLine l() => supporting_line;

        public KineticTriangle incident_triangle() => incident_triangle_;

        public bool is_dead() => is_dead_;

        public WavefrontVertex vertex(int i)
        {
            assert(!is_dead_);
            assert(i <= 1);
            return vertices[i];
        }

        public WavefrontVertex initial_vertex(int i)
        {
            assert(is_initial);
            assert(i <= 1);
            return initial_vertices[i];
        }

        public void set_wavefrontedge_vertex(int i, WavefrontVertex v)
        {
            assert(!is_dead_);
            assert(i <= 1);
            assert(v != null);
            vertices[i] = v;
            invalidate_collapse_spec();
        }

        internal void set_initial_vertices()
        {
            assert(is_initial);
            for (int i = 0; i <= 1; ++i)
            {
                assert(vertices[i] != null);
                assert(initial_vertices[i] == null);
                initial_vertices[i] = vertices[i];
            };
        }

        internal void set_initial_incident_triangle(KineticTriangle incident_triangle)
        {
            incident_triangle_ = incident_triangle;
            int idx = incident_triangle.index(this);
            assert(incident_triangle.vertex(KineticTriangle.ccw(idx)) == vertices[0]);
            assert(incident_triangle.vertex(KineticTriangle.cw(idx)) == vertices[1]);
        }
        public void set_incident_triangle(KineticTriangle incident_triangle)
        {
            assert(!is_dead_);
            assert(incident_triangle.has_wavefront(this));
            int idx = incident_triangle.index(this);

            incident_triangle_ = incident_triangle;

            assert(incident_triangle.vertex(KineticTriangle.ccw(idx)) == vertices[0]);
            assert(incident_triangle.vertex(KineticTriangle.cw(idx)) == vertices[1]);
            invalidate_collapse_spec();
        }

        public CollapseSpec get_collapse(int component, double time_now, int collapsing_edge)
        {
            assert_edge_sane(collapsing_edge);
            return new CollapseSpec(component, get_edge_collapse(time_now), collapsing_edge);
        }

        public EdgeCollapseSpec get_edge_collapse(double time_now)
        {
            assert(!is_dead_);
            assert(vertices[0] != null);
            assert(vertices[1] != null);
            if (!collapse_spec_valid)
            {
                collapse_spec = compute_collapse(time_now);
                collapse_spec_valid = true;
#if !SURF_NDEBUG
                set_to_cur_wf_vertices(collapse_spec_computed_with_vertices);
#endif
            };
            return get_cached_edge_collapse();
        }

        public partial EdgePtrPair split(WavefrontEdgeList wavefront_edges);

        public bool parallel_endpoints(double time_now)
        {
            EdgeCollapseSpec e = get_edge_collapse(time_now);
            switch (e.type())
            {
                case EdgeCollapseType.UNDEFINED:
                    assert(false);
                    return false;

                case EdgeCollapseType.PAST:
                case EdgeCollapseType.FUTURE:
                    return false;

                case EdgeCollapseType.ALWAYS:
                case EdgeCollapseType.NEVER:
                    return true;
            }
            assert(false);
            return false;
        }

#if !SURF_NDEBUG

        private partial void assert_edge_sane(int collapsing_edge);

#else
private void assert_edge_sane(int collapsing_edge) {};
#endif

        private EdgeCollapseSpec get_cached_edge_collapse()
        {
            assert(collapse_spec_valid);
#if !SURF_NDEBUG
            assert_cur_wf_vertices(collapse_spec_computed_with_vertices);
#endif
            return collapse_spec;
        }

        private void invalidate_collapse_spec()
        {
            collapse_spec_valid = false;
#if !SURF_NDEBUG
            invalidte_cur_wf_vertices(collapse_spec_computed_with_vertices);
#endif
        }

#if !SURF_NDEBUG

        private void set_to_cur_wf_vertices(WavefrontVertex[] arr)
        {
            for (int i = 0; i < 2; ++i)
            {
                arr[i] = vertices[i];
            }
        }

        private void invalidte_cur_wf_vertices(WavefrontVertex[] arr)
        {
            for (int i = 0; i < 2; ++i)
            {
                arr[i] = null;
            }
        }

        private void assert_cur_wf_vertices(WavefrontVertex[] arr)
        {
            for (int i = 0; i < 2; ++i)
            {
                assert(arr[i] != null);
                assert(arr[i] == vertices[i]);
            }
        }

#endif

        private partial EdgeCollapseSpec compute_collapse(double time_now);
    }
}