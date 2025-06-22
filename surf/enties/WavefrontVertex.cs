using System.Diagnostics;

namespace SurfNet
{
    using static DebugLog;
    using static Mathex;
    using VertexAngle = OrientationEnum;
    public enum InfiniteSpeedType { NONE, OPPOSING, WEIGHTED };


    public partial class WavefrontVertex
    {

        

        public Point2 pos_zero;
        public Point2 pos_start;
        public WavefrontVertex NextInLAV { get; set; }
        public WavefrontVertex PrevInLAV { get; set; }

        public double time_start;
        private WavefrontEdge?[] incident_wavefront_edges;
        VertexAngle angle;
        public bool is_initial;
        public bool is_beveling;
        public bool is_infinite;
        public InfiniteSpeedType infinite_speed; /** This wavefront vertex is
      either between parallel, opposing wavefront elements that have crashed
      into each other and become collinear, or it is between neighboring
      wavefront edges that have become collinear yet have different weights. */
        public Vector2 velocity;
        private Polynomial1D px_, py_;
        private bool has_stopped_ = false;
        private double time_stop_;
        private Point2 pos_stop_;

        private bool is_degenerate_ = false; /* if pos_stop == pos_start */
        private SkeletonDCELHalfedge?[] skeleton_dcel_halfedge_;

        /* wavefront vertices form a doubly-linked edge list to represent
         * the boundary of their left(0) and right(1) incident faces.
         *
         * prev points to the wavefront vertex earlier in time, next to
         * the one later in time, so when traversing a face, care needs
         * to be taken at each arc (i.e. wavefront-vertex) wrt direction.
         */
        private WavefrontVertex?[] next_vertex_ = new WavefrontVertex[2];
        private WavefrontVertex?[] prev_vertex_ = new WavefrontVertex[2];

        

        public int Id { get; private set; }

        private Point2 point_2;
        public override string ToString()
        {
            return $"wv:{Id} ";
        }
        public WavefrontVertex(
                        int id,
                    Point2 p_pos_zero,
                    Point2 p_pos_start,
                    double p_time_start,
                    WavefrontEdge? a = null,
                    WavefrontEdge? b = null,
                    bool p_is_initial = false,
                    bool p_is_beveling = false,
                    bool p_is_infinite = false)

        {


            Id = id;

            pos_zero = (p_pos_zero);
            pos_start = (p_pos_start);
            time_start = (p_time_start);
            incident_wavefront_edges = new WavefrontEdge?[] { a, b };
            angle = (a != null && b != null) ? orientation(a.l().l.to_vector(), b.l().l.to_vector()) : OrientationEnum.STRAIGHT;
            is_initial = (p_is_initial);
            is_beveling = (p_is_beveling);
            is_infinite = (p_is_infinite);
            infinite_speed = get_infinite_speed_type(a, b, angle);
            velocity = ((a != null && b != null && (infinite_speed == InfiniteSpeedType.NONE)) ? compute_velocity(pos_zero, a.l(), b.l(), angle) : Vector2.NaN);
            px_ = ((infinite_speed != InfiniteSpeedType.NONE) ? new Polynomial1D(0) : new Polynomial1D(pos_zero.X, velocity.X));
            py_ = ((infinite_speed != InfiniteSpeedType.NONE) ? new Polynomial1D(0) : new Polynomial1D(pos_zero.Y, velocity.Y));
            skeleton_dcel_halfedge_ = new SkeletonDCELHalfedge?[] { null, null };
            next_vertex_ = new WavefrontVertex?[] { null, null };
            prev_vertex_ = new WavefrontVertex?[] { null, null };


            // assert(!!a == !!b);
        }

        public WavefrontVertex(int idx, Point2 point_2) : this(idx, point_2, point_2, 0, null, null)
        {
            this.Id = idx;

        }


        /** type of intersections for two lines */
        protected enum LineIntersectionType
        {
            ONE,   /* the two lines intersect in one point */
            ALL,   /* the two lines are parallel and coincide */
            NONE,  /* the two liens are parallel and distinct */
        };
        //  friend std.ostream& operator<<(std.ostream& os, const WavefrontVertex.LineIntersectionType t);

        protected static partial InfiniteSpeedType get_infinite_speed_type(WavefrontEdge a, WavefrontEdge b, OrientationEnum angle);
        protected static Intersection compute_intersection(Line2 a, Line2 b)
        {


            return Mathex.intersection(a, b);
        }

        public static Vector2 compute_velocity(
             Point2 pos_zero,
             WavefrontSupportingLine a,
             WavefrontSupportingLine b,
             OrientationEnum angle)
        {
            Vector2 result;

            if (angle != OrientationEnum.STRAIGHT)
            {
                Line2 la = a.line_at_one();
                Line2 lb = b.line_at_one();

                var lit = compute_intersection(la, lb);

                if (lit.Result != Intersection.Intersection_results.POINT)
                {
                    // CANNOTHAPPEN_MSG << "No point intersection between WavefrontEmittingEdges at offset 1.  Bad.";
                    Debug.Assert(false);
                    throw new Exception("No point intersection between WavefrontEmittingEdges at offset 1.  Bad.");
                }
                var intersection = lit.Points[0];
                result = new Vector2(pos_zero, intersection);
            }
            else
            {
                //DBG(DBG_KT) << "a:" << CGAL_vector(a.normal);
                //DBG(DBG_KT) << "b:" << CGAL_vector(b.normal);

                if (orientation(a.l.to_vector(), b.l.to_vector().perpendicular(OrientationEnum.CLOCKWISE)) == OrientationEnum.RIGHT_TURN)
                {
                    /* They are in the same direction */
                    if (a.normal == b.normal)
                    {
                        result = a.normal;
                    }
                    else
                    {
                        throw new Exception("collinear incident wavefront edges with different speeds.");

                    }
                }
                else
                {
                    throw new Exception("This is an infinitely fast vertex.  We should not be in compute_velocity().");
                }
            }
            return result;
        }

        public bool is_reflex_or_straight() { return angle != OrientationEnum.CONVEX; }
        public bool is_convex_or_straight() { return angle != OrientationEnum.REFLEX; }
        // public bool is_straight() const { return angle == STRAIGHT; }
        public bool has_stopped() { return has_stopped_; }
        public double time_stop() { return time_stop_; }
        public Point2 pos_stop() { return pos_stop_; }

        //KineticTriangle const * const * triangles() const { return incident_triangles; };
        public WavefrontEdge[] wavefronts() { return incident_wavefront_edges; }


        public Point2 p_at_dbg(double t)
        {
            if(has_stopped() &&  t >= time_stop_) return pos_stop_;
             return (Point2)(pos_zero + (Point2)(velocity * t));
        }

        public Point2 p_at(double t)
        {
            assert(!has_stopped_ || t <= time_stop_);
            assert(!is_infinite);
            return (Point2)(pos_zero + (Point2)(velocity * t));
        }

        /** return the position of this vertex for drawing purposes.
         *
         * If the time to draw is later than the stop position, return
         * the stop position.
         *
         * If the time to draw is prior to the start position, no such
         * special handling is done and we return the location where
         * the vertex would have been such that it is at the start position
         * at the start time given its velocity.
         */
        public Point2 p_at_draw(double t)
        {
            assert(!is_infinite);
            bool return_stop_pos;

            if (has_stopped_)
            {
                if (t < time_stop_)
                {
                    return_stop_pos = false;
                }
                else
                {
                    return_stop_pos = true;
                };
            }
            else
            {
                return_stop_pos = false;
            }

            return return_stop_pos ? pos_stop_
                                   : (pos_zero + (Point2)velocity * t);
        }

        public void stop(double t)
        {
            assert(!has_stopped_);
            assert(infinite_speed == InfiniteSpeedType.NONE);
            time_stop_ = t;
            pos_stop_ = p_at(t);
            has_stopped_ = true;

            if (pos_stop_ .AreNear(pos_start))
            {
                is_degenerate_ = true;
                assert(time_stop_.AreNear(time_start));
            };
        }

        public void stop(double t, Point2 p)
        {
            assert(!has_stopped_);
            assert(infinite_speed != InfiniteSpeedType.NONE);
            time_stop_ = t;
            pos_stop_ = p;
            has_stopped_ = true;

            assert(time_stop_.AreNear(time_start));

            if (pos_stop_.AreNear(pos_start))
            {
                is_degenerate_ = true;
            };
        }

        public WavefrontEdge incident_wavefront_edge(int i)
        {
            assert(i <= 1);
            return incident_wavefront_edges[i];
        }
        public void set_incident_wavefront_edge(int i, WavefrontEdge e)
        {
            assert(i <= 1);
            assert(incident_wavefront_edges[i] != null);
            assert(e != null);
            assert(incident_wavefront_edges[i].l() == e.l());
            incident_wavefront_edges[i] = e;
        }

        public WavefrontVertex next_vertex(int side)
        {
            assert(side <= 1);
            return next_vertex_[side];
        }
        public WavefrontVertex prev_vertex(int side)
        {
            assert(side <= 1);
            return prev_vertex_[side];
        }

        /*
        public
            friend inline std.ostream& operator<<(std.ostream& os, const WavefrontVertex * const kv) {
              if (kv) {
                os << "kv";
        DEBUG_STMT(os << kv.id);
        os << (kv.angle == CONVEX ? "c" :
               kv.angle == REFLEX ? "r" :
               kv.angle == STRAIGHT ? "=" :
                                       "XXX-INVALID-ANGLE")
           << kv.infinite_speed
           << (kv.has_stopped_ ? "s" : "");
              } else
        {
            os << "kv*";
        }
        return os;
            }
        */
        //  std.string details() const;

        public Polynomial1D px()
        {
            assert(!is_infinite);
            assert(!has_stopped_);
            return px_;
        }
        public Polynomial1D py()
        {
            assert(!is_infinite);
            assert(!has_stopped_);
            return py_;
        }

#if !SURF_NDEBUG
        internal void assert_valid()
        {
            assert(is_initial || !is_beveling); // !initial => !beveling   <=>  !!initial v !beveling
            for (int i = 0; i < 2; ++i)
            {
                if (is_initial)
                {
                    assert(prev_vertex_[i] == null || is_beveling);
                }
                else
                {
                    assert(prev_vertex_[i] != null);
                }
                assert(!has_stopped() ^ next_vertex_[i] != null);
            }
        }
#else
void assert_valid() const {};
#endif

        // ==================== functions maintaining the DCEL =====================

        /** set the successor in the DCEL
         *
         * Also update their prev (or next) pointer depending on whether we have
         * head_to_tail set to true or not.
         */
        internal void set_next_vertex(int side, WavefrontVertex next, bool head_to_tail = true)
        {
            assert(side <= 1);
            assert(next != null);
            assert(has_stopped());
            //DBG(DBG_KT) << " for " << this << " next_vertex_[" << side << "] is " << next_vertex_[side];
            //DBG(DBG_KT) << " for " << this << " next_vertex_[" << side << "] := " << next;
            assert(next_vertex_[side] == null);
            next_vertex_[side] = next;

            if (head_to_tail)
            {
                //DBG(DBG_KT) << " +for " << next << " prev_vertex_[" << side << "] is " << next.prev_vertex_[side];
                //DBG(DBG_KT) << " +for " << next << " prev_vertex_[" << side << "] := " << this;
                assert(next.prev_vertex_[side] == null);
                next.prev_vertex_[side] = this;
            }
            else
            {
                /* head to head */
                //DBG(DBG_KT) << " +for " << next << " next_vertex_[" << 1 - side << "] is " << next.next_vertex_[1 - side];
                //DBG(DBG_KT) << " +for " << next << " next_vertex_[" << 1 - side << "] := " << this;
                assert(next.next_vertex_[1 - side] == null);
                next.next_vertex_[1 - side] = this;
            };
        }

        /** join two wavefront vertices, tail-to-tail.
         *
         * This is used after a split event.  Is called at a, where a's left (ccw)
         * side is towards the split edge, i.e. where prev_vertex[0] is null.
         *
         * This is also used to link initial vertices while beveling.  Going
         * clockwise about a vertex, this is called at each kinetc vertex with
         * the previous one as an argument.
         */
        internal void link_tail_to_tail(WavefrontVertex other)
        {
            assert(prev_vertex_[0] == null);
            assert(other.prev_vertex_[1] == null);
            prev_vertex_[0] = other;
            other.prev_vertex_[1] = this;
            //DBG(DBG_KT) << " For " << this << " prev_vertex_[0] := " << other;
            //DBG(DBG_KT) << " For " << other << " prev_vertex_[1] := " << this;
        }

        internal bool is_degenerate() { return is_degenerate_; }

     public   SkeletonDCELHalfedge skeleton_dcel_halfedge(int i)
        {
            assert(i <= 1);

            assert(!is_degenerate());
            return skeleton_dcel_halfedge_[i];
        }

        internal void set_skeleton_dcel_halfedge(int i, SkeletonDCELHalfedge he)
        {
            assert(i <= 1);

            assert(!is_degenerate());
            skeleton_dcel_halfedge_[i] = he;
        }
    }
}