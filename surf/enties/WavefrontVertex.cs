using System.Diagnostics;
using System.Security.Cryptography.X509Certificates;
using System.Text;

namespace SurfNet
{
    using static DebugLog;
    using static Mathex;
    using VertexAngle = OrientationEnum;

    public enum InfiniteSpeedType
    { NONE, OPPOSING, WEIGHTED };

    public partial class WavefrontVertex
    {

        public int dbg_orientation = 0;
        public Vector2 dbg_vector_velocity = new Vector2(0, 0);
        public InfiniteSpeedType infinite_speed;
        public bool is_beveling;
        public bool is_infinite;
        public bool is_initial;
        public Point2 pos_zero;
        // This wavefront vertex is either between parallel, opposing wavefront elements that have crashed
        // into each other and become collinear, or it is between neighboring
        // wavefront edges that have become collinear yet have different weights.
        public Vector2 velocity;

        private SkeletonDCELHalfedge _halfedge = null;
        private VertexAngle angle;
        private bool has_stopped_ = false;
        // wavefront vertices form a doubly-linked edge list to represent
        // the boundary of their left(0) and right(1) incident faces.
        //
        // prev points to the wavefront vertex earlier in time, next to
        // the one later in time, so when traversing a face, care needs
        // to be taken at each arc (i.e. wavefront-vertex) wrt direction.
        private WavefrontEdge?[] incident_wavefront_edges;

        private bool is_degenerate_ = false;
        private WavefrontVertex?[] next_vertex_ = new WavefrontVertex[2];
        private Point2 pos_stop_;
        private WavefrontVertex?[] prev_vertex_ = new WavefrontVertex[2];
        private Polynomial1D px_, py_;
        // if pos_stop == pos_start 
        private SkeletonDCELHalfedge?[] skeleton_dcel_halfedge_;

        private double time_stop_;
        public SkeletonDCELHalfedge Halfedge
        {
            get => _halfedge != null ? _halfedge : Vertex.Halfedge;
            set
            {
                _halfedge = value;
                if (Vertex.Halfedge == null)
                { Vertex.Halfedge = value; }
            }
        }

        public int Id { get; private set; }
        public WavefrontVertex NextInLAV { get; set; }
        public Point2 pos_start => Vertex.Point;
        public WavefrontVertex PrevInLAV { get; set; }
        public double time_start => Vertex.Time;
        public SkeletonDCELVertex Vertex { get; private set; }
        public WavefrontVertex(
                        int id,
                    Point2 p_pos_zero,
                    SkeletonDCELVertex vertex,
                    WavefrontEdge? a = null,
                    WavefrontEdge? b = null,
                    bool p_is_initial = false,
                    bool p_is_beveling = false,
                    bool p_is_infinite = false)

        {
            Id = id;

            pos_zero = (p_pos_zero);
            Vertex = vertex;

            incident_wavefront_edges = new WavefrontEdge?[] { a, b };
            angle = (a != null && b != null) ? orientation(a.l().l.ToVector(), b.l().l.ToVector()) : OrientationEnum.STRAIGHT;
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

        protected enum LineIntersectionType
        {
            ONE,   /* the two lines intersect in one point */
            ALL,   /* the two lines are parallel and coincide */
            NONE,  /* the two liens are parallel and distinct */
        };

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
                    System.Diagnostics.Debug.Assert(false);
                    throw new Exception("No point intersection between WavefrontEmittingEdges at offset 1.  Bad.");
                }
                var intersection = lit.Points[0];
                result = new Vector2(pos_zero, intersection);
            }
            else
            {
                //DBG(DBG_KT) << "a:" << CGAL_vector(a.normal);
                //DBG(DBG_KT) << "b:" << CGAL_vector(b.normal);

                if (orientation(a.l.ToVector(), b.l.ToVector().perpendicular(OrientationEnum.CLOCKWISE)) == OrientationEnum.RIGHT_TURN)
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

        public string Debug()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine($" kv:{Id:##0} start:{pos_start.Debug()} t:{time_start.Debug(4)} stop: {pos_stop_} t:{time_stop_} ");
            sb.AppendLine($"    iwe0: {incident_wavefront_edges[0]} iwe1: {incident_wavefront_edges[1]}");
            //sb.AppendLine($"    he[0]:   {skeleton_dcel_halfedge_[0]}");
            //sb.AppendLine($"    he[1]:   {skeleton_dcel_halfedge_[1]}");
            sb.AppendLine($"   f0 prev: {prev_vertex_[0]} next: {next_vertex_[0]}");
            sb.AppendLine($"   f1 prev: {prev_vertex_[1]} next: {next_vertex_[1]}");
            return sb.ToString();


            return null;
        }

        // public bool is_straight() const { return angle == STRAIGHT; }
        public bool has_stopped()
        { return has_stopped_; }

        public WavefrontEdge incident_wavefront_edge(int i)
        {
            assert(i <= 1);
            return incident_wavefront_edges[i];
        }

        public bool is_convex_or_straight()
        { return angle != OrientationEnum.REFLEX; }

        public bool is_reflex_or_straight()
        { return angle != OrientationEnum.CONVEX; }

        public WavefrontVertex next_vertex(int side)
        {
            assert(side <= 1);
            return next_vertex_[side];
        }

        public Point2 p_at(double t)
        {
            assert(!has_stopped_ || t <= time_stop_);
            assert(!is_infinite);
            return (Point2)(pos_zero + (Point2)(velocity * t));
        }

        public Point2 p_at_dbg(double t)
        {
            if (has_stopped() && t >= time_stop_)
                return pos_stop_;
            return (Point2)(pos_zero + (Point2)(velocity * t));
        }

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

        public Point2 pos_stop()
        { return pos_stop_; }

        public WavefrontVertex prev_vertex(int side)
        {
            assert(side <= 1);
            return prev_vertex_[side];
        }

        public Polynomial1D px()
        {
            assert(!is_infinite);
            assert(!has_stopped_);
            return px_;
        }

        //  std.string Debug() const;
        public Polynomial1D py()
        {
            assert(!is_infinite);
            assert(!has_stopped_);
            return py_;
        }

        public void set_incident_wavefront_edge(int i, WavefrontEdge e)
        {
            assert(i <= 1);
            assert(incident_wavefront_edges[i] != null);
            assert(e != null);
            assert(incident_wavefront_edges[i].l() == e.l());
            incident_wavefront_edges[i] = e;
        }

        public SkeletonDCELHalfedge skeleton_dcel_halfedge(int i)
        {
            assert(i <= 1);

            assert(!is_degenerate());
            return skeleton_dcel_halfedge_[i];
        }

        public void stop(double t)
        {
            assert(!has_stopped_);
            assert(infinite_speed == InfiniteSpeedType.NONE);
            time_stop_ = t;
            pos_stop_ = p_at(t);
            has_stopped_ = true;

            if (pos_stop_.AreNear(pos_start))
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

        public double time_stop()
        { return time_stop_; }

        public override string ToString()
        {
            return $"wv:{Id} ";
        }
        //public WavefrontVertex(int idx, Point2 point_2) : this(idx, point_2, point_2, 0, null, null)
        //{
        //    this.Id = idx;
        //}

        /** type of intersections for two lines */
        //  friend std.ostream& operator<<(std.ostream& os, const WavefrontVertex.LineIntersectionType t);

        //KineticTriangle const * const * triangles() const { return incident_triangles; };
        public WavefrontEdge[] wavefronts()
        { return incident_wavefront_edges; }

        [Conditional("DEBUG")]
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

        internal bool is_degenerate()
        { return is_degenerate_; }

        internal void link_tail_to_tail(WavefrontVertex other)
        {
            assert(prev_vertex_[0] == null);
            assert(other.prev_vertex_[1] == null);
            prev_vertex_[0] = other;
            other.prev_vertex_[1] = this;
            //DBG(DBG_KT) << " For " << this << " prev_vertex_[0] := " << other;
            //DBG(DBG_KT) << " For " << other << " prev_vertex_[1] := " << this;
        }

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

        /// <summary>
        /// join two wavefront vertices, tail-to-tail.
        ///
        /// This is used after a split event.  Is called at a, where a's left (ccw)
        /// side is towards the split edge, i.e. where prev_vertex[0] is null.
        ///
        /// This is also used to link initial vertices while beveling.  Going
        /// clockwise about a vertex, this is called at each kinetc vertex with
        /// the previous one as an argument.
        ///
        /// </summary>
        internal void set_skeleton_dcel_halfedge(int i, SkeletonDCELHalfedge he)
        {
            assert(i <= 1);

            assert(!is_degenerate());
            skeleton_dcel_halfedge_[i] = he;
        }

        protected static Intersection compute_intersection(Line2 a, Line2 b)
        {
            return Mathex.intersection(a, b);
        }

        protected static partial InfiniteSpeedType get_infinite_speed_type(WavefrontEdge a, WavefrontEdge b, OrientationEnum angle);
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
        /*
        public
            friend inline std.ostream& operator<<(std.ostream& os, const WavefrontVertex * const kv) {
              if (kv) {
                os << "kv";
        DEBUG_STMT(os << kv.Id);
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
        // ==================== functions maintaining the DCEL =====================

        /** set the successor in the DCEL
         *
         * Also update their prev (or next) pointer depending on whether we have
         * head_to_tail set to true or not.
         */
        protected static partial InfiniteSpeedType get_infinite_speed_type(WavefrontEdge a, WavefrontEdge b, OrientationEnum angle)
        {
            if (a != null && b != null && angle == (int)OrientationEnum.STRAIGHT)
            {
                if (orientation(a.l().l.ToVector(), b.l().l.ToVector().perpendicular(OrientationEnum.CLOCKWISE)) == OrientationEnum.LEFT_TURN)
                {
                    return InfiniteSpeedType.OPPOSING;
                }
                else if (a.l().weight != b.l().weight)
                {
                    return InfiniteSpeedType.WEIGHTED;
                }
                else
                {
                    assert(a.l().normal == b.l().normal);
                    return InfiniteSpeedType.NONE;
                }
            }
            else
            {
                return InfiniteSpeedType.NONE;
            }
        }

        
    }
}