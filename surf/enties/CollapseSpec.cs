namespace SurfNet
{
    using static DebugLog;
    using static Mathex;
    public enum EdgeCollapseType
    {
        /// <summary>
        /// Undefined
        /// </summary>
        UNDEFINED = 1,
        /// <summary>
        /// endpoints moving away from another
        /// </summary>
        PAST,
        /// <summary>
        /// endpoints moving towards one another
        /// </summary>
        FUTURE,
        /// <summary>
        /// endpoints moving moving in parrallel and conincident
        /// </summary>
        ALWAYS,
        /// <summary>
        /// endpoints moving moving in parrallel but not conincident
        /// </summary>
        NEVER,
    };

    /*
std.ostream & operator <<(std.ostream& os, const EdgeCollapseType a);
    */





    public enum CollapseType
    {
        UNDEFINED = 1,
        /// <summary>
        /// This triangle has a vertex
        /// which is between parallel, opposing wavefront elements
        /// that have crashed into each other and their intersection
        /// is now a line segment.
        /// </summary>
        FACE_HAS_INFINITELY_FAST_VERTEX_OPPOSING,
        /// <summary>
        /// </summary>
        TRIANGLE_COLLAPSE,
        // NEIGHBORING_TRIANGLE_COLLAPSES,  /* do we ever need to handle this */
        /// <summary>
        /// </summary>
        CONSTRAINT_COLLAPSE,
        /// <summary>
        /// two non-incident vertices become incident,
        /// splitting the wavefront here.
        ///  UNUSED except in get_generic
        /// </summary>
        SPOKE_COLLAPSE,
        /// <summary>
        /// vertex moves onto supporting line of constraint,
        /// can refine event type when it comes to it.
        /// </summary>
        SPLIT_OR_FLIP_REFINE,
        /// <summary>
        /// This triangle has a vertex which is
        /// between parallel adjacent wavefront elements that have
        /// different weights but move in the same direction.
        /// </summary>
        FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED,
        /// <summary>
        /// vertex moves into spoke (triangulation edge interior), flip event
        /// </summary>
        VERTEX_MOVES_OVER_SPOKE,
        /// <summary>
        /// the ccw vertex of the infinite vertex in an
        /// unbounded triangle leaves the convex hull of the
        /// wavefront polygon

        /// </summary>
        CCW_VERTEX_LEAVES_CH,
        // GENERIC_FLIP_EVENT,
        /// <summary>
        /// The triangle will collapse at this time, but we should
        /// never see this as prior events should have rebuilt the
        /// triangulation in some way.  If this is the next event,
        /// something went wrong.
        /// </summary>
        INVALID_EVENT,
        /// <summary>
        /// Leave this one last.  It serves also as a counter!
        /// </summary>
        NEVER,

        RefineTriangle
    };
    /*
std.ostream & operator <<(std.ostream& os, const CollapseType a);
    */

    public class CollapseSpec : IComparable<CollapseSpec>
    {
        public static int COUNTER_NT_cmp;

        private CollapseType type_;
        protected double time_;

        public override string ToString()
        {
            return $"{type_} time:{time_.Debug(4)} e:{(relevant_edge_ != -1 ? relevant_edge_.ToString():"")} sk:{(secondary_key_ != double.NegativeInfinity? secondary_key_.Debug(4):"")} ";
        }

       

        /* extra info */
        // for all collapses listed in requires_relevant_edge(), such as CONSTRAINT_COLLAPSE
        private int relevant_edge_ = -1;
        // for VERTEX_MOVES_OVER_SPOKE
        private double secondary_key_ = double.NegativeInfinity; // Higher number is more important
                                       // - for flip events: longest spoke: squared length
                                       // - for weighted infinite one: higher edge speed (if faster edge wins)



        public void update_collapse(CollapseSpec other)
        {
            type_ = other.type_;
            time_ = other.time_;
            relevant_edge_ = other.relevant_edge_;
            secondary_key_ = other.secondary_key_;
        }
        private static bool requires_relevant_edge(CollapseType type)
        {
            return
               type == CollapseType.FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED ||
               type == CollapseType.CONSTRAINT_COLLAPSE ||
               type == CollapseType.SPOKE_COLLAPSE ||
               type == CollapseType.SPLIT_OR_FLIP_REFINE ||
               type == CollapseType.VERTEX_MOVES_OVER_SPOKE ||
               type == CollapseType.CCW_VERTEX_LEAVES_CH ||
               false;
        }
        public static bool requires_relevant_edge_plus_secondary_key(CollapseType type)
        {
            return
               type == CollapseType.FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED ||
               type == CollapseType.VERTEX_MOVES_OVER_SPOKE ||
               false;
        }
        public readonly int component;

        public bool requires_relevant_edge()
        {
            return requires_relevant_edge(type_);
        }
        public bool requires_relevant_edge_plus_secondary_key()
        {
            return requires_relevant_edge_plus_secondary_key(type_);
        }

        /* TODO */
        /*
      public CollapseSpec(const CollapseSpec&) = default;
    public CollapseSpec(CollapseSpec &&) = default;

    CollapseSpec & operator = (CollapseSpec&& o)
    {
        assert(component == o.component);

        type_ = std.move(o.type_);
        time_ = std.move(o.time_);
        relevant_edge_ = std.move(o.relevant_edge_);
        secondary_key_ = std.move(o.secondary_key_);
        return *this;
    }

    CollapseSpec & operator =(const CollapseSpec& o)
    {
        assert(component == o.component);
        type_ = o.type_;
        time_ = o.time_;
        relevant_edge_ = o.relevant_edge_;
        secondary_key_ = o.secondary_key_;
        return *this;
    }
        */

       
        protected CollapseSpec(CollapseSpec other)
        {
            type_ = other.type_;
            time_ = other.time_;
            relevant_edge_ = other.relevant_edge_;
            secondary_key_ = other.secondary_key_;
        }
        public CollapseSpec(int p_component, CollapseType type = CollapseType.UNDEFINED)

        {
            type_ = type;
            component = p_component;
            //    assert(type_ == CollapseType.NEVER);
        }
        public CollapseSpec(int p_component, CollapseType type, double time)
        {
            type_ = type;
            time_ = time;
            component = p_component;

            assert(type_ != CollapseType.UNDEFINED);
            assert(type_ != CollapseType.NEVER);
            assert(!requires_relevant_edge(type_));
        }
        public CollapseSpec(int p_component, CollapseType type, double time, int relevant_edge)

        {
            type_ = type;
            time_ = time;
            relevant_edge_ = relevant_edge;

            component = p_component;

            assert(requires_relevant_edge(type_));
            assert(!requires_relevant_edge_plus_secondary_key(type_));
            assert(0 <= relevant_edge_ && relevant_edge_ < 3);
        }

        public CollapseSpec(int p_component,
                        CollapseType type,
                        double time,
                        int relevant_edge,
                        double secondary_key)
        {
            type_ = type;
            time_ = time;
            relevant_edge_ = relevant_edge;
            secondary_key_ = secondary_key;

            component = p_component;

            assert(requires_relevant_edge(type_));
            assert(requires_relevant_edge_plus_secondary_key(type_));
            assert(0 <= relevant_edge_ && relevant_edge_ < 3);
        }
        public CollapseSpec(int p_component,
                         EdgeCollapseSpec edge_collapse,
                         int relevant_edge)

        {
            type_ = edge_collapse.type() == EdgeCollapseType.FUTURE ? CollapseType.CONSTRAINT_COLLAPSE :
                    edge_collapse.type() == EdgeCollapseType.ALWAYS ? CollapseType.CONSTRAINT_COLLAPSE :
                                                                     CollapseType.NEVER;
            time_ = type_ == CollapseType.CONSTRAINT_COLLAPSE ? edge_collapse.time() : CORE_ZERO;
            if (type_ == CollapseType.CONSTRAINT_COLLAPSE)
                assert(true);
            relevant_edge_ = type_ == CollapseType.CONSTRAINT_COLLAPSE ? relevant_edge : 0;
            component = (p_component);

            assert(edge_collapse.type() == EdgeCollapseType.FUTURE ||
                   edge_collapse.type() == EdgeCollapseType.ALWAYS ||
                   edge_collapse.type() == EdgeCollapseType.NEVER ||
                   edge_collapse.type() == EdgeCollapseType.PAST);
            assert(0 <= relevant_edge_ && relevant_edge_ < 3);
        }
        /*
        CollapseSpec(const CollapseSpec &o)
          : type_(o.type_)
          , time_(o.time_)
        {}
        */

        public CollapseType type() { return type_; }
        public double time() { return time_; }
        public double get_printable_time() { return time_; }
        public double get_printable_secondary_key() { return secondary_key_; }
        public int relevant_edge()
        {
            assert(requires_relevant_edge());
            assert(0 <= relevant_edge_ && relevant_edge_ < 3);
            return relevant_edge_;
        }
        public bool allows_refinement_to(CollapseSpec o)
        {
            assert(time_ == o.time_);
            if (type_ == CollapseType.SPLIT_OR_FLIP_REFINE)
            {
                if (o.type_ == CollapseType.VERTEX_MOVES_OVER_SPOKE ||
                    o.type_ == CollapseType.SPOKE_COLLAPSE)
                {
                    if (relevant_edge_ != o.relevant_edge_)
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        private static int compare_NT(double a, double b)
        {
            ++COUNTER_NT_cmp;
            if (a.AreNear(b)) return EQUAL;
        
            if (a < b)
            {
                return SMALLER;
            }
            else
            {
                return LARGER;
            }
           
        }
        public override bool Equals(object? obj)
        {
            bool result = false;
            if (obj is CollapseSpec c)
                result = CompareTo(c) == EQUAL;
            return result;
        }


        public static bool operator ==(CollapseSpec? a, CollapseSpec? b)
        {
            if (((object)a) == null)
            {
                if (((object)b) == null) return true;
            }
            if (b == null) return false;

            return a.Equals(b);
          
        }
        public static bool operator !=(CollapseSpec? a, CollapseSpec? b)
        {
            return !(a==b);
        }

        public int CompareTo(CollapseSpec? o)
        {
            assert(type_ != CollapseType.UNDEFINED);
            assert(o.type_ != CollapseType.UNDEFINED);

            if (type_ == CollapseType.NEVER)
            {
                if (o.type_ == CollapseType.NEVER)
                {
                    return EQUAL;
                }
                else
                {
                    return LARGER;
                }
            }
            else if (o.type_ == CollapseType.NEVER)
            {
                return SMALLER;
            }

            //if (component < o.component)
            //{
            //    return SMALLER;
            //}
            //else if (component > o.component)
            //{
            //    return LARGER;
            //}

            var c = compare_NT(time_, o.time_);
            if (c == EQUAL)
            {
                if (type_ < o.type_)
                {
                    c = SMALLER;
                }
                else if (type_ > o.type_)
                {
                    c = LARGER;
                }
                else if (requires_relevant_edge_plus_secondary_key())
                {
                    c = -compare_NT(secondary_key_, o.secondary_key_);
                }
            }
            return c;
        }

        /*
        public bool operator <(CollapseSpec o) { return compare(o) == SMALLER; }
        public bool operator >(CollapseSpec o) const { return compare(o) == LARGER; }
    public bool operator >=(CollapseSpec o) const { return !(*this < o); };
    public bool operator <=(CollapseSpec o) { return !(*this > o); };
    public bool operator ==( CollapseSpec o) { return compare(o) == SurfNet.EQUAL; }
    public bool operator !=( CollapseSpec o) { return !(*this == o); };
        */
    };
    // std.ostream & operator <<(std.ostream& os, const CollapseSpec& s);
}