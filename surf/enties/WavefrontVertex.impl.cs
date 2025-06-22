namespace SurfNet
{
    using static DebugLog;
    using static Mathex;

    public partial class WavefrontVertex
    {



        protected static partial InfiniteSpeedType get_infinite_speed_type(WavefrontEdge a, WavefrontEdge b, OrientationEnum angle)
        {
            if (a != null && b != null && angle == (int)OrientationEnum.STRAIGHT)
            {
                if (orientation(a.l().l.to_vector(), b.l().l.to_vector().perpendicular(OrientationEnum.CLOCKWISE)) == OrientationEnum.LEFT_TURN)
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





        //private static partial WavefrontVertex make_vertex(
        //  Point2 pos,
        //  double time,
        //  WavefrontEdge a,
        //  WavefrontEdge b,
        //  bool from_split
        //)
        //{

        //    //DBG_FUNC_BEGIN(DBG_KT);
        //    //DBG(DBG_KT) << "a:" << *a << " " << CGAL_line(a.l().l);
        //    //DBG(DBG_KT) << "b:" << *b << " " << CGAL_line(b.l().l);

        //    if (!from_split)
        //    {
        //        assert(a.vertex(1) != null && a.vertex(1).has_stopped());
        //        assert(b.vertex(0) != null && b.vertex(0).has_stopped());
        //    }

        //    Point2 pos_zero;
        //    LineIntersectionType lit;
        //    (lit, pos_zero) = compute_intersection(a.l().l, b.l().l);
        //    switch (lit)
        //    {
        //        case LineIntersectionType.ALL:
        //            pos_zero = pos -( (Point2)compute_velocity(Point2.ORIGIN, a.l(), b.l(), OrientationEnum.STRAIGHT))*time;
        //            break;
        //        // fall through
        //        case LineIntersectionType.ONE:
        //            pos_zero = pos - ( (Point2)compute_velocity(Point2.ORIGIN, a.l(), b.l(), OrientationEnum.STRAIGHT))*time;
        //            break;
        //        case LineIntersectionType.NONE:
        //            //  DBG(DBG_KT) << "No intersection at time 0 between supporting lines of wavefrontedges.  Parallel wavefronts crashing (or wavefronts of different speeds becoming collinear).";
        //            pos_zero = pos;
        //            break;
        //        default:
        //            // CANNOTHAPPEN_MSG << "Fell through switch which should cover all cases.";
        //            assert(false);
        //            throw new Exception("Fell through switch which should cover all cases.");

        //    }

        //    WavefrontVertex v = new WavefrontVertex(pos_zero, pos, time, a, b);
        //    assert((lit == LineIntersectionType.NONE) == (v.infinite_speed != InfiniteSpeedType.NONE));

        //    assert(v.p_at(time) == pos);
        //    //DBG_FUNC_END(DBG_KT);
        //    return v;
        //}



        string details()
        {
            //  std.ostringstream oss;
            //DEBUG_STMT(oss << "kv" << id);
            //oss << "(";
            //if (is_infinite)
            //{
            //    oss << "inf";
            //}
            //else
            //{
            //    oss << "wf: " << *incident_wavefront_edges[0]
            //        << "; " << *incident_wavefront_edges[1];
            //    oss << "; o: " << CGAL_point(pos_zero);
            //    oss << "; v" << infinite_speed;
            //    if (infinite_speed != InfiniteSpeedType.NONE)
            //    {
            //        oss << ": (" << CGAL_vector(velocity);
            //    }
            //    oss << "; s: (" << CGAL_point(pos_start) << ") @ " << SurfNet.to_double(time_start);
            //    if (has_stopped_)
            //    {
            //        oss << "; e: (" << CGAL_point(pos_stop_) << ") @ " << SurfNet.to_double(time_stop_);
            //    };
            //}
            //oss << ")";
            //return oss.str();
            return null;
        }

        //std.ostream & operator <<(std.ostream& os, const WavefrontVertex.LineIntersectionType t)
        //{
        //    switch (t)
        //    {
        //        case WavefrontVertex.LineIntersectionType.ONE: return os << "ONE";
        //        case WavefrontVertex.LineIntersectionType.ALL: return os << "ALL";
        //        case WavefrontVertex.LineIntersectionType.NONE: return os << "NONE";
        //    }
        //    CANNOTHAPPEN_MSG << "Fell through switch which should cover all cases.";
        //    assert(false);
        //    abort();
        //}

        //std.ostream &
        //operator <<(std.ostream& os, const InfiniteSpeedType &a)
        //{
        //    switch (a)
        //    {
        //        case InfiniteSpeedType.NONE:
        //            os << "";
        //            break;
        //        case InfiniteSpeedType.OPPOSING:
        //            os << "^";
        //            break;
        //        case InfiniteSpeedType.WEIGHTED:
        //            os << "~";
        //            break;
        //    }
        //    return os;
        //}


    }
}
