namespace SurfNet
{
    using static DebugLog;
    using static Mathex;
    public class WavefrontVertexList : List<WavefrontVertex>
    {




        public WavefrontVertex make_initial_vertex(
           Point2 pos_zero,
           WavefrontEdge a,
           WavefrontEdge b,
           bool is_beveling =false
         )
        {
            assert(a.l().l.has_on(pos_zero));
            assert(b.l().l.has_on(pos_zero));

            var v = new WavefrontVertex(this.Count, pos_zero, pos_zero, CORE_ZERO, a, b, true, is_beveling);
            Add(v);
            return v;
        }
        public static WavefrontVertex make_infinite_vertex()
        {
            return new WavefrontVertex(-1, Point2.ORIGIN, Point2.ORIGIN, CORE_ZERO, null, null, true, false, true);
        }



        /** create a new kinetic vertex during the propagation period.
         *  link the vertex to its ancestors on the wavefront sides
         *  that are still propagating.
         *
         *  Linking the predecessors on the no-longer-propagating side
         *  is the responsibility of the caller, as it's not clear
         *  how they should be linked.  (together in the case of
         *  a constraint collapse, but not necessarily in spoke collapses.)
         */
        public WavefrontVertex make_vertex(
          Point2 pos,
          double time,
          WavefrontEdge a,
          WavefrontEdge b,
      bool from_split = false
    )
        {

            if (!from_split)
            {
                assert(a.vertex(1)!=null && a.vertex(1).has_stopped() );
                assert(b.vertex(0) != null && b.vertex(0).has_stopped() );
            }

            Line2 la = a.l().l;
            Line2 lb = b.l().l;

            var lit = Mathex.intersection(la, lb);


            Point2 pos_zero;
            switch (lit.Result)
            {
                case Intersection.Intersection_results.LINE:
                    pos_zero = pos - WavefrontVertex.compute_velocity(Point2.ORIGIN, a.l(), b.l(), OrientationEnum.STRAIGHT) * time;
                    break;
                // fall through
                case Intersection.Intersection_results.POINT:
                    pos_zero = lit.Points[0];
                    break;

                case Intersection.Intersection_results.NO_INTERSECTION:
                    //DBG(DBG_KT) << "No intersection at time 0 between supporting lines of wavefrontedges.  Parallel wavefronts crashing (or wavefronts of different speeds becoming collinear).";
                    pos_zero = pos;
                    break;
                default:
                    //CANNOTHAPPEN_MSG << "Fell through switch which should cover all cases.";
                    assert(false);

                    throw new Exception("Fell through switch which should cover all cases.");

            }

            WavefrontVertex v = new WavefrontVertex(this.Count, pos_zero, pos, time, a, b);
            Add(v);
            assert((lit.Result == Intersection.Intersection_results.NO_INTERSECTION) == (v.infinite_speed != InfiniteSpeedType.NONE));

            assert(v.p_at(time).AreNear(pos));
            //   DBG_FUNC_END(DBG_KT);
            return v;




           
        }
    }
}