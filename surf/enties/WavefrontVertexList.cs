namespace SurfNet
{
    using static DebugLog;
    using static Mathex;
    public class WavefrontVertexList : List<WavefrontVertex>
    {




        public WavefrontVertex make_initial_vertex(SkeletonDCEL skeleton, int vId,
           Point2 pos_zero,
           WavefrontEdge a,
           WavefrontEdge b,
           bool is_beveling =false
         )
        {
            assert(a.l().l.has_on(pos_zero));
            assert(b.l().l.has_on(pos_zero));

            
            var v = new WavefrontVertex(vId, pos_zero, skeleton.vertices[vId], a, b, true, is_beveling);


            this[vId]=v;
            return v;
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
            SkeletonDCEL skeleton,
          Point2 pos,
          double time,
          WavefrontEdge a,
          WavefrontEdge b,
      bool from_split = false
    )
        {
          var v =  make_vertex(skeleton.new_vertex(pos, time), a, b, from_split);
            v.Vertex.Id = v.Id;

            return v;
           


           
        }

        public WavefrontVertex make_vertex( SkeletonDCELVertex vertex,
         WavefrontEdge a,
         WavefrontEdge b,
     bool from_split = false
   )
        {
          
            if (!from_split)
            {
                assert(a.vertex(1) != null && a.vertex(1).has_stopped());
                assert(b.vertex(0) != null && b.vertex(0).has_stopped());
            }

            Line2 la = a.l().l;
            Line2 lb = b.l().l;

            var lit = Mathex.intersection(la, lb);


            Point2 pos_zero;
            switch (lit.Result)
            {
                case Intersection.Intersection_results.LINE:
                    pos_zero = vertex.Point - WavefrontVertex.compute_velocity(Point2.ORIGIN, a.l(), b.l(), OrientationEnum.STRAIGHT) * vertex.Time;
                    break;
                // fall through
                case Intersection.Intersection_results.POINT:
                    pos_zero = lit.Points[0];
                    break;

                case Intersection.Intersection_results.NO_INTERSECTION:
                    //DBG(DBG_KT) << "No intersection at time 0 between supporting lines of wavefrontedges.  Parallel wavefronts crashing (or wavefronts of different speeds becoming collinear).";
                    pos_zero = vertex.Point;
                    break;
                default:
                    //CANNOTHAPPEN_MSG << "Fell through switch which should cover all cases.";
                    assert(false);

                    throw new Exception("Fell through switch which should cover all cases.");

            }

            WavefrontVertex v = new WavefrontVertex(this.Count, pos_zero, vertex, a, b);
            Add(v);
            assert((lit.Result == Intersection.Intersection_results.NO_INTERSECTION) == (v.infinite_speed != InfiniteSpeedType.NONE));

            assert(v.p_at(vertex.Time).AreNear(vertex.Point));
            //   DBG_FUNC_END(DBG_KT);
            return v;





        }
    }
}