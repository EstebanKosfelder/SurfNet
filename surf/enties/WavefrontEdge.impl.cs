namespace SurfNet
{
    using static DebugLog;
    using static Mathex;
    public partial class WavefrontEdge
    {


        //std::ostream&
        //operator<<(std::ostream& os, const WavefrontEdge& e) {
        //  os << "wfe";
        //  DEBUG_STMT(os << "#" << e.id );
        //        os << "("
        //    << e.vertices[0] << ","
        //    << e.vertices[1] << ")";
        //  return os;
        //}

        /** returns when this edge will collapse.
         *
         * If the two vertices are parallel or moving away from one another,
         * they will NEVER collapse.  Even if the two end-points are
         * conincident right now but are moving away from another, we will consider
         * this as not collapsing.
         *
         * Otherwise, they will collapse at some point in time.  Since
         * we are asking, we assume (and we assert() during debugging),
         * that this will be in the future (or now).
         */

        private partial EdgeCollapseSpec compute_collapse(double time_now)
        {
            //DBG_FUNC_BEGIN(DBG_KT);
            //DBG(DBG_KT) << "Computing edge collapse time for " << * this;

            EdgeCollapseSpec res;
            WavefrontVertex wfv0 = vertices[0];
            WavefrontVertex wfv1 = vertices[1];
            assert(wfv0 != null);
            assert(wfv1 != null);
            Vector2 v0 = wfv0.velocity;
            Vector2 v1 = wfv1.velocity;

            //DBG(DBG_KT) << "v0" << CGAL_vector(v0);
            // DBG(DBG_KT) << "v1" << CGAL_vector(v1);
            var o = orientation(v0, v1);

            if (o != OrientationEnum.LEFT_TURN)
            {
                /* If the two wavefront vertices move away from each other
                 * or in parallel, this edge will never collapse.
                 */
                if (o == OrientationEnum.RIGHT_TURN)
                {
                    //  DBG(DBG_KT) << "Orientation is right turn";
                    // let's consider two points that are identical right now but moving away from another as not collapsing.
                    res = new EdgeCollapseSpec(EdgeCollapseType.PAST);
                }
                else
                {
                    // DBG(DBG_KT) << "Orientation is collinear";
                    assert(o == OrientationEnum.COLLINEAR);

                    /* Previously we computed the distance at time_now.  I wonder why.
                     * If we then claim that the edge is always collapsing, then it should
                     * suffice to compute the distance at t=0. */
                    double sqdist = squared_distance(wfv0.pos_zero, wfv1.pos_zero);
                    //DBG(DBG_KT) << "sqdist zero: " << CGAL::to_double(sqdist);

                    {
                        Point2 p0 = (wfv0.p_at(time_now));
                        Point2 p1 = (wfv1.p_at(time_now));
                        double sqdistnow = squared_distance(p0, p1);
                        //     DBG(DBG_KT) << "sqdist now : " << sqdistnow;
                        if (sqdist.AreNear(CORE_ZERO))
                        {
                            assert(sqdist.AreNear(sqdistnow), $"{sqdist} == {sqdistnow}");
                        }
                        else
                        {
                            assert(sqdistnow > CORE_ZERO);
                        }
                    }

                    if (sqdist == CORE_ZERO)
                    {
                        // DBG(DBG_KT) << "Distance is zero now.";
                        res = new EdgeCollapseSpec(EdgeCollapseType.ALWAYS, time_now);
                    }
                    else
                    {
                        //  DBG(DBG_KT) << "Distance is not zero.";
                        res = new EdgeCollapseSpec(EdgeCollapseType.NEVER);
                    }
                }
            }
            else
            {
                // DBG(DBG_KT) << "Orientation is left turn";
                /* Note that, by construction, if you start on the wavefront edge, go out v0,
                 * and go back v1, you end up on the wavefrong edge again.  Or, in other words,
                 * v0-v1 is collinear with the direction of the wavefront edge e.
                 *
                 * Now, we want to know when e=AB collapses.  So we can restrict ourselves to
                 * consider the projection of A+t*v0 and B+t*v1 to e itself.  Once they meet,
                 * the edge collapses.  Equivalently, once the length of projected t*v0 and
                 * projected t*v1 equals the length of e, the edge collapses.
                 *
                 * Let d be the vector AB (i.e., B-A).  The dot product v0.d is the length of
                 * the projected v0 times the length of d.  So v0.d/|d| is the length of the
                 * projected v0.  Likewise, v1.d/|d| is the length of the projected v1.
                 * Thus e will collapse at time t := |d| / ( v0.d/|d| - v1.d/|d| ) ==
                 * == |d|^2 / ( v0.d - v1.d ) == |d|^2 / ( (v0 - v1).d) ==
                 * == d.d /  ( (v0 - v1).d).
                 *
                 * Isn't that interesting?  Note how we compare the length of d projected
                 * onto d with the length of (v0 - v1) projected onto d?  Remember that, as
                 * previously mentioned, by construction, (v0 - v1) is in the same direction
                 * as d.  Thus, this quotient will have the same value regardless of
                 * which line we project both d and (v0 - v1) -- as long as the line is not
                 * orthogonal to (v0 - v1).  So we might just as well project onto (1,0) and
                 * only consider their x-coordinates (if d is not exactly vertical).
                 * So t == d_x / (v0_x - v1_x)
                 *
                 * This also makes sense when looking at it another way.  Consider the points
                 * of the projection of A+t*v0 and B+t*v1 onto e.  Since they are on e, they
                 * will become incident when and only when their x-coordinates is the same.
                 */
                double edge_delta;
                double wfvs_delta;
                if (!supporting_line.l.is_vertical())
                {
                    edge_delta = wfv1.pos_zero.X - wfv0.pos_zero.X;
                    wfvs_delta = v0.X - v1.X;
                }
                else
                { /* l is vertical, do the same with y-coordinates */
                    edge_delta = wfv1.pos_zero.Y - wfv0.pos_zero.Y;
                    wfvs_delta = v0.Y - v1.Y;
                }

                assert(edge_delta != 0);
                assert(wfvs_delta != 0);
                double time = edge_delta / wfvs_delta;

                Point2 p0 = (wfv0.p_at(time));
                Point2 p1 = (wfv1.p_at(time));
                double sqdistnow = squared_distance(p0, p1);

                if (!sqdistnow.IsZero())
                {
                    sqdistnow.IsZero();
                }

                //DBG(DBG_KT) << "future edge collapse: " << CGAL::to_double(time);
                //DBG(DBG_KT) << "time_now            : " << CGAL::to_double(time_now);
                assert(time > time_now);
                res = new EdgeCollapseSpec(EdgeCollapseType.FUTURE, time);
            }
            //DBG(DBG_KT) << "returning " << res;
            //DBG_FUNC_END(DBG_KT);
            return res;
        }

#if !SURF_NDEBUG

        private partial void assert_edge_sane(int collapsing_edge)
        {

            assert(0 <= collapsing_edge && collapsing_edge < 3);
            assert(incident_triangle_ != null);
            assert(incident_triangle_.wavefront(collapsing_edge) == this);
            assert(vertices[0] != null);
            assert(vertices[1] != null);
        }

#endif


        /** Duplicate this edge in the course of a split event.
         *
         * This is just a simple helper that creates two copies of this edge and marks
         * the original as dead.
         *
         * It is the job of the caller (KineticTriangulation) to then give us new vertices.
         */

        public partial EdgePtrPair split(WavefrontEdgeList wavefront_edges)
        {
            assert(vertices[0] != null);
            assert(vertices[1] != null);
            assert(vertices[0].incident_wavefront_edge(1) == this);
            assert(vertices[1].incident_wavefront_edge(0) == this);
            set_dead();

            wavefront_edges.Add(new WavefrontEdge(vertices[0], null, supporting_line, incident_triangle_, skeleton_face, is_beveling));
            var pea = wavefront_edges.back();
            wavefront_edges.Add(new WavefrontEdge(null, vertices[1], supporting_line, incident_triangle_, skeleton_face, is_beveling));
            var peb = wavefront_edges.back();

            return new EdgePtrPair(pea, peb);
        }
    }
}