namespace SurfNet
{
    using static DebugLog;
    using static Mathex;
    public partial class WavefrontPropagator
    {
        private SkeletonStructure sk;
        private EventQueue eq;
        private uint event_ctr_ = 0;
        private bool finalized = false;

        private double time = CORE_ZERO;
        private double last_event_time = CORE_ZERO;
        private double increment = 0.0005;
        private int current_component = -1;
        private int last_event_component = -1;


        //public WavefrontPropagator enter_contour(IEnumerable<Point_2> aPoly, bool aCheckValidity = true, double eps = Tools.EPS)
        //{
        //    if (aCheckValidity)
        //    {
        //        List<Point_2> nPoly = new List<Point_2>();

        //        int idx = 0;
        //        foreach (var p in aPoly)
        //        {
        //            if (nPoly.Count == 0 || !Tools.AreNear(nPoly.Last(), p, eps))
        //            {
        //                nPoly.Add(p);
        //            }
        //            else
        //            {
        //                throw new ArgumentOutOfRangeException(nameof(aPoly), $"point {nPoly[idx]} at index:{idx} is near at previuos point {nPoly[idx - 1]}  ");
        //            }
        //            idx++;
        //        }
        //        if (nPoly.Count > 0 && Tools.AreNear(nPoly.First(), nPoly.Last(), eps))
        //        {
        //            nPoly.RemoveAt(nPoly.Count - 1);
        //        }
        //        if (nPoly.Count < 3)
        //            throw new ArgumentOutOfRangeException(nameof(aPoly), "Degenerate contour (less than 3 non-degenerate vertices).");
        //        aPoly = nPoly;
        //    };

        //    sk.EnterValidContour(aPoly);
        //    return this;
        //}


        private bool no_more_events()
        {

            assert(eq != null);
            return eq.empty() || peak().get_priority().type() == CollapseType.NEVER;
        }

        public WavefrontPropagator(SkeletonStructure p_sk)
        {
            //  eq = new EventQueue();
            sk = p_sk;
        }

        public partial void setup_queue(KineticTriangulation kt);

        public bool propagation_complete() => finalized;


        public double get_time() => time;



        /// <summary>
        ///  Get the current component we are working on.
        /// This is for display purposes only.
        ///    -1 if we should show all. if >=0 we are currently working on this component. 
        /// </summary>
        public int get_current_component() => current_component;

        //public void set_time(double t) { time = t; };
        public void set_increment(double i) => increment = i;

        /// <summary>
        /// Move backwards in time
        /// </summary>

        public void reverse_time() { time -= increment; }

        /// <summary>
        /// Move forward in time, but ignore any event that may have happened
        /// </summary>
        public void advance_time_ignore_event() { time += increment; }

        public void advance_time_ignore_event(double t) { time = t; }

        /// <summary>
        /// Move forward in time by the increment, or until the next event and handle it
        /// </summary>
        public partial void advance_time();

        /// <summary>
        /// Move forward in time to the next event but do not handle it yet.
        /// </summary>
        public partial void advance_time_next();

        /// <summary>
        /// Move forward in time to the next event and handle it.
        /// </summary>
        public partial void advance_step();

        /// <summary>
        /// Process events until done.
        /// </summary>
        public partial void advance_to_end();

        public void reset_time_to_last_event()
        {
            time = last_event_time;
            current_component = last_event_component;
        }

        public CollapseEvent peak()
        {
            assert(eq != null);
            assert(!eq.empty());
            return eq.peak().priority;
        }

        public uint event_ctr() => event_ctr_;

        public partial void do_initial_skips(bool skip_all, int skip_to, double skip_until_time);

        /// <summary>
        /// finish up and create the dcel and whatever else is necessary once the propagation is done.
        /// </summary>
        public partial void finalize();
    }
}