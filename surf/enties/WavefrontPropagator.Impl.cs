namespace SurfNet
{
    using static DebugLog;
    using static Mathex;
    public partial class WavefrontPropagator
    {
        public partial void setup_queue(KineticTriangulation kt)
        {
            eq = new EventQueue(kt.triangles);
            kt.set_queue(eq);
        }

        public partial void advance_time()
        {
            time += increment;

            if (!propagation_complete())
            {
                if (no_more_events())
                {
                    advance_step();
                }
                else
                {
                    double want_time = time;
                    while (!propagation_complete() && want_time > eq.peak().priority.time())
                    {
                        advance_step();
                    }
                    time = want_time;
                }
            }
        }

        public partial void advance_time_next()
        {
            throw new NotImplementedException();
            //if (!no_more_events())
            //{
            //    EventQueueItem next = eq.peak();
            //    time = next.get_priority().time();
            //  current_component = peak().get_priority().t.component;

            //}
        }

        public partial void advance_step()
        {


            ///DBG_INDENT_LEVEL_STORE;
            ///DBG_FUNC_BEGIN(///DBG_PROP);

            if (!no_more_events())
            {
                CollapseEvent next = eq.peak().priority;
                time = next.time();
                if (sk.get_kt().restrict_component() != 0)
                {
                    current_component = peak().t.component;
                }
                ++event_ctr_;
                //VLOG(2) << " event#" << event_ctr_ << " @ " << CGAL::to_double(time);
                sk.kt.handle_event(next);
                ///DBG(///DBG_PROP) << " event#" << event_ctr_ << " handling done.  Processing pending PQ updates.";
                eq.process_pending_updates(time);
                ///DBG(///DBG_PROP) << " event#" << event_ctr_ << " PQ updates done.  Time is now " << CGAL::to_double(time);
                //LOG(INFO) << " event#" << event_ctr_ << " done.  Time is now " << CGAL::to_double(time);

                last_event_time = time;
                last_event_component = current_component;
            }

            if (no_more_events())
            {
                //LOG(INFO) << "All done.";
                finalize();
            }
            else
            {
                CollapseEvent next_event = eq.peak().priority;
                ///DBG(///DBG_PROP) << " event#" << (event_ctr_+1) << " will be " << next_event;
                assert((next_event.t.component == current_component && next_event.time() >= last_event_time) ||
                        next_event.t.component > current_component);
                if (eq.Count >= 2)
                {
                    ///DBG(///DBG_PROP) << "   event child in heap: " << eq.peak(1).get_priority();
                    if (eq.Count >= 3)
                    {
                        ///DBG(///DBG_PROP) << "   event child in heap: " << eq.peak(2).get_priority();
                    }
                }
            }
            ///DBG_FUNC_END(///DBG_PROP);
            ///DBG_INDENT_LEVEL_CHECK;
        }

        public partial void advance_to_end()
        {
            while (!propagation_complete())
            {
                advance_step();
            }
        }

        public partial void do_initial_skips(bool skip_all, int skip_to, double skip_until_time)
        {
            if (skip_all)
            {
                advance_to_end();
            }
            else
            {
                while (!propagation_complete() &&
                       skip_to > event_ctr() + 1)
                {
                    advance_step();
                };
            }
            if (skip_until_time > CORE_ZERO)
            {
                while (!propagation_complete() &&
                       (no_more_events() || skip_until_time > peak().time()))
                {
                    advance_step();
                }
                if (skip_until_time > get_time())
                {
                    advance_time_ignore_event(skip_until_time);
                }
            }
        }

        public partial void finalize()
        {
            ///DBG_FUNC_BEGIN(///DBG_PROP);
            assert(no_more_events());
            if (!finalized)
            {
                ///DBG(///DBG_PROP) << "Calling create_remaining_skeleton_dcel()";
                sk.kt.create_remaining_skeleton_dcel();
                finalized = true;
                current_component = -1;
                ///DBG(///DBG_PROP) << "Finalized.";

                sk.kt.update_event_timing_stats(-1);
            }
            ///DBG_FUNC_END(///DBG_PROP);
        }
    }
}