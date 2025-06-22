using System.Runtime.ConstrainedExecution;

namespace SurfNet
{




public class SkeletonDCELHalfedge
    {

        internal bool is_emanating_input_ = false;
        ///<summary>
        /// Is this an input edge that is emmanating a wavefront.
        /// True for input edges on the side where we do a wavefront propagation.

        ///</summary>
        public  bool is_emanating_input()=>is_emanating_input_; 
        /** Is this an input edge for which we do a wavefront propagation on at least one side.
         *
         * May only be called when class is specialized by Arr_halfedge */
        bool is_input() { return false; }


        public SkeletonDCELHalfedge() { id = (ctr++); }
        private static int ctr;
        public int id;



        public bool is_on_outer_ccb()
        {
            throw new NotImplementedException();
        }

        public SkeletonDCELHalfedge next()
        {
            throw new NotImplementedException();
        }

        public  SkeletonDCELCbb outer_ccb()
        {
            throw new NotImplementedException();
        }

        public  void set_outer_ccb(SkeletonDCELCbb ccb)
        {
            throw new NotImplementedException();
        }

        public SkeletonDCELVertex vertex()
        {
            throw new NotImplementedException();
        }

        public void set_curve(X_monotone_curve p)
        {
            throw new NotImplementedException();
        }

       

        internal bool has_null_curve()
        {
            throw new NotImplementedException();
        }

        internal SkeletonDCELHalfedge opposite()
        {
            throw new NotImplementedException();
        }

        internal void set_vertex(SkeletonDCELVertex new_v)
        {
            throw new NotImplementedException();
        }

        internal SkeletonDCELHalfedge prev()
        {
            throw new NotImplementedException();
        }

        internal void set_next(SkeletonDCELHalfedge next)
        {
            throw new NotImplementedException();
        }

        internal void set_prev(SkeletonDCELHalfedge prev_he)
        {
            throw new NotImplementedException();
        }
    }
}
