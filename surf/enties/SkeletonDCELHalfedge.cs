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

        public SkeletonDCELCbb OutCbb { get; set; }

        public SkeletonDCELVertex Vertex { get; set; }

        public SkeletonDCELHalfedge Prev { get; set; }
        public SkeletonDCELHalfedge Next { get; set; }

        public bool is_on_outer_ccb()
        {
            return OutCbb != null;
        }

     
        public  SkeletonDCELCbb outer_ccb()
        {
           return OutCbb;
        }

        public  void set_outer_ccb(SkeletonDCELCbb ccb)
        {
           OutCbb = ccb;
        }

        public SkeletonDCELVertex vertex()
        {
            return Vertex;
        }

        public void set_curve(X_monotone_curve p)
        {
            throw new NotImplementedException();
        }

       

        internal bool has_null_curve()
        {
            throw new NotImplementedException();
        }

        public SkeletonDCELHalfedge Opposite { get; set; }

        internal SkeletonDCELHalfedge opposite()
        {
            return Opposite;
        }

        internal void set_vertex(SkeletonDCELVertex new_v)
        {
            Vertex = new_v;
        }

        internal SkeletonDCELHalfedge prev()
        {
            return this.Prev;
        }

        internal SkeletonDCELHalfedge next()
        {
            return this.Next;
        }
        public void set_prev(SkeletonDCELHalfedge he)
        {
            
            this.Prev = he;
          
            he.Next = this;
        }
        internal void set_next(SkeletonDCELHalfedge he)
        {
            this.Next = he;

            he.Prev = this;
        }

       
    }
}
