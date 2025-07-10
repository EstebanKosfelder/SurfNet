using System.Runtime.ConstrainedExecution;
using TriangleNet.Topology.DCEL;

namespace SurfNet
{



    using static DebugLog; 
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

     

        SkeletonDCELHalfedge _next;
        SkeletonDCELHalfedge _prev;

        public SkeletonDCELHalfedge Next
        {
            get => _next;

            set
            {
                if (value != null && _next!=null)
                {
                    DebugLog.Warning($"v{Vertex.Id}.Next != null ant:{_next.Vertex.Id}  new: {value.Vertex.Id}    ");
                }
                
                _next = value;


                if (this.Face != null && value != null && value.Face != Face)
                {
                    if (value.Face != null && value.Face != Face)
                    {
                        DebugLog.Warning($"Face != null ant:{value.Face}  new: {Face}    ");
                    }

                    value.Face = Face;
                }
                else if(this.Face ==null && value != null  && value.Face!=null)
                {
                    Face = value.Face;
                }   

                if (value != null && value.Prev!=this)
                {
                    value.Prev = this;
                    
                }
            }
        }

        public SkeletonDCELHalfedge Prev
        {
            get => _prev;

            set
            {
                if (value != null && _prev != null)
                {
                    DebugLog.Warning($"v{Vertex.Id}.Prev != null ant:{_prev.Vertex.Id} new: {value.Vertex.Id}    ");
                }
                _prev = value;
                if (value != null && value.Next != this)
                {
                    value.Next = this;
                }
            }
        }

        public IEnumerable<SkeletonDCELHalfedge> Circulation(int max = 100000)
        {
            SkeletonDCELHalfedge curr = this;
            SkeletonDCELFace face = curr.Face;
            int i = 0;
            if (curr == null)
            {
                throw new InvalidOperationException($"{this} has not Halfedge asigned");
            }
            do
            {
                yield return curr;
                curr = curr.Next;

                if (curr == null)
                {
                    throw new InvalidOperationException($"{this} has not Halfedge asigned");
                }
                if (face != curr.Face)
                {
      //              throw new InvalidOperationException($"{this} has not a only Face {curr}");
                }
                if (max < i++)
                {
                    throw new GFLException($" {this} iterate more than {max} "); break;
                }
            } while (curr != this);
        }
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
        public SkeletonDCELFace Face { get; internal set; }

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
