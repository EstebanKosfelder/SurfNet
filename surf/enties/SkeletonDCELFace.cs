

namespace SurfNet
{
    public class SkeletonDCELFace 
    {
        private bool is_beveling_face_;
        public void set_is_beveling_face(bool is_beveling_face) { is_beveling_face_ = is_beveling_face; }
        public bool is_beveling_face() { return is_beveling_face_; }

        public string ToString() { return $"Id:{Id}"; }

        public SkeletonDCELFace() { Id = (ctr++); }
        private static int ctr;
        public int Id { get; private set; }

        public SkeletonDCELHalfedge Halfedge { get; set; }

        public SkeletonDCELFace(DcelVertex generator) 
        {
        }

        public SkeletonDCELFace(DcelVertex generator, DcelHalfEdge edge) 
        {
        }

        internal void add_outer_ccb(SurfNet.SkeletonDCELCbb ccb, SurfNet.SkeletonDCELHalfedge halfedge)
        {
        // TODO       throw new NotImplementedException();
        }

        public  SkeletonDCELHalfedge outer_ccbs_begin()
        {
            return null;// throw new NotImplementedException();
        }

        internal SkeletonDCELHalfedge outer_ccbs_end()
        {
            return null;// throw new NotImplementedException();
        }
    }
}