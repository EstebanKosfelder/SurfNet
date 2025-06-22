

namespace SurfNet
{
    using System.Collections.Generic;
    public class SkeletonDCEL
    {
        public IEnumerable<SkeletonDCELHalfedge> halfedges;

        internal SkeletonDCELFace setup_new_input_edge(WavefrontEdge? buddy_wavefront)
        {
            throw new NotImplementedException();
            //SkeletonDCELFace face = new_face();
            //SkeletonDCELCbb ccb = new_outer_ccb();
            //SkeletonDCELHalfedge halfedge = null;

            //ccb.set_face(face);

            //if (buddy_wavefront!=null)
            //{ /* buddy wavefront already set up */
            //    halfedge = (buddy_wavefront.skeleton_face.outer_ccbs_begin()).opposite();
            //    assert(halfedge.is_on_outer_ccb());
            //    assert(halfedge.outer_ccb() == null);
            //}
            //else
            //{
            //    halfedge = new_edge();
            //}
            //halfedge.is_emanating_input_ = true;
            //halfedge.set_outer_ccb(ccb);

            //face.add_outer_ccb(ccb, halfedge);
            //return face;
        }

        internal SkeletonDCELCbb new_outer_ccb()
        {
            throw new NotImplementedException();
        }

        internal SkeletonDCELFace new_face()
        {
            throw new NotImplementedException();
        }

        internal SkeletonDCELHalfedge new_edge()
        {
            throw new NotImplementedException();
        }

        internal void set_num_v_skew(int v)
        {
            throw new NotImplementedException();
        }

        internal void set_number_of_points_and_curves()
        {
            throw new NotImplementedException();
        }

        internal X_monotone_curve new_segment(Segment_3 segment_3)
        {
            throw new NotImplementedException();
        }

        internal SkeletonDCELHalfedge halfedges_begin()
        {
            throw new NotImplementedException();
        }

        internal SkeletonDCELHalfedge halfedges_end()
        {
            throw new NotImplementedException();
        }

        internal void assert_sane()
        {
            throw new NotImplementedException();
        }

        internal X_monotone_curve new_ray(Ray3 ray_3)
        {
            throw new NotImplementedException();
        }

        internal SkeletonDCELVertex new_vertex()
        {
            throw new NotImplementedException();
        }

        internal Point_3 new_point(Point_3 point_3)
        {
            throw new NotImplementedException();
        }
    }


    public class SkeletonDCELFace 
    {
        private bool is_beveling_face_;
        public void set_is_beveling_face(bool is_beveling_face) { is_beveling_face_ = is_beveling_face; }
        public bool is_beveling_face() { return is_beveling_face_; }



        public SkeletonDCELFace() { id = (ctr++); }
        private static int ctr;
        public int id;

        public SkeletonDCELFace(DcelVertex generator) 
        {
        }

        public SkeletonDCELFace(DcelVertex generator, DcelHalfEdge edge) 
        {
        }

        internal void add_outer_ccb(SurfNet.SkeletonDCELCbb ccb, SurfNet.SkeletonDCELHalfedge halfedge)
        {
            throw new NotImplementedException();
        }

        public  SkeletonDCELHalfedge outer_ccbs_begin()
        {
            throw new NotImplementedException();
        }

        internal SkeletonDCELHalfedge outer_ccbs_end()
        {
            throw new NotImplementedException();
        }
    }


    public class SkeletonDCELVertex
    {
        public SkeletonDCELVertex() { id = ctr++; }
        private int ctr;
        public readonly int id;
        internal bool has_null_point()
        {
            throw new NotImplementedException();
        }

        internal Point_3 Point()
        {
            throw new NotImplementedException();
        }

        internal Point_3 point()
        {
            throw new NotImplementedException();
        }

        internal void set_halfedge(SkeletonDCELHalfedge start)
        {
            throw new NotImplementedException();
        }

        internal void set_point(Point_3 pp)
        {
            throw new NotImplementedException();
        }

        
    }
}