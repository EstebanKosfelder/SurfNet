

namespace SurfNet
{
    using System.Collections.Generic;
    using static DebugLog;
    public class SkeletonDCEL
    {
        public List<SkeletonDCELHalfedge> halfedges = new List<SkeletonDCELHalfedge>();

        public List<SkeletonDCELVertex> vertices = new List<SkeletonDCELVertex>();

        public List<SkeletonDCELFace> faces = new List<SkeletonDCELFace>();                // The faces container.
        public List<SkeletonDCELCbb> out_ccbs = new List<SkeletonDCELCbb>();             // The outer CCBs.
        private int num_v_skew;

        //internal SkeletonDCELFace setup_new_input_edge(WavefrontEdge? buddy_wavefront)
        //{
            
        //    //SkeletonDCELFace face =new_face();
        //    //SkeletonDCELCbb ccb = new_outer_ccb();
        //    //SkeletonDCELHalfedge halfedge = null;

        //    //ccb.set_face(face);

        //    //if (buddy_wavefront != null)
        //    //{ /* buddy wavefront already set up */
        //    //    halfedge = (buddy_wavefront.skeleton_face.outer_ccbs_begin()).opposite();
        //    //    assert(halfedge.is_on_outer_ccb());
        //    //    assert(halfedge.outer_ccb() == null);
        //    //}
        //    //else
        //    //{
        //    //    halfedge = new_edge();
        //    //}
        //    //halfedge.is_emanating_input_ = true;
        //    //halfedge.set_outer_ccb(ccb);

        //    //face.add_outer_ccb(ccb, halfedge);
        //    //return face;
        //}

        internal SkeletonDCELCbb new_outer_ccb()
        {
            out_ccbs.Add(new SkeletonDCELCbb());
            return out_ccbs.Last();
        }

        internal SkeletonDCELFace new_face(SkeletonDCELHalfedge he )

        {
            
            faces.Add(new SkeletonDCELFace());

            var f =faces.Last();
            f.Halfedge=he;
            he.Face = f;    
            return f;
        }

        internal SkeletonDCELHalfedge new_edge()
        {
            var ha = new SkeletonDCELHalfedge();
            var hb = new SkeletonDCELHalfedge();
            ha.Opposite = hb;
            hb.Opposite = ha;
            halfedges.Add(ha);
            halfedges.Add(hb);

            return ha;
        }

        internal void set_num_v_skew(int skew)
        {
           // assert(points.size() == 0);
            num_v_skew = skew;
        }

        internal void set_number_of_points_and_curves()
        {
           // throw new NotImplementedException();
        }

        internal X_monotone_curve new_segment(Segment_3 segment_3)
        {
           return new X_monotone_curve();
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
           // throw new NotImplementedException();
        }

        internal X_monotone_curve new_ray(Ray3 ray_3)
        {
            return new X_monotone_curve();
        }


        internal SkeletonDCELVertex new_vertex(Point2 p , double time)
        {
            vertices.Add(new SkeletonDCELVertex(vertices.Count,p, time,new_edge()));
            SkeletonDCELVertex.IncId();
            return vertices.Last();
        }

        internal SkeletonDCELVertex new_vertex(int id, Point2 p, double time)
        {
         var  v = new SkeletonDCELVertex(id, p, time, new_edge());
            return v;
            
        }


        internal Point_3 new_point(Point_3 point_3)
        {
            return point_3;
        }
    }
}