using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Numerics;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;

namespace SurfNet
{
 

//  public struct Point_3
//    {

//    }
//  public  class SkeletonDCELVertexBase : List<Point_3> {
 
//};
///** A halfedge.  The Segment_3/Ray_3 is undirected. */
//public class SkeletonDCELHalfedgeBase
//    {

//   // public SurfNet.Arr_halfedge_base<boost.variant<Segment_3, Ray_3>>
//            internal bool is_emanating_input_ = false;

//        ///<summary>
//        /// Is this an input edge that is emmanating a wavefront.
//        /// True for input edges on the side where we do a wavefront propagation.
//        ///</summary>
//          public bool is_emanating_input() { return is_emanating_input_; }

//        ///<summary>
//        /// Is this an input edge for which we do a wavefront propagation on at least one side.
//        /// May only be called when class is specialized by Arr_halfedge */
//        ///</summary>
//        internal partial bool is_input() ;

//};
//public class SkeletonDCELFaceBase : List<FaceBase>
//{
//  private  bool is_beveling_face_;
//public void set_is_beveling_face(bool is_beveling_face) { is_beveling_face_ = is_beveling_face; };
//        public  bool is_beveling_face() { return is_beveling_face_; };

//};

//public class SkeletonDCEL : public SurfNet.Arr_dcel_base<SkeletonDCELVertexBase, SkeletonDCELHalfedgeBase, SkeletonDCELFaceBase>
//{
    
////private:
////    using Base = SurfNet.Arr_dcel_base<SkeletonDCELVertexBase, SkeletonDCELHalfedgeBase, SkeletonDCELFaceBase>;
////    using Curve = Halfedge.X_monotone_curve;

//private List<Point_3> points;
//    private List<Curve> curves;
//    private int num_v_skew = -1;
//private int get_expected_num_v() {
//      var num_he = size_of_halfedges();
//assert(num_he % 2 == 0);
//var num_f = size_of_faces();
//return 2 + num_he / 2 - num_f;
//    }

//  private Face setup_new_input_edge(WavefrontEdge buddy_wavefront);

///* skew normally is -1, since we do not need a point for the infinite vertex.
// * However, if we don't have an infinite vertex (because we only do an
// * interior SK), then the skew does not hold and needs to be one.
// */
//void set_num_v_skew(int skew)
//{
//    assert(points.size() == 0);
//    num_v_skew = skew;
//}
//void set_number_of_points_and_curves()
//{
//    var num_he = size_of_halfedges();
//    var num_f = size_of_faces();
//    assert_sane();
//    points.reserve(get_expected_num_v() + num_v_skew);
//    curves.reserve(num_he / 2);
//}

//Point_3 new_point(Point_3 p)
//{
//    points.emplace_back(std.forward<Point_3>(p));
//    Point_3* pp = &points.back();
//    return pp;
//}
//Curve new_segment(Segment_3 s)
//{
//    curves.emplace_back(std.forward<Segment_3>(s));
//    Curve p =&curves.back();
//    return p;
//}
//Curve new_ray(Ray_3 s)
//{
//    curves.emplace_back(std.forward<Ray_3>(s));
//    Curve p = curves.back();
//    return p;
//}

////private    enum class SegmentType { INPUT, ARC, OFFSET };
////static void write_ipe_segment(Segment_2 s, SegmentType t);
////public void write_ipe(std.ostream& os, std.string ffset_spec) ;
////void write_obj(std.ostream& os) const;

////public:
// /  protected:


//    /// <summary>

//    /// </summary>
//    /// <param name="offsetting_distance"></param>
//    /// <param name="he"></param>
//    /// <returns>sign</returns>
//public  int  arc_offset_do_intersect(double offsetting_distance, Halfedge he) ;
//    public Point_2 arc_offset_get_intersect( Plane_3 offset_plane, Halfedge he) ;
//public (bool, Segment_2) make_one_offset_segment(double offsetting_distance, Halfedge start, HashSet<Halfedge> visited) ;

//public  List<Segment_2>  make_offset(const NT& offsetting_distance) const;
//public static List<double> parse_offset_spec( offset_spec);

//    public void assert_sane() { }
//};
//using SkeletonDCELVertex = SkeletonDCEL.Vertex;
//using SkeletonDCELHalfedge = SkeletonDCEL.Halfedge;
//using SkeletonDCELCbb = SkeletonDCEL.Outer_ccb;
//using SkeletonDCELFace = SkeletonDCEL.Face;

//inline bool
//SkeletonDCELHalfedgeBase.
//is_input() const {
//  static_assert(std.is_base_of<SkeletonDCELHalfedgeBase, SkeletonDCELHalfedge>.value);
//const SkeletonDCELHalfedge* const he = static_cast<const SkeletonDCELHalfedge*>(this);
//assert(he);
//assert(he->opposite());
//return he->is_emanating_input() || he->opposite()->is_emanating_input();
//}
}
