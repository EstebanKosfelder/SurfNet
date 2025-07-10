

using TriangleNet.Topology;

namespace SurfNet
{
    public class SkeletonDCELVertex
    {

        public static readonly SkeletonDCELVertex NaN = new SkeletonDCELVertex(-1,Point2.NaN,double.NaN,null);
        public override string ToString()
        {
            return $"Id:{Id}";
        }
        public SkeletonDCELVertex( int id, Point2 point, double time, SkeletonDCELHalfedge halfedge )
        { 
            
            this.Id = id;
            this.Point = point;
            this.Time = time;
            Halfedge = halfedge;
            if (halfedge != null)
            {
                Halfedge.Vertex = this;
            }


        }


        public static void IncId()
        {
            ctr++;
        }
        private static  int ctr;
        public  int Id;

        public bool IsNaN => this == NaN;
       
        internal bool has_null_point()
        {
            throw new NotImplementedException();
        }

        public Point2 Point {  get; set; }
        public double Time {  get; set; }
        

        public SkeletonDCELHalfedge Halfedge
        {
            get;set;
        }


        public Point_3 point()
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