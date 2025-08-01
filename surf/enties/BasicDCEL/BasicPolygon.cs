
namespace SurfNet
{

    public class BasicPolygon
    {
     
        public virtual bool IsHole { get; set; }
        public int Id { get; set; }
        public BasicPolygon(int id, bool isHole) 
        {
            Id = id;
            IsHole = isHole;
        }
        public BasicHalfedge Halfedge { get; set; }
        public BasicPolygonFace PolygonFace { get; set; }
    }


    public class BasicPolygonFace : BasicFace
    {
        public BasicPolygon Polygon { get; set; }

        public BasicPolygonFace(int id, BasicPolygon polygon): base(id)
        {
          Polygon = polygon;
        }
    }




    public class BasicContour : BasicFace
    {
        public BasicContour(int id) : base(id)
        {
        }
    }
}