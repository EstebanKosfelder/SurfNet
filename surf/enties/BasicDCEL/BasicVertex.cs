namespace SurfNet
{
    public  class BasicVertex
    {


        public BasicVertex(int id, Point2 point,int originalId,int pointId  )
        {
            this.Id = id;
            this.OriginalId = id;
            this.Point = point;
            this.PointId = pointId;

        }
        public BasicVertex(int id, Point2 point, int pointId):this(id,point, id,pointId)
        {
        
        }   

        public int Id { get; private set; } 
        public int OriginalId { get; internal set; }
        public int PointId { get; set; }
        public Point2 Point { get; internal set; }

        private BasicHalfedge _halfedge;
        public BasicHalfedge Halfedge { get => _halfedge; set {
                _halfedge = value; } 
        }
        public int degree { get; internal set; } = 0;
        public int reflex_beveling_add { get; internal set; } = 0;

        public override string ToString()
        {
            return $"{Id:D3} {OriginalId:D3} {PointId:D3} {Point}";
        }
        
    }
}