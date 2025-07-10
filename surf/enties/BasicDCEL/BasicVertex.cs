namespace SurfNet
{
    public  class BasicVertex
    {
        public BasicVertex(int id, Point2 point, int degree)
        {
            this.Id = id;
            this.p = point;
            this.degree = 0;
        }
        public int Id { get; private set; } 
        public Point2 p { get; internal set; }
        public int degree { get; internal set; }
        public int reflex_beveling_add { get; internal set; } = 0;
        
    }
}