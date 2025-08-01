
namespace SurfNet
{
    public partial  class BasicTriangle:BasicFace
    {
        public BasicTriangle(int id):base(id)
        {
            this.Vertices = new BasicVertex3();
            this.Neighbors = new BasicTriangle3();
        }
        public BasicVertex3 Vertices = new BasicVertex3();
        public BasicTriangle3 Neighbors = new BasicTriangle3();
        public override string ToString() => $"{Id} v:{Vertices} n:{Neighbors}";

        public  int index(BasicTriangle needle)
        {
            
                var  idx =
                  (Neighbors[0] == needle) ? 0 :
                  (Neighbors[1] == needle) ? 1 :
                  (Neighbors[2] == needle) ? 2 : throw new Exception( $"no existe Neighbors {needle} en {this} ") ;

                return idx;
            
        }

        public Vector2[] Normals
        {
            get
            {

                return new[] {
                     new Vector2(Vertices[1].Point, Vertices[0].Point).Normal()
                    ,new Vector2(Vertices[2].Point, Vertices[1].Point).Normal()
                    ,new Vector2(Vertices[0].Point, Vertices[2].Point).Normal()

                };

            }
        }
    }
}