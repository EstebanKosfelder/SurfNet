
namespace SurfNet
{
    using static DebugLog;

    public class BasicTriangle3 : Array3Nulleable<BasicTriangle>
    {
        public override string ToString()
        {
            return $"n:[{string.Join(", ", array.Select(e => (e != null ? e.Id.ToString() : " *").PadLeft(3)))}]";
        }
        public BasicTriangle3()
        {
        }

        public BasicTriangle3(params BasicTriangle?[] values) : base(values)
        {
        }
    }
    public class BasicVertex3 : Array3<BasicVertex>
    {
        public override string ToString()
        {
            return $"v:[{string.Join(", ", array.Select(e => (e != null ? e.Id.ToString() : " *").PadLeft(3)))}]";
        }
        public BasicVertex3()
        {
        }

        public BasicVertex3(params BasicVertex[] values) : base(values)
        {
        }
    }
    public  class BasicFace
    {
       // public BasicFaceInfo info() => new BasicFaceInfo();

        internal bool is_constrained(int i)
        {
            throw new NotImplementedException();
        }
    }

    public  class BasicTriangle
    {
        public BasicTriangle(int id)
        {
            Id = id;
        }
        public int Id = -1;
        public BasicVertex3 Vertices = new BasicVertex3();
        public BasicTriangle3 Neighbors = new BasicTriangle3();

        public  int index(BasicTriangle needle)
        {
            
                var  idx =
                  (Neighbors[0] == needle) ? 0 :
                  (Neighbors[1] == needle) ? 1 :
                  (Neighbors[2] == needle) ? 2 : throw new Exception( $"no existe Neighbors {needle} en {this} ") ;

                return idx;
            
        }
    }
}