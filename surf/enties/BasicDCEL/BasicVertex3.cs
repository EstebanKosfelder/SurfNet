
namespace SurfNet
{
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
}