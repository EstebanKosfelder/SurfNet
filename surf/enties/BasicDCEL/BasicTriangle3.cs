
namespace SurfNet
{
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
}