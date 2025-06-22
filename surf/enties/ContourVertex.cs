namespace SurfNet
{
    public class ContourVertex
    {

        public int ID { get; private set; }
        public Point2 Point = Point2.NaN;
        public ContourVertex(int id)
        {
            ID = id;
        }

        public ContourVertex(int id, Point2 point)
        {
            ID = id;
            Point = point;
        }

        public static ContourVertex NULL = new ContourVertex(-1);
        public ContourVertex NextInLAV { get; internal set; }
        public ContourVertex PrevInLAV { get; internal set; }

        public override string ToString() => $"C{ID} {Point}";

    }

}
