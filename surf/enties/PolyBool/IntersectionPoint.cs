namespace SurfNetPolyBool
{
	using SurfNet;
	public class IntersectionPoint
    {
        public Point2 Pt { get; set; }
        public int AlongA { get; set; }
        public int AlongB { get; set; }

        public IntersectionPoint(int alongA, int alongB, Point2 pt)
        {
            AlongA = alongA;
            AlongB = alongB;
            Pt = pt;
        }
    }
}