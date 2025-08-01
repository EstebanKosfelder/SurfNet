namespace SurfNetPolyBool
{
	using SurfNet;
	public class Segment
    {
        public Segment()
        {
        }
        public bool IsInverted;
        public BasicVertex End
        {
            get => Halfedge.Opposite.Vertex;
            set
            {
                    Halfedge.Opposite.Vertex = value;
            }
        }
        public BasicVertex Start
        {
            get => Halfedge.Vertex ;
            set
            {
                    Halfedge.Vertex = value;
              
            }
        }
        public override string ToString()
        {
            return $"{Halfedge} ";
		}
		public Fill MyFill { get; set; }
        public Fill OtherFill { get; set; }
        public BasicHalfedge Halfedge { get; internal set; }
        public bool IsOpposite { get; internal set; }
    }
}