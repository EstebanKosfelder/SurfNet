namespace SurfNetPolyBool
{
	using SurfNet;
	public class Node : IEquatable<Node>
    {
        public override string ToString()
        {
            if (Seg==null) return $"* {this.Pt} - {this.Other.Pt}";
			return $"{ Seg.Halfedge} ";
        }
        public Node Status { get; set; }

    public Node Other { get; set; }

    public Node Ev { get; set; }

    public Node Prev { get; set; }

    public Node Next { get; set; }

    public bool IsRoot { get; set; }

    public Action Remove { get; set; }

    public bool IsStart { get; set; }
    Point2 _pt;
	public Point2 Pt 
        {  get=>_pt; 
            set { 
                _pt = value;
            }       
        }

    public Segment Seg { get; set; }

    public bool Primary { get; set; }

    public bool Equals(Node other)
    {
        if (ReferenceEquals(null, other)) return false;
        if (ReferenceEquals(this, other)) return true;
        return Equals(Prev, other.Prev) && Equals(Next, other.Next) && IsRoot == other.IsRoot;
    }

    public override bool Equals(object obj)
    {
        if (ReferenceEquals(null, obj)) return false;
        if (ReferenceEquals(this, obj)) return true;
        if (obj.GetType() != GetType()) return false;
        return Equals((Node)obj);
    }
    public void LinkPrev(Node other)
        {
            this.Prev = other;
            if (other != null)
            {
                other.Next = this;
			}
		}
        public void LinkNext(Node other)
        {
            this.Next = other;
            if (other != null)
            {
                other.Prev = this;
            }
		}

		public override int GetHashCode()
    {
        unchecked
        {
            var hashCode = (Prev != null ? Prev.GetHashCode() : 0);
            hashCode = (hashCode * 397) ^ (Next != null ? Next.GetHashCode() : 0);
            hashCode = (hashCode * 397) ^ IsRoot.GetHashCode();
            return hashCode;
        }
    }
}
}