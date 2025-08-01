namespace SurfNet
{
    [Flags]
    public enum BasicHalfedgeType : byte
    {
        None = 0,
        Countour = 1 << 0,
        Constraint = 1 << 1,
        IsHole = 1 << 2,
        IsOutside = 1 << 3,
        Interior = 1 << 4,
        Deleted = 1 << 7,
    }

    public class BasicHalfedge
    {
        public int Id;
        public BasicHalfedgeType Type { get; set; } = BasicHalfedgeType.None;

        public BasicHalfedge(int id, BasicHalfedgeType type)
        {
            this.Id = id;
            this.Type = type;
        }
        public int TagPoint { get; set; }

        BasicVertex _vertex;
        public BasicVertex Vertex {
            get => _vertex;
            set {if (_vertex != null)
                { 
                }
              _vertex = value;
                 } 
        }

        private BasicFace face;

        public BasicFace Face
        {
            get => face;
            set
            {
                if (face != null)
                {
                }
                face = value;
            }
        }

        public BasicHalfedge Next { get; set; }
        public BasicHalfedge Prev { get; set; }
        public BasicHalfedge Opposite { get; set; }
        string PointId()
        {
            if (Vertex == null) return $"({TagPoint})";
            return $"{Vertex.PointId:D3}";
        }
        string VertexId()
        {
            if (Vertex == null) return $"*";
            return $"{{{Vertex.Id:D3}}}";
        }
        public override string ToString()
        {
            //  return $"[{Id:D3}] {Prev?.Vertex?.Id:D3}<- v:{Vertex?.Id:D3} ->{Next?.Vertex?.Id:D3} f:{Face?.Id:D3}    | [{Opposite?.Id:D3}]  {Opposite?.Prev?.Vertex?.Id:D3}<- v:{Opposite?.Vertex?.Id:D3} ->{Opposite?.Next?.Vertex?.Id:D3} f:{Opposite?.Face?.Id:D3} ";
            //  return $"[{Id:D3}] {Prev?.Vertex?.PointId:D3}<- v:{Vertex?.PointId:D3} ->{Next?.Vertex?.PointId:D3} f:{Face?.Id:D3}    | [{Opposite?.Id:D3}]  {Opposite?.Prev?.Vertex?.PointId:D3}<- v:{Opposite?.Vertex?.PointId:D3} ->{Opposite?.Next?.Vertex?.PointId:D3} f:{Opposite?.Face?.Id:D3} ";
            string letra = IsContour ? "c" :
                IsOutside ? "o" :
                IsInterior ? "i" : "*";
             string r=  $"{letra}[{Id:D3}]  {PointId()}->{ Opposite?.PointId()} {VertexId()}->{Opposite?.VertexId()}";
            // string r = $"{PointId()}-{Opposite?.PointId()}";

            if (IsOutside) return $"[{r}]";
            return $" {r} ";

        }
        public IEnumerable<BasicHalfedge> AroundVerticeCCW(int max = 100000)
        {
            BasicHalfedge curr = this;
            //    BasicFace face = curr.Face;
            int i = 0;
            if (curr == null)
            {
                throw new InvalidOperationException($"{this} has not Halfedge asigned");
            }
            do
            {
                yield return curr;
                curr = curr.Opposite.Next;

                if (curr == null)
                {
                    throw new InvalidOperationException($"{this} has not Halfedge asigned");
                }
                //if (face != curr.Face)
                //{
                //    //              throw new InvalidOperationException($"{this} has not a only Face {curr}");
                //}
                if (max < i++)
                {
                    throw new GFLException($" {this} iterate more than {max} "); break;
                }
            } while (curr != this);
        }
        public IEnumerable<BasicHalfedge> AroundVerticeCW (int max = 100000)
        {
            BasicHalfedge curr = this;
            //    BasicFace face = curr.Face;
            int i = 0;
            if (curr == null)
            {
                throw new InvalidOperationException($"{this} has not Halfedge asigned");
            }
            do
            {
                yield return curr;
                curr = curr.Prev.Opposite;

                if (curr == null)
                {
                    throw new InvalidOperationException($"{this} has not Halfedge asigned");
                }
                //if (face != curr.Face)
                //{
                //    //              throw new InvalidOperationException($"{this} has not a only Face {curr}");
                //}
                if (max < i++)
                {
                    throw new GFLException($" {this} iterate more than {max} "); break;
                }
            } while (curr != this);
        }
        public IEnumerable<BasicHalfedge> Circulation(int max = 100000)
        {
            BasicHalfedge curr = this;
            //    BasicFace face = curr.Face;
            int i = 0;
            if (curr == null)
            {
                throw new InvalidOperationException($"{this} has not Halfedge asigned");
            }
            do
            {
                yield return curr;
                curr = curr.Next;

                if (curr == null)
                {
                    throw new InvalidOperationException($"{this} has not Halfedge asigned");
                }
                //if (face != curr.Face)
                //{
                //    //              throw new InvalidOperationException($"{this} has not a only Face {curr}");
                //}
                if (max < i++)
                {
                    throw new GFLException($" {this} iterate more than {max} "); break;
                }
            } while (curr != this);
        }

        public void FullLinkNext(BasicHalfedge other)
        {
            this.Next = other;
            this.Opposite.Prev = other.Opposite;
            this.Opposite.Vertex = other.Vertex;

            other.Prev = this;
            other.Opposite.Next = this.Opposite;
        }

        public BasicHalfedge LinkOpposite(BasicHalfedge other)
        {
            this.Opposite = other;
            other.Opposite = this;

            return this;
        }

        public void LinkNext(BasicHalfedge other)
        {
            this.Next = other;
            other.Prev = this;
        }

        public bool IsHole
        {
            get => Type.HasFlag(BasicHalfedgeType.IsHole);
            set
            {
                if (value)
                    Type |= BasicHalfedgeType.IsHole;
                else Type &= ~BasicHalfedgeType.IsHole;
            }
        }

        public bool IsContour
        {
            get => Type.HasFlag(BasicHalfedgeType.Countour);
            set
            {
                if (value)
                    Type |= BasicHalfedgeType.Countour;
                else
                    Type &= ~BasicHalfedgeType.Countour;
            }
        }

        public bool IsConstraint
        {
            get => Type.HasFlag(BasicHalfedgeType.Constraint);
            set
            {
                if (value)
                    Type |= BasicHalfedgeType.Constraint;
                else
                    Type &= ~BasicHalfedgeType.Constraint;
            }
        }

        public bool IsOutside
        {
            get => Type.HasFlag(BasicHalfedgeType.IsOutside);
            set
            {
                if (value)
                    Type |= BasicHalfedgeType.IsOutside;
                else
                    Type &= ~BasicHalfedgeType.IsOutside;
            }
        }

        public bool IsDeleted
        {
            get => Type.HasFlag(BasicHalfedgeType.Deleted);
            set
            {
                if (value)
                    Type |= BasicHalfedgeType.Deleted;
                else
                    Type &= ~BasicHalfedgeType.Deleted;
            }
        }
        
           

        public bool IsInterior
        {
            get => Type.HasFlag(BasicHalfedgeType.Interior);
            set
            {
                if (value)
                    Type |= BasicHalfedgeType.Interior;
                else
                    Type &= ~BasicHalfedgeType.Interior;
            }
        }
    }
}