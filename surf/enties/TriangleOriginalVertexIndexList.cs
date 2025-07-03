namespace SurfNet
{
    public class TriangleOriginalVertexIndexList : List<int>
    {
        public TriangleOriginalVertexIndexList()
        {
        }

        public TriangleOriginalVertexIndexList(IEnumerable<int> collection) : base(collection)
        {
        }

        public TriangleOriginalVertexIndexList(int capacity) : base(capacity)
        {
        }
    }
}
