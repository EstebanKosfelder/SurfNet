using TriangleNet.Topology.DCEL;


namespace SurfNet
{
    public class DcelFace
    {
        #region Static initialization of "Outer Space" face

        /// <summary>
        /// A face representing "outer space".
        /// </summary>
        public static readonly DcelFace Empty;

        static DcelFace()
        {
            Empty = new DcelFace(null);
            Empty.id = -1;
        }

        #endregion

        public int id;
        public int mark;

        // If the face is a Voronoi cell, this is the point that generates the cell.
        public DcelVertex generator;

        public DcelHalfEdge edge;
        public bool bounded;

        /// <summary>
        /// If part of a Voronoi diagram, returns the generator vertex
        /// of the face. Otherwise <c>null</c>.
        /// </summary>
        public DcelVertex Generator => generator;

        /// <summary>
        /// Gets or sets the face Id.
        /// </summary>
        public int Id
        {
            get { return id; }
            set { id = value; }
        }

        /// <summary>
        /// Gets or sets a half-edge connected to the face.
        /// </summary>
        public DcelHalfEdge Edge
        {
            get { return edge; }
            set { edge = value; }
        }

        /// <summary>
        /// Gets or sets a value, indicating if the face is bounded (for Voronoi diagram).
        /// </summary>
        public bool Bounded
        {
            get { return bounded; }
            set { bounded = value; }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Face" /> class.
        /// </summary>
        /// <param name="generator">The generator of this face (for Voronoi diagram)</param>
        public DcelFace(DcelVertex generator)
            : this(generator, null)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Face" /> class.
        /// </summary>
        /// <param name="generator">The generator of this face (for Voronoi diagram)</param>
        /// <param name="edge">The half-edge connected to this face.</param>
        public DcelFace(DcelVertex generator, DcelHalfEdge edge)
        {
            this.generator = generator;
            this.edge = edge;

            bounded = true;

            if (generator != null)
            {
                id = generator.ID;
            }
        }

        /// <summary>
        /// Enumerates all half-edges of the face boundary.
        /// </summary>
        /// <returns></returns>
        public IEnumerable<DcelHalfEdge> EnumerateEdges()
        {
            var edge = Edge;
            int first = edge.ID;

            do
            {
                yield return edge;

                edge = edge.Next;
            } while (edge.ID != first);
        }

        /// <inheritdoc />
        public override string ToString()
        {
            return string.Format("F-ID {0}", id);
        }

        internal bool is_constrained(int i)
        {
            throw new NotImplementedException();
        }

        internal DcelHalfEdge outer_ccbs_begin()
        {
            throw new NotImplementedException();
        }
    }

}
