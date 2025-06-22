using TriangleNet.Topology.DCEL;


namespace SurfNet
{
    using static DebugLog;
    public class DcelHalfEdge
    {


        public enum HalfEdgeType
        {
            OutContour,
            InContour,
            Kinetic
        }

        public int id;
        public int mark;

        public DcelVertex origin;
        public DcelFace face;
        public DcelHalfEdge twin;
        public DcelHalfEdge next;

        /// <summary>
        /// Gets or sets the half-edge id.
        /// </summary>
        public int ID
        {
            get { return id; }
            set { id = value; }
        }

        /// <summary>
        /// Gets or sets a boundary marker.
        /// </summary>
        public int Boundary
        {
            get { return mark; }
            set { mark = value; }
        }

        /// <summary>
        /// Gets or sets the origin of the half-edge.
        /// </summary>
        public DcelVertex Origin
        {
            get { return origin; }
            set { origin = value; }
        }

        /// <summary>
        /// Gets or sets the face connected to the half-edge.
        /// </summary>
        public DcelFace Face
        {
            get { return face; }
            set { face = value; }
        }

        /// <summary>
        /// Gets or sets the twin of the half-edge.
        /// </summary>
        public DcelHalfEdge Twin
        {
            get { return twin; }
            set { twin = value; }
        }

        /// <summary>
        /// Gets or sets the next pointer of the half-edge.
        /// </summary>
        public DcelHalfEdge Next
        {
            get { return next; }
            set { next = value; }
        }

        public bool IsContour => this.Face == null || Face.Id < 0;

        /// <summary>
        /// Initializes a new instance of the <see cref="HalfEdge" /> class.
        /// </summary>
        /// <param name="origin">The origin of this half-edge.</param>
        public DcelHalfEdge(DcelVertex origin)
        {
            this.origin = origin;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="HalfEdge" /> class.
        /// </summary>
        /// <param name="origin">The origin of this half-edge.</param>
        /// <param name="face">The face connected to this half-edge.</param>
        public DcelHalfEdge(DcelVertex origin, DcelFace? face)
        {
            this.origin = origin;
            this.face = face;

            // IMPORTANT: do not remove the (face.edge == null) check!
            if (face != null && face.edge == null)
            {
                face.edge = this;
            }
        }

        /// <inheritdoc />
        public override string ToString()
        {
            return string.Format("HE-ID {0} (Origin = VID-{1})", id, origin.ID);
        }

        internal bool is_constrained()
        {
            return this.twin.Face == null || this.Twin.Face.Id < 0;
        }

        internal SkeletonDCELHalfedge? opposite()
        {
            throw new NotImplementedException();
        }
    }


    public abstract class Arr_halfedge_base<X_monotone_curve>
    {

        ///*! \struct
        // * An auxiliary structure for rebinding the halfedge with a new curve class.
        // */
        //template<typename XCV>
        //struct rebind { typedef Arr_halfedge_base<XCV> other; };


        protected Arr_halfedge_base<X_monotone_curve> p_opp;  // The opposite halfedge.
        protected Arr_halfedge_base<X_monotone_curve> p_prev; // The previous halfedge in the component boundary.
        protected Arr_halfedge_base<X_monotone_curve> p_next; // The next halfedge in the component boundary.

        protected DcelVertex p_v;    // The incident vertex (the target of the halfedge).
                                // The LSB of this pointer is used to store the
                                // direction of the halfedge.
        protected object p_comp; // The component this halfedge belongs to: the incident
                               // face for outer CCBs and the inner CCB information for
                               // inner CCBs. The LSB of the pointer indicates whether
                               // the halfedge lies on the boundary of an inner CCB.

        protected X_monotone_curve p_cv; // The associated x-monotone curve.


        /*! Default constructor */
        public  Arr_halfedge_base() 
        { }


        /*! Check if the curve pointer is null. */
        public bool has_null_curve() { return (p_cv == null); }

   

    /*! Obtain the x-monotone curve (non-const version). */
    public X_monotone_curve  curve()
  {
    assert(p_cv != null);
    return (p_cv);
}

    /*! Set the x-monotone curve. */
    public void set_curve(X_monotone_curve c)
{
    p_cv = c;

    // Set the curve for the opposite halfedge as well.
    Arr_halfedge_base<X_monotone_curve> opp = p_opp;

    opp.p_cv = c;
}

/*! Assign from another halfedge. */
public virtual void assign(Arr_halfedge_base<X_monotone_curve> he)
{ p_cv = he.p_cv; }
};


}
