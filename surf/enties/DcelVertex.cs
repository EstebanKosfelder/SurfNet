namespace SurfNet
{
    using System.Runtime.ConstrainedExecution;
    using static DebugLog;

    public class DcelVertex
    {
        public Point2 DcelPoint;
        internal int id;
        internal int label;
        internal DcelHalfEdge leaving;

        public int ID { get => id; set => id = value; }
        public int Label { get => label; set => label = value; }
        internal DcelHalfEdge Leaving { get => leaving; set => leaving = value; }

        public DcelVertex(double x, double y)
        {
            DcelPoint = new Point2(x, y);
        }
    }

    public enum Arr_parameter_space : byte
    {
        LEFT_BOUNDARY = 0,
        RIGHT_BOUNDARY,
        BOTTOM_BOUNDARY,
        TOP_BOUNDARY,
        INTERIOR,
        EXTERIOR,
        ARR_INTERIOR = INTERIOR,
        ARR_EXTERIOR = EXTERIOR
    };

    public class SkeletonDCELVertexBase
    {
        /*! \struct
         * An auxiliary structure for rebinding the vertex with a new point class.
         */
        //  template<typename PNT> struct rebind { typedef Arr_vertex_base<PNT> other; };

        protected object p_inc;  // An incident halfedge pointing at the vertex,
                                 // or the isolated vertex information (in case it is
                                 // isolated). The LSB of the pointer indicates whether
                                 // the vertex is isolated.
        protected Point_3 p_pt = Point_3.NaN;  // The point associated with the vertex.
        Arr_parameter_space[] pss = new Arr_parameter_space[2];  // The x and y parameter spaces (condensed in two bytes).

        /*! Default constructor. */
        private static int ctr;

        public SkeletonDCELVertexBase()
        {
            id = ctr++;
            pss[0] = pss[1] = Arr_parameter_space.ARR_INTERIOR;
        }

        
        public readonly int id;

        // Access/modification for pointer squatting
        public object inc() => p_inc;
        public void set_inc(object inc) { p_inc = inc; }

        /*! Check if the point pointer is null. */
        public bool has_null_point() =>p_pt.IsNaN(); 

        /*! Obtain the point (const version). */
        public Point_3 point()
        {
            assert(!p_pt.IsNaN( ));
            return (p_pt);
        }

        /*! Set the point (may be a null point). */
        public void set_point(Point_3 p) { p_pt = p; }

        /*! Obtain the boundary type in x. */
        public Arr_parameter_space parameter_space_in_x() { return pss[0]; }

        /*! Obtain the boundary type in y. */
        public Arr_parameter_space parameter_space_in_y() { return pss[1]; }

        /*! Set the boundary conditions of the vertex. */
        public void set_boundary(Arr_parameter_space ps_x, Arr_parameter_space ps_y)
        {
            pss[0] = ps_x;
            pss[1] = ps_y;
            return;
        }

        /*! Assign from another vertex. */
        public virtual void assign(SkeletonDCELVertexBase v)
        {
            p_pt = v.p_pt;
            pss[0] = v.pss[0];
            pss[1] = v.pss[1];
        }
    };
}