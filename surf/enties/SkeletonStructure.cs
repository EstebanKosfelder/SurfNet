namespace SurfNet
{
    using static DebugLog;




    public partial class SkeletonStructure
    {
        public DcelMesh input;
        public WavefrontEdgeList wavefront_edges;
        public KineticTriangulation kt;

        public WavefrontPropagator wp;

        public SkeletonStructure(DcelMesh input_ = null)
        {
            input = input_;
            wavefront_edges = new WavefrontEdgeList();
            kt = new KineticTriangulation();
            wp = new WavefrontPropagator(this);
        }

        /** copy vertices and edges from a Boost Graph. */

        public void initialize(int restrict_component = -1)
        {
            kt.initialize(input, wavefront_edges, restrict_component);
            assert(kt.triangles_size() > 0);
            wp.setup_queue(kt);
        }


        public DcelMesh get_input() => input;

        public KineticTriangulation get_kt() => kt;

        public SkeletonDCEL get_skeleton() => kt.get_skeleton();
    }
}