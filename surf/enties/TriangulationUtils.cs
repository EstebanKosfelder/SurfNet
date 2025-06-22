namespace SurfNet
{
    public static class TriangulationUtils
    {
        public static readonly int[] _cw = { 2, 0, 1 };
        public static readonly int[] _ccw = { 1, 2, 0 };
        public static readonly int[] _mod3 = { 0, 1, 2, 0, 1 };

        public static int cw(this int i) => _cw[i];

        public static int ccw(this int i) => _ccw[i];

        public static int mod3(this int i) => _mod3[i];
    };
}
