namespace SurfNet
{
    using System.Runtime.CompilerServices;
    using static DebuggerInfo;
    using FT = double;
    using RT = double;
    using ST = int;
    using Vector_2 = Vector2;
    using We = double;

    public static partial class Mathex
    {
        public const FT EPS = 1e-9;
        public const FT EPS2 = EPS * EPS;
        public static FT PI => Math.PI;

        public const int SMALLER = -1;
        public const int EQUAL = 0;
        public const int LARGER = 1;

        public const int NEGATIVE = -1;
        public const int ZERO = 0;
        public const int POSITIVE = 1;


        public const double CORE_ZERO = 0.0;

        public const double FASTER_EDGE_WINS_IN_COLLINEAR_CASES = 1.0;



        public static bool are_near(FT a, FT b, FT eps = EPS)
        {
            if (!FT.IsRealNumber(a) || !FT.IsRealNumber(b)) throw new Exception();
            return (Math.Abs(a - b) < eps);
        }

        public static bool are_near(Point2 a, Point2 b, FT eps = EPS)
        {
            return are_near(a.X, b.X, eps) && are_near(a.Y, b.Y, eps);
        }

        public static bool AreNear(this FT a, FT b, FT eps = EPS)
        {
            if (!FT.IsRealNumber(a) || !FT.IsRealNumber(b)) throw new Exception();
            return (Math.Abs(a - b) < eps);
        }

        public static void Assert(bool v)
        {
            if (!v)
                throw new Exception();
        }

        public static void Assert(bool v, string message)
        {
            if (!v)
                throw new Exception(message);
        }

        public static int compare(FT a, FT b, FT eps = EPS)
        {
            return are_near(a, b, eps) ? 0 : a.CompareTo(b);
        }

        //CompareResultEnum?
        public static UncertainCompareResult compare_isec_anglesC2(in Vector_2 aBV1, in Vector_2 aBV2,
                                                  in Vector_2 aLV, in Vector_2 aRV)
        {
            UncertainCompareResult rResult = UncertainCompareResult.indeterminate;
            Vector_2 lBisectorDirection = aBV2 - aBV1;
            FT lLNorm = inexact_sqrt(compute_scalar_product_2(aLV, aLV));
            FT lRNorm = inexact_sqrt(compute_scalar_product_2(aRV, aRV));

            if ((bool)(Not(certified_is_positive(lLNorm))
                   .Or(Not(certified_is_positive(lRNorm)))))
                return rResult;

            var aLVn = aLV / lLNorm;
            var aRVn = aRV / lRNorm;

            FT lLSp = compute_scalar_product_2(lBisectorDirection, aLVn);
            FT lRSp = compute_scalar_product_2(lBisectorDirection, aRVn);

            // Smaller if the scalar product is larger, so swapping
            rResult = certified_compare(lRSp, lLSp);

            return rResult;
        }

        public static UncertainCompareResult compare_ss_event_angles_2(in Vector_2 aBV1, in Vector_2 aBV2,
                                                  in Vector_2 aLV, in Vector_2 aRV)
        {
            return compare_isec_anglesC2(aBV1, aBV2, aLV, aRV);
        }

        public static int compare_x(Point2 p, Point2 q, FT eps = EPS)
        {
            return compare(p.X, q.X, eps);
        }

        public static int compare_xy(Point2 p, Point2 q, FT eps = EPS)
        {
            var phx = p.hx();
            var phy = p.hy();
            var phw = p.hw();
            var qhx = q.hx();
            var qhy = q.hy();
            var qhw = q.hw();

            var pV = phx * qhw;
            var qV = qhx * phw;
            if (are_near(pV, qV))
            {
                pV = phy * qhw;
                qV = qhy * phw;
            }
            return compare(pV, qV);
        }

        public static int compare_y(Point2 p, Point2 q, FT eps = EPS)
        {
            return compare(p.Y, q.Y, eps);
        }

        public static Polynomial1D compute_determinant(Polynomial1D x0, Polynomial1D y0,
                          Polynomial1D x1, Polynomial1D y1,
                          Polynomial1D x2, Polynomial1D y2)
        {
            var result = (
              (x0 * y1
              + x1 * y2
              + x2 * y0
              )
              -
              (y0 * x1
              + y1 * x2
              + y2 * x0
              )
            );
            return result;
        }

        public static double compute_determinant(double x0, double y0,
                         double x1, double y1,
                         double x2, double y2)
        {
            var result = (
              (x0 * y1
              + x1 * y2
              + x2 * y0
              )
              -
              (y0 * x1
              + y1 * x2
              + y2 * x0
              )
            );
            return result;
        }

        public static FT compute_scalar_product_2(Vector2 a, Vector2 b) => a.X * b.X + a.Y * b.Y;

        public static Vector2 construct_opposite_vector_2(Vector2 v) => -v;

        public static Direction2 construct_perpendicular_direction_2(Direction2 d, OrientationEnum o)
        {
            DebuggerInfo.CGAL_kernel_precondition(o != OrientationEnum.COLLINEAR);
            if (o == OrientationEnum.COUNTERCLOCKWISE)
                return new Direction2(-d.dy(), d.dx());
            else
                return new Direction2(d.dy(), -d.dx());
        }

        public static FT determinant(FT a00, FT a01, FT a10, FT a11)
        {
            // First compute the det2x2
            FT m01 = a00 * a11 - a10 * a01;
            return m01;
        }

        public static FT determinant(FT a00, FT a01, FT a02, FT a10, FT a11, FT a12, FT a20, FT a21, FT a22)
        {
            // First compute the det2x2
            FT m01 = a00 * a11 - a10 * a01;
            FT m02 = a00 * a21 - a20 * a01;
            FT m12 = a10 * a21 - a20 * a11;
            // Now compute the minors of rank 3
            FT m012 = m01 * a22 - m02 * a12 + m12 * a02;
            return m012;
        }

        public static FT determinant(FT a00, FT a01, FT a02, FT a03, FT a10, FT a11, FT a12, FT a13, FT a20, FT a21, FT a22, FT a23, FT a30, FT a31, FT a32, FT a33)
        {
            // First compute the det2x2
            FT m01 = a10 * a01 - a00 * a11;
            FT m02 = a20 * a01 - a00 * a21;
            FT m03 = a30 * a01 - a00 * a31;
            FT m12 = a20 * a11 - a10 * a21;
            FT m13 = a30 * a11 - a10 * a31;
            FT m23 = a30 * a21 - a20 * a31;
            // Now compute the minors of rank 3
            FT m012 = m12 * a02 - m02 * a12 + m01 * a22;
            FT m013 = m13 * a02 - m03 * a12 + m01 * a32;
            FT m023 = m23 * a02 - m03 * a22 + m02 * a32;
            FT m123 = m23 * a12 - m13 * a22 + m12 * a32;
            // Now compute the minors of rank 4
            FT m0123 = m123 * a03 - m023 * a13 + m013 * a23 - m012 * a33;
            return m0123;
        }

        public static FT determinant(FT a00, FT a01, FT a02, FT a03, FT a04, FT a10, FT a11, FT a12, FT a13, FT a14, FT a20, FT a21, FT a22, FT a23, FT a24, FT a30, FT a31, FT a32, FT a33, FT a34, FT a40, FT a41, FT a42, FT a43, FT a44)
        {
            // First compute the det2x2
            FT m01 = a10 * a01 - a00 * a11;
            FT m02 = a20 * a01 - a00 * a21;
            FT m03 = a30 * a01 - a00 * a31;
            FT m04 = a40 * a01 - a00 * a41;
            FT m12 = a20 * a11 - a10 * a21;
            FT m13 = a30 * a11 - a10 * a31;
            FT m14 = a40 * a11 - a10 * a41;
            FT m23 = a30 * a21 - a20 * a31;
            FT m24 = a40 * a21 - a20 * a41;
            FT m34 = a40 * a31 - a30 * a41;
            // Now compute the minors of rank 3
            FT m012 = m12 * a02 - m02 * a12 + m01 * a22;
            FT m013 = m13 * a02 - m03 * a12 + m01 * a32;
            FT m014 = m14 * a02 - m04 * a12 + m01 * a42;
            FT m023 = m23 * a02 - m03 * a22 + m02 * a32;
            FT m024 = m24 * a02 - m04 * a22 + m02 * a42;
            FT m034 = m34 * a02 - m04 * a32 + m03 * a42;
            FT m123 = m23 * a12 - m13 * a22 + m12 * a32;
            FT m124 = m24 * a12 - m14 * a22 + m12 * a42;
            FT m134 = m34 * a12 - m14 * a32 + m13 * a42;
            FT m234 = m34 * a22 - m24 * a32 + m23 * a42;
            // Now compute the minors of rank 4
            FT m0123 = m123 * a03 - m023 * a13 + m013 * a23 - m012 * a33;
            FT m0124 = m124 * a03 - m024 * a13 + m014 * a23 - m012 * a43;
            FT m0134 = m134 * a03 - m034 * a13 + m014 * a33 - m013 * a43;
            FT m0234 = m234 * a03 - m034 * a23 + m024 * a33 - m023 * a43;
            FT m1234 = m234 * a13 - m134 * a23 + m124 * a33 - m123 * a43;
            // Now compute the minors of rank 5
            FT m01234 = m1234 * a04 - m0234 * a14 + m0134 * a24 - m0124 * a34 + m0123 * a44;
            return m01234;
        }

        public static FT determinant(FT a00, FT a01, FT a02, FT a03, FT a04, FT a05, FT a10, FT a11, FT a12, FT a13, FT a14, FT a15, FT a20, FT a21, FT a22, FT a23, FT a24, FT a25, FT a30, FT a31, FT a32, FT a33, FT a34, FT a35, FT a40, FT a41, FT a42, FT a43, FT a44, FT a45, FT a50, FT a51, FT a52, FT a53, FT a54, FT a55)
        {
            // First compute the det2x2
            FT m01 = a00 * a11 - a10 * a01;
            FT m02 = a00 * a21 - a20 * a01;
            FT m03 = a00 * a31 - a30 * a01;
            FT m04 = a00 * a41 - a40 * a01;
            FT m05 = a00 * a51 - a50 * a01;
            FT m12 = a10 * a21 - a20 * a11;
            FT m13 = a10 * a31 - a30 * a11;
            FT m14 = a10 * a41 - a40 * a11;
            FT m15 = a10 * a51 - a50 * a11;
            FT m23 = a20 * a31 - a30 * a21;
            FT m24 = a20 * a41 - a40 * a21;
            FT m25 = a20 * a51 - a50 * a21;
            FT m34 = a30 * a41 - a40 * a31;
            FT m35 = a30 * a51 - a50 * a31;
            FT m45 = a40 * a51 - a50 * a41;
            // Now compute the minors of rank 3
            FT m012 = m01 * a22 - m02 * a12 + m12 * a02;
            FT m013 = m01 * a32 - m03 * a12 + m13 * a02;
            FT m014 = m01 * a42 - m04 * a12 + m14 * a02;
            FT m015 = m01 * a52 - m05 * a12 + m15 * a02;
            FT m023 = m02 * a32 - m03 * a22 + m23 * a02;
            FT m024 = m02 * a42 - m04 * a22 + m24 * a02;
            FT m025 = m02 * a52 - m05 * a22 + m25 * a02;
            FT m034 = m03 * a42 - m04 * a32 + m34 * a02;
            FT m035 = m03 * a52 - m05 * a32 + m35 * a02;
            FT m045 = m04 * a52 - m05 * a42 + m45 * a02;
            FT m123 = m12 * a32 - m13 * a22 + m23 * a12;
            FT m124 = m12 * a42 - m14 * a22 + m24 * a12;
            FT m125 = m12 * a52 - m15 * a22 + m25 * a12;
            FT m134 = m13 * a42 - m14 * a32 + m34 * a12;
            FT m135 = m13 * a52 - m15 * a32 + m35 * a12;
            FT m145 = m14 * a52 - m15 * a42 + m45 * a12;
            FT m234 = m23 * a42 - m24 * a32 + m34 * a22;
            FT m235 = m23 * a52 - m25 * a32 + m35 * a22;
            FT m245 = m24 * a52 - m25 * a42 + m45 * a22;
            FT m345 = m34 * a52 - m35 * a42 + m45 * a32;
            // Now compute the minors of rank 4
            FT m0123 = m012 * a33 - m013 * a23 + m023 * a13 - m123 * a03;
            FT m0124 = m012 * a43 - m014 * a23 + m024 * a13 - m124 * a03;
            FT m0125 = m012 * a53 - m015 * a23 + m025 * a13 - m125 * a03;
            FT m0134 = m013 * a43 - m014 * a33 + m034 * a13 - m134 * a03;
            FT m0135 = m013 * a53 - m015 * a33 + m035 * a13 - m135 * a03;
            FT m0145 = m014 * a53 - m015 * a43 + m045 * a13 - m145 * a03;
            FT m0234 = m023 * a43 - m024 * a33 + m034 * a23 - m234 * a03;
            FT m0235 = m023 * a53 - m025 * a33 + m035 * a23 - m235 * a03;
            FT m0245 = m024 * a53 - m025 * a43 + m045 * a23 - m245 * a03;
            FT m0345 = m034 * a53 - m035 * a43 + m045 * a33 - m345 * a03;
            FT m1234 = m123 * a43 - m124 * a33 + m134 * a23 - m234 * a13;
            FT m1235 = m123 * a53 - m125 * a33 + m135 * a23 - m235 * a13;
            FT m1245 = m124 * a53 - m125 * a43 + m145 * a23 - m245 * a13;
            FT m1345 = m134 * a53 - m135 * a43 + m145 * a33 - m345 * a13;
            FT m2345 = m234 * a53 - m235 * a43 + m245 * a33 - m345 * a23;
            // Now compute the minors of rank 5
            FT m01234 = m0123 * a44 - m0124 * a34 + m0134 * a24 - m0234 * a14 + m1234 * a04;
            FT m01235 = m0123 * a54 - m0125 * a34 + m0135 * a24 - m0235 * a14 + m1235 * a04;
            FT m01245 = m0124 * a54 - m0125 * a44 + m0145 * a24 - m0245 * a14 + m1245 * a04;
            FT m01345 = m0134 * a54 - m0135 * a44 + m0145 * a34 - m0345 * a14 + m1345 * a04;
            FT m02345 = m0234 * a54 - m0235 * a44 + m0245 * a34 - m0345 * a24 + m2345 * a04;
            FT m12345 = m1234 * a54 - m1235 * a44 + m1245 * a34 - m1345 * a24 + m2345 * a14;
            // Now compute the minors of rank 6
            FT m012345 = m01234 * a55 - m01235 * a45 + m01245 * a35 - m01345 * a25
                             + m02345 * a15 - m12345 * a05;
            return m012345;
        }

        public static (int, int, int) indirect_sort_3<T>(T[] t) where T : IComparable
        {
            int i0 = 0;
            int i1 = 1;
            int i2 = 2;
            if (t[i0].CompareTo(t[i1]) > 0) swap(ref i0, ref i1);
            if (t[i1].CompareTo(t[i2]) > 0) swap(ref i1, ref i2);
            if (t[i0].CompareTo(t[i1]) > 0) swap(ref i0, ref i1);
            return (i0, i1, i2);
        }
        public static Point2 validate(Point2? value)
        {
            if (value != null) return value.Value;
            throw new ArithmeticOverflowException();
        }

        public static T validate<T>(T? v) where T : class
        {
            if (v == null) throw new NullReferenceException();
            return v;
        }

        public static T validate<T>(T? v) where T : struct
        {
            if (!v.HasValue) throw new ArithmeticOverflowException();
            return v.Value;
        }

        public static bool handle_assigned<T>(T? value)
        {
            return value != null;
        }


        public static bool is_negative(FT x) => sign(x) == -1;

        public static bool is_one(FT x) => are_near(x, 1.0);

        public static bool is_zero(this FT value, double eps = EPS)
        {
            return Math.Abs(value) < eps;
        }

        public static Point2 Min(Point2 a, Point2 b)
        {
            return new Point2(Math.Min(a.X, b.X), Math.Min(a.Y, b.Y));
        }

        public static Point2 Max(Point2 a, Point2 b)
        {
            return new Point2(Math.Max(a.X, b.X), Math.Max(a.Y, b.Y));
        }
        public static FT inexact_sqrt(FT a) => Math.Sqrt(a);

        public static (Point2 c, double r) InnerCircle(params Point2[] v)
        {
            var a = distance(v[1], v[0]);
            var b = distance(v[2], v[1]);
            var c = distance(v[0], v[2]);

            //var t = (v[0] + v[1] + v[2]);
            //var center = t / 3.0;
            // Calculate the semiperimeter
            var s = (a + b + c) / 2.0;

            // Calculate the area using Heron's formula
            var area = Math.Sqrt(s * (s - a) * (s - b) * (s - c));

            // Calculate the inradius
            var inradius = area / s;

            // Calculate the incenter using weighted average of vertices
            var center = (v[2] * a + v[0] * b + v[1] * c) / (a + b + c);

            return (center, inradius);
        }

        public static Intersection intersection(Ray2 ray, Segment2 seg) => new Ray2Segment2Intersection(ray, seg);

        public static Intersection intersection(Line2 a, Line2 b) => new Line_2_Line_2_pair(a, b);

        public static bool is_finite(FT a) => FT.IsFinite(a);

        public static bool is_positive(FT a) => sign(a) == 1;

        public static bool IsFinite(this FT a) => is_finite(a);

        public static bool IsZero(this FT a, double eps = EPS) => is_zero(a, eps);

        public static Point2 midpoint(in Point2 a, in Point2 b) => (a + b) * 0.5;

        public static int orientation(Point2 u, Point2 v)
        {
            return sign_of_determinant(u.X, u.Y, v.X, v.Y);
        }

        public static OrientationEnum orientation(Vector2 u, Vector2 v)
        {
            return (OrientationEnum)sign_of_determinant(u.X, u.Y, v.X, v.Y);
        }

        public static int oriented(FT a, FT eps = EPS)
        {
            if (are_near(a, 0, eps)) return 1;
            return Math.Sign(a);
        }

        //    (OrientationEnum)sign_of_determinant(q.X - p.X, q.Y - p.Y, r.X - p.X, r.Y - p.Y);
        public static void Postcondition(bool v)
        {
            if (!v)
                throw new Exception(); ;
        }

        //public static OrientationEnum orientation(in Point2 p, in Point2 q, in Point2 r) =>
        public static void Precondition_msg(bool v1, string v2)
        {
            Assert(v1, v2);
        }

        public static double Proporcional(this double from, double to, double time)
        {
            return (to - from) * time + from; // is de same (1 - time) * from + time * to;
        }

        public static int side_of(FT a, FT eps = EPS)
        {
            if (are_near(a, 0, eps)) return 0;
            return Math.Sign(a);
        }

        public static int sign(FT vaule, FT eps = EPS) => is_zero(vaule) ? 0 : Math.Sign(vaule);

        public static int sign_of_determinant(FT a00, FT a01, FT a10, FT a11)
        {
            return compare((a00 * a11), (a10 * a01));
        }

        public static int sign_of_determinant(FT a00, FT a01, FT a02, FT a10, FT a11, FT a12, FT a20, FT a21, FT a22)
        {
            return sign(determinant(a00, a01, a02, a10, a11, a12, a20, a21, a22));
        }

        public static int sign_of_determinant(FT a00, FT a01, FT a02, FT a03, FT a10, FT a11, FT a12, FT a13, FT a20, FT a21, FT a22, FT a23, FT a30, FT a31, FT a32, FT a33)
        {
            return sign(determinant(a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33));
        }

        public static int sign_of_determinant(FT a00, FT a01, FT a02, FT a03, FT a04, FT a10, FT a11, FT a12, FT a13, FT a14, FT a20, FT a21, FT a22, FT a23, FT a24, FT a30, FT a31, FT a32, FT a33, FT a34, FT a40, FT a41, FT a42, FT a43, FT a44)
        {
            return sign(determinant(a00, a01, a02, a03, a04,
                                                   a10, a11, a12, a13, a14,
                                                   a20, a21, a22, a23, a24,
                                                   a30, a31, a32, a33, a34,
                                                   a40, a41, a42, a43, a44));
        }

        public static int sign_of_determinant(FT a00, FT a01, FT a02, FT a03, FT a04, FT a05, FT a10, FT a11, FT a12, FT a13, FT a14, FT a15, FT a20, FT a21, FT a22, FT a23, FT a24, FT a25, FT a30, FT a31, FT a32, FT a33, FT a34, FT a35, FT a40, FT a41, FT a42, FT a43, FT a44, FT a45, FT a50, FT a51, FT a52, FT a53, FT a54, FT a55)
        {
            return sign(determinant(a00, a01, a02, a03, a04, a05, a10, a11, a12, a13, a14, a15, a20, a21, a22, a23, a24, a25, a30, a31, a32, a33, a34, a35, a40, a41, a42, a43, a44, a45, a50, a51, a52, a53, a54, a55));
        }

        public static (bool has_real_roots, bool is_square) solve_quadratic(Polynomial1D f, out double x0, out double x1) => f.solve_quadratic(out x0, out x1);

        public static FT square(FT value) => value * value;

        public static void swap<T>(ref T a, ref T b)
        { (b, a) = (a, b); }

        public static FT validate(FT value)
        {
            if (!is_finite(value))
                throw new ArithmeticOverflowException();
            return value;
        }

        public static Point2 validate(Point2 value)
        {
            validate(value.X);
            validate(value.Y);
            return value;
        }

        public static Vector2 validate(Vector2 value)
        {
            validate(value.X);
            validate(value.Y);
            return value;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void
       midpoint(FT px, FT py,
                   FT qx, FT qy,
                   out FT x, out FT y)
        {
            x = (px + qx) * 0.5;
            y = (py + qy) * 0.5;
        }

        ////CGAL_KERNEL_LARGE_INLINE
        public static void
        circumcenter_translate(FT dqx, FT dqy,
                                 FT drx, FT dry,
                                       out FT dcx, out FT dcy)
        {
            // Given 3 points P, Q, R, this function takes as input:
            // qx-px, qy-py, rx-px, ry-py.  And returns cx-px, cy-py,
            // where (cx, cy) are the coordinates of the circumcenter C.

            // What we do is intersect the bisectors.
            FT r2 = square(drx) + square(dry);
            FT q2 = square(dqx) + square(dqy);
            FT den = 2 * determinant(dqx, dqy, drx, dry);

            // The 3 points aren't collinear.
            // Hopefully, this is already checked at the upper level.
            CGAL_kernel_assertion(!is_zero(den));

            // One possible optimization here is to precompute 1/den, to avoid one
            // division.  However, we lose precision, and it's maybe not worth it (?).
            dcx = determinant(dry, dqy, r2, q2) / den;
            dcy = -determinant(drx, dqx, r2, q2) / den;
        }

        //
        public static void
         circumcenter(FT px, FT py,
                         FT qx, FT qy,
                         FT rx, FT ry,
                         out FT x, out FT y)
        {
            circumcenter_translate(qx - px, qy - py, rx - px, ry - py, out x, out y);
            x += px;
            y += py;
        }

        public static void
        barycenter(FT p1x, FT p1y, FT w1,
                     FT p2x, FT p2y,
                     out FT x, out FT y)
        {
            FT w2 = 1 - w1;
            x = w1 * p1x + w2 * p2x;
            y = w1 * p1y + w2 * p2y;
        }

        public static void
        barycenter(FT p1x, FT p1y, FT w1,
                     FT p2x, FT p2y, FT w2,
                     out FT x, out FT y)
        {
            FT sum = w1 + w2;
            CGAL_kernel_assertion(sum != 0);
            x = (w1 * p1x + w2 * p2x) / sum;
            y = (w1 * p1y + w2 * p2y) / sum;
        }

        public static void
        barycenter(FT p1x, FT p1y, FT w1,
                     FT p2x, FT p2y, FT w2,
                     FT p3x, FT p3y,
                     out FT x, out FT y)
        {
            FT w3 = 1 - w1 - w2;
            x = w1 * p1x + w2 * p2x + w3 * p3x;
            y = w1 * p1y + w2 * p2y + w3 * p3y;
        }

        public static void
        barycenter(FT p1x, FT p1y, FT w1,
                     FT p2x, FT p2y, FT w2,
                     FT p3x, FT p3y, FT w3,
                     out FT x, out FT y)
        {
            FT sum = w1 + w2 + w3;
            CGAL_kernel_assertion(sum != 0);
            x = (w1 * p1x + w2 * p2x + w3 * p3x) / sum;
            y = (w1 * p1y + w2 * p2y + w3 * p3y) / sum;
        }

        public static void
        barycenter(FT p1x, FT p1y, FT w1,
                     FT p2x, FT p2y, FT w2,
                     FT p3x, FT p3y, FT w3,
                     FT p4x, FT p4y,
                     out FT x, out FT y)
        {
            FT w4 = 1 - w1 - w2 - w3;
            x = w1 * p1x + w2 * p2x + w3 * p3x + w4 * p4x;
            y = w1 * p1y + w2 * p2y + w3 * p3y + w4 * p4y;
        }

        public static void
        barycenter(FT p1x, FT p1y, FT w1,
                     FT p2x, FT p2y, FT w2,
                     FT p3x, FT p3y, FT w3,
                     FT p4x, FT p4y, FT w4,
                     out FT x, out FT y)
        {
            FT sum = w1 + w2 + w3 + w4;
            CGAL_kernel_assertion(sum != 0);
            x = (w1 * p1x + w2 * p2x + w3 * p3x + w4 * p4x) / sum;
            y = (w1 * p1y + w2 * p2y + w3 * p3y + w4 * p4y) / sum;
        }

        //
        public static void
        centroid(FT px, FT py,
                    FT qx, FT qy,
                    FT rx, FT ry,
                    out FT x, out FT y)
        {
            x = (px + qx + rx) / 3;
            y = (py + qy + ry) / 3;
        }

        //
        public static void
        centroid(FT px, FT py,
                    FT qx, FT qy,
                    FT rx, FT ry,
                    FT sx, FT sy,
                    out FT x, out FT y)
        {
            x = (px + qx + rx + sx) / 4;
            y = (py + qy + ry + sy) / 4;
        }

        //    [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void
           line_from_points(FT px, FT py,
                              FT qx, FT qy,
                              out FT a, out FT b, out FT c)
        {
            // The horizontal and vertical line get a special treatment
            // in order to make the intersection code robust for doubles
            if (py == qy)
            {
                a = 0;
                if (qx > px)
                {
                    b = 1;
                    c = -py;
                }
                else if (qx == px)
                {
                    b = 0;
                    c = 0;
                }
                else
                {
                    b = -1;
                    c = py;
                }
            }
            else if (qx == px)
            {
                b = 0;
                if (qy > py)
                {
                    a = -1;
                    c = px;
                }
                else if (qy == py)
                {
                    a = 0;
                    c = 0;
                }
                else
                {
                    a = 1;
                    c = -px;
                }
            }
            else
            {
                a = py - qy;
                b = qx - px;
                c = -px * a - py * b;
            }
        }

        //   [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void
        line_from_point_direction(FT px, FT py,
                                    FT dx, FT dy,
                                    out FT a, out FT b, out FT c)
        {
            a = -dy;
            b = dx;
            c = px * dy - py * dx;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void
        bisector_of_points(FT px, FT py,
                             FT qx, FT qy,
                             out FT a, out FT b, out FT c)
        {
            a = 2 * (px - qx);
            b = 2 * (py - qy);
            c = square(qx) + square(qy) -
                square(px) - square(py);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void
        bisector_of_lines(FT pa, FT pb, FT pc,
                            FT qa, FT qb, FT qc,
                            out FT a, out FT b, out FT c)
        {
            // We normalize the equations of the 2 lines, and we then add them.
            FT n1 = Math.Sqrt((square(pa) + square(pb)));
            FT n2 = Math.Sqrt((square(qa) + square(qb)));
            a = n2 * pa + n1 * qa;
            b = n2 * pb + n1 * qb;
            c = n2 * pc + n1 * qc;

            // Care must be taken for the case when this produces a degenerate line.
            if (a == 0 && b == 0)
            {
                a = n2 * pa - n1 * qa;
                b = n2 * pb - n1 * qb;
                c = n2 * pc - n1 * qc;
            }
        }

        //   [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT
        line_y_at_x(FT a, FT b, FT c, FT x)
        {
            return (-a * x - c) / b;
        }

        //   [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void
        line_get_point(FT a, FT b, FT c, FT i,
                         out FT x, out FT y)
        {
            if (is_zero(b))
            {
                x = -c / a;
                y = 1 - i * a;
            }
            else
            {
                x = 1 + i * b;
                y = -(a + c) / b - i * a;
            }
        }

        public static void
            perpendicular_through_point(FT la, FT lb,
                                          FT px, FT py,
                                          out FT a, out FT b, out FT c)
        {
            a = -lb;
            b = la;
            c = lb * px - la * py;
        }

        //
        public static void
        line_project_point(FT la, FT lb, FT lc,
                             FT px, FT py,
                             out FT x, out FT y)
        {
            if ((is_zero(la))) // horizontal line
            {
                x = px;
                y = -lc / lb;
            }
            else if ((is_zero(lb))) // vertical line
            {
                x = -lc / la;
                y = py;
            }
            else
            {
                FT a2 = square(la);
                FT b2 = square(lb);
                FT d = a2 + b2;
                x = (b2 * px - la * lb * py - la * lc) / d;
                y = (-la * lb * px + a2 * py - lb * lc) / d;
            }
        }

        //
        public static FT
        squared_radius(FT px, FT py,
                         FT qx, FT qy,
                         FT rx, FT ry,
                         out FT x, out FT y)
        {
            circumcenter_translate(qx - px, qy - py, rx - px, ry - py, out x, out y);
            FT r2 = square(x) + square(y);
            x += px;
            y += py;
            return r2;
        }

        //
        public static FT
        squared_radius(FT px, FT py,
                         FT qx, FT qy,
                         FT rx, FT ry)
        {
            FT x, y;
            circumcenter_translate(qx - px, qy - py, rx - px, ry - py, out x, out y);
            return square(x) + square(y);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT squared_distance(FT px, FT py, FT qx, FT qy)
        {
            return square(px - qx) + square(py - qy);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT squared_distance(in Point2 p, in Point2 q)
        {
            return squared_distance(p.X, p.Y, q.X, q.Y);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT distance(in Point2 p, in Point2 q)
        {
            return Math.Sqrt(squared_distance(p.X, p.Y, q.X, q.Y));
        }

        //    [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT squared_radius(FT px, FT py,
                     FT qx, FT qy)
        {
            return squared_distance(px, py, qx, qy) / 4;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT scaled_distance_to_line(FT la, FT lb, FT lc,
                                   FT px, FT py)
        {
            // for comparisons, use distance_to_directionsC2 instead
            // since lc is irrelevant
            return la * px + lb * py + lc;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT
        scaled_distance_to_direction(FT la, FT lb,
                                        FT px, FT py)
        {
            // scalar product with direction
            return la * px + lb * py;
        }

        //
        public static FT
        scaled_distance_to_line(FT px, FT py,
                                   FT qx, FT qy,
                                   FT rx, FT ry)
        {
            return determinant(px - rx, py - ry, qx - rx, qy - ry);
        }

        public static void
         weighted_circumcenter_translate(RT dqx, RT dqy, RT dqw,
                                           RT drx, RT dry, RT drw,
                                           out RT dcx, out RT dcy)
        {
            // Given 3 points P, Q, R, this function takes as input:
            // qx-px, qy-py,qw-pw,  rx-px, ry-py, rw-pw.  And returns cx-px, cy-py,
            // where (cx, cy) are the coordinates of the circumcenter C.

            // What we do is intersect the radical axis
            RT r2 = square(drx) + square(dry) - drw;
            RT q2 = square(dqx) + square(dqy) - dqw;

            RT den = (2) * determinant(dqx, dqy, drx, dry);

            // The 3 points aren't collinear.
            // Hopefully, this is already checked at the upper level.
            CGAL_assertion(den != 0);

            // One possible optimization here is to precompute 1/den, to avoid one
            // division.  However, we lose precision, and it's maybe not worth it (?).
            dcx = determinant(dry, dqy, r2, q2) / den;
            dcy = -determinant(drx, dqx, r2, q2) / den;
        }

        /*
        //template < class RT >
        public static void
        weighted_circumcenter( RT px, RT py, We pw,
                                 RT qx, RT qy, We qw,
                                 RT rx, RT ry, We rw,
                                 out RT x, out RT y )
        {
          RT dqw = (RT)(qw-pw);
          RT drw = (RT)(rw-pw);

          weighted_circumcenter_translateC2<RT>(qx-px, qy-py, dqw,rx-px, ry-py,drw,x, y);
          x += px;
          y += py;
        }
        */

        public static FT
        power_product(FT px, FT py, FT pw,
                        FT qx, FT qy, FT qw)
        {
            // computes the power product of two weighted points
            FT qpx = qx - px;
            FT qpy = qy - py;
            FT qp2 = square(qpx) + square(qpy);
            return qp2 - pw - qw;
        }

        /*
        public static void
        radical_axis(RT px, RT py,  We pw,
                       RT qx, RT qy, We qw,
                       out RT a, out RT b, out RT c )
        {
          a =  (RT)(2)*(px - qx);
          b =  (RT)(2)*(py - qy);
          c = - square(px) - square(py)
              + square(qx) + square(qy)
              + (RT)(pw) - (RT)(qw);
        }

        */

        public static FT
        squared_radius_orthogonal_circle(FT px, FT py, FT pw,
                                           FT qx, FT qy, FT qw,
                                           FT rx, FT ry, FT rw)
        {
            FT FT4 = (4);
            FT dpx = px - rx;
            FT dpy = py - ry;
            FT dqx = qx - rx;
            FT dqy = qy - ry;
            FT dpp = square(dpx) + square(dpy) - pw + rw;
            FT dqq = square(dqx) + square(dqy) - qw + rw;

            FT det0 = determinant(dpx, dpy, dqx, dqy);
            FT det1 = determinant(dpp, dpy, dqq, dqy);
            FT det2 = determinant(dpx, dpp, dqx, dqq);

            return (square(det1) + square(det2)) /
                                            (FT4 * square(det0)) - rw;
        }

        //
        public static FT
        squared_radius_smallest_orthogonal_circle(FT px, FT py, FT pw,
                                                    FT qx, FT qy, FT qw)
        {
            FT FT4 = (4);
            FT dpz = square(px - qx) + square(py - qy);
            return (square(dpz - pw + qw) / (FT4 * dpz) - qw);
        }
    }

    public static partial class Mathex
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool parallel(FT l1a, FT l1b, FT l2a, FT l2b)
        {
            return sign_of_determinant(l1a, l1b, l2a, l2b) == 0;
        }

        public static bool parallel(FT s1sx, FT s1sy, FT s1tx, FT s1ty, FT s2sx, FT s2sy, FT s2tx, FT s2ty)
        {
            return sign_of_determinant(s1tx - s1sx, s1ty - s1sy, s2tx - s2sx, s2ty - s2sy) == 0;
        }

        public static bool equal_line(FT l1a, FT l1b, FT l1c,
             FT l2a, FT l2b, FT l2c)
        {
            if (sign_of_determinant(l1a, l1b, l2a, l2b) != 0)
                return false; // Not parallel.
            ST s1a = sign(l1a);
            if (s1a != 0)
                return s1a == sign(l2a) && sign_of_determinant(l1a, l1c, l2a, l2c) == 0;
            return sign(l1b) == sign(l2b) && sign_of_determinant(l1b, l1c, l2b, l2c) == 0;
        }

        public static ST
        compare_x(FT px,
            FT la, FT lb, FT lc,
            FT ha, FT hb, FT hc)
        {
            // The abscissa of the intersection point is num/den.
            FT num = determinant(lb, lc, hb, hc);
            FT den = determinant(la, lb, ha, hb);
            ST s = sign(den);
            CGAL_kernel_assertion(s != 0);
            return s * compare(px * den, num);
        }

        public static ST
        compare_x(FT la, FT lb, FT lc,
            FT h1a, FT h1b, FT h1c,
            FT h2a, FT h2b, FT h2c)
        {
            /*
            FT num1 = determinant( lb, lc, h1b, h1c);
            FT den1 = determinant( la, lb, h1a, h1b);
            FT num2 = determinant( lb, lc, h2b, h2c);
            FT den2 = determinant( la, lb, h2a, h2b);
            Sign s = Sign (sign(den1) * sign(den2));
            CGAL_kernel_assertion( s != 0 );
            return s * sign_of_determinant(num1, num2, den1, den2);
            */
            FT num1 = determinant(la, lc, h1a, h1c);
            FT num2 = determinant(la, lc, h2a, h2c);
            FT num = determinant(h1a, h1c, h2a, h2c) * lb
                      + determinant(num1, num2, h1b, h2b);
            FT den1 = determinant(la, lb, h1a, h1b);
            FT den2 = determinant(la, lb, h2a, h2b);
            return sign(lb) *
                   sign(num) *
                   sign(den1) *
                   sign(den2);
        }

        public static ST
        compare_x(FT l1a, FT l1b, FT l1c,
            FT h1a, FT h1b, FT h1c,
            FT l2a, FT l2b, FT l2c,
            FT h2a, FT h2b, FT h2c)
        {
            FT num1 = determinant(l1b, l1c, h1b, h1c);
            FT den1 = determinant(l1a, l1b, h1a, h1b);
            FT num2 = determinant(l2b, l2c, h2b, h2c);
            FT den2 = determinant(l2a, l2b, h2a, h2b);
            ST s = sign(den1) * sign(den2);
            CGAL_kernel_assertion(s != 0);
            return s * sign_of_determinant(num1, num2, den1, den2);
        }

        public static ST compare_y_at_x(FT px, FT py,
                 FT la, FT lb, FT lc)
        {
            var s = sign(lb);
            CGAL_kernel_assertion(s != 0);
            return s * sign(la * px + lb * py + lc);
        }

        public static ST compare_y_at_x(FT px,
                 FT l1a, FT l1b, FT l1c,
                 FT l2a, FT l2b, FT l2c)
        {
            ST s = sign(l1b) * sign(l2b);
            CGAL_kernel_assertion(s != 0);
            return s * sign_of_determinant(l2a * px + l2c, l2b, l1a * px + l1c, l1b);
        }

        ////CGAL_KERNEL_LARGE_INLINE
        public static ST compare_y_at_x(
                   FT l1a, FT l1b, FT l1c,
                 FT l2a, FT l2b, FT l2c,
                 FT ha, FT hb, FT hc)
        {
            ST s = sign(hb) * sign_of_determinant(l1a, l1b, l2a, l2b);
            CGAL_kernel_assertion(s != 0);
            return s * sign_of_determinant(l1a, l1b, l1c,
                                              l2a, l2b, l2c,
                                              ha, hb, hc);
        }

        ////CGAL_KERNEL_LARGE_INLINE
        public static ST compare_y_at_x(FT l1a, FT l1b, FT l1c,
                         FT l2a, FT l2b, FT l2c,
                         FT h1a, FT h1b, FT h1c,
                         FT h2a, FT h2b, FT h2c)
        {
            // The abscissa of the intersection point is num/den.
            FT num = determinant(l1b, l1c, l2b, l2c);
            FT den = determinant(l1a, l1b, l2a, l2b);
            ST s = sign(h1b) * sign(h2b) * sign(den);
            CGAL_kernel_assertion(s != 0);
            return s * sign_of_determinant(h2a * num + h2c * den, h2b,
                                                  h1a * num + h1c * den, h1b);
        }

        // forward-declaration of orientationC2, used in compare_y_at_xC2

        public static bool are_ordered(FT a, FT b, FT c)
        {
            FT min = Math.Min(a, c);
            FT max = Math.Max(a, c);
            return min <= b && b <= max;
        }

        // //CGAL_KERNEL_LARGE_INLINE
        public static ST compare_y_at_x(FT px, FT py,
                     FT ssx, FT ssy,
                     FT stx, FT sty)
        {
            // compares the y-coordinates of p and the vertical projection of p on s.
            // Precondition : p is in the x-range of s.

            CGAL_kernel_precondition(are_ordered(ssx, px, stx));

            if (ssx < stx)
                return orientation(px, py, ssx, ssy, stx, sty);
            else if (ssx > stx)
                return orientation(px, py, stx, sty, ssx, ssy);
            else
            {
                if (py < Math.Min(sty, ssy))
                    return -1;
                if (py > Math.Max(sty, ssy))
                    return 1;
                return 0;
            }
        }

        ////CGAL_KERNEL_LARGE_INLINE
        public static ST compare_y_at_x_segment_(FT px,
                          FT s1sx, FT s1sy,
                          FT s1tx, FT s1ty,
                          FT s2sx, FT s2sy,
                          FT s2tx, FT s2ty)
        {
            // compares the y-coordinates of the vertical projections of p on s1 and s2
            // Precondition : p is in the x-range of s1 and s2.
            // - if one or two segments are vertical :
            //   - if the segments intersect, return 0
            //   - if not, return the obvious -1/1.

            CGAL_kernel_precondition(are_ordered(s1sx, px, s1tx));
            CGAL_kernel_precondition(are_ordered(s2sx, px, s2tx));

            if (s1sx != s1tx && s2sx != s2tx)
            {
                FT s1stx = s1sx - s1tx;
                FT s2stx = s2sx - s2tx;

                return compare(s1sx, s1tx) *
                       compare(s2sx, s2tx) *
                       compare(-(s1sx - px) * (s1sy - s1ty) * s2stx,
                                        (s2sy - s1sy) * s2stx * s1stx
                                        - (s2sx - px) * (s2sy - s2ty) * s1stx);
            }
            else
            {
                if (s1sx == s1tx)
                { // s1 is vertical
                    ST c1, c2;
                    c1 = compare_y_at_x(px, s1sy, s2sx, s2sy, s2tx, s2ty);
                    c2 = compare_y_at_x(px, s1ty, s2sx, s2sy, s2tx, s2ty);
                    if (c1 == c2)
                        return c1;
                    return 0;
                }
                // s2 is vertical
                ST c3, c4;
                c3 = compare_y_at_x(px, s2sy, s1sx, s1sy, s1tx, s1ty);
                c4 = compare_y_at_x(px, s2ty, s1sx, s1sy, s1tx, s1ty);
                if (c3 == c4)
                    return -c3;
                return 0;
            }
        }

        public static bool equal_direction(FT dx1, FT dy1,
                          FT dx2, FT dy2)
        {
            return (sign(dx1) == sign(dx2) &&
                               sign(dy1) == sign(dy2) &&
                               sign_of_determinant(dx1, dy1, dx2, dy2) == 0);
        }

        public static ST compare_angle_with_x_axis(FT dx1, FT dy1, FT dx2, FT dy2)
        {
            // angles are in [-pi,pi], and the angle between Ox and d1 is compared
            // with the angle between Ox and d2
            int quadrant_1 = (dx1 >= 0) ? (dy1 >= 0 ? 1 : 4)
                                        : (dy1 >= 0 ? 2 : 3);
            int quadrant_2 = (dx2 >= 0) ? (dy2 >= 0 ? 1 : 4)
                                        : (dy2 >= 0 ? 2 : 3);
            // We can't use compare(quadrant_1,quadrant_2) because in case
            // of tie, we need additional computation
            if (quadrant_1 > quadrant_2)
                return 1;
            else if (quadrant_1 < quadrant_2)
                return -1;
            return -sign_of_determinant(dx1, dy1, dx2, dy2);
        }

        public static ST compare_angle_with_x_axis(Direction2 d1, Direction2 d2) => compare_angle_with_x_axis(d1.DX, d1.DY, d2.DX, d2.DY);

        public static ST compare_slopes(FT l1a, FT l1b, FT l2a, FT l2b)
        {
            if (is_zero(l1a))  // l1 is horizontal
                return is_zero(l2b) ? -1 : sign(l2a) * sign(l2b);
            if (is_zero(l2a)) // l2 is horizontal
                return is_zero(l1b) ? 1 : -sign(l1a) * sign(l1b);
            if (is_zero(l1b)) return is_zero(l2b) ? 0 : 1;
            if (is_zero(l2b)) return -1;
            var l1_sign = -sign(l1a) * sign(l1b);
            var l2_sign = -sign(l2a) * sign(l2b);

            if (l1_sign < l2_sign) return -1;
            if (l1_sign > l2_sign) return 1;

            if (l1_sign > 0)
                return compare(Math.Abs(l1a * l2b), Math.Abs(l2a * l1b));

            return compare(Math.Abs(l2a * l1b),
                                      Math.Abs(l1a * l2b));
        }

        public static ST compare_slopes(FT s1_src_x, FT s1_src_y, FT s1_tgt_x,
                 FT s1_tgt_y, FT s2_src_x, FT s2_src_y,
                 FT s2_tgt_x, FT s2_tgt_y)
        {
            var cmp_y1 = compare(s1_src_y, s1_tgt_y);
            if (cmp_y1 == 0) // horizontal
            {
                if (compare(s2_src_x, s2_tgt_x) == 0) return -1;
                return -sign(s2_src_y - s2_tgt_y) * sign(s2_src_x - s2_tgt_x);
            }

            var cmp_y2 = compare(s2_src_y, s2_tgt_y);
            if (cmp_y2 == 0)
            {
                if (compare(s1_src_x, s1_tgt_x) == 0) return 1;
                return sign(s1_src_y - s1_tgt_y) * sign(s1_src_x - s1_tgt_x);
            }

            var cmp_x1 = compare(s1_src_x, s1_tgt_x);
            var cmp_x2 = compare(s2_src_x, s2_tgt_x);

            if (cmp_x1 == 0) return cmp_x2 == 0 ? 0 : 1;

            if (cmp_x2 == 0) return -1;

            FT s1_xdiff = s1_src_x - s1_tgt_x;
            FT s1_ydiff = s1_src_y - s1_tgt_y;
            FT s2_xdiff = s2_src_x - s2_tgt_x;
            FT s2_ydiff = s2_src_y - s2_tgt_y;
            var s1_sign = sign(s1_ydiff) * sign(s1_xdiff);
            var s2_sign = sign(s2_ydiff) * sign(s2_xdiff);

            if (s1_sign < s2_sign) return -1;
            if (s1_sign > s2_sign) return 1;

            if (s1_sign > 0)
                return compare(Math.Abs(s1_ydiff * s2_xdiff), Math.Abs(s2_ydiff * s1_xdiff));

            return compare(Math.Abs(s2_ydiff * s1_xdiff), Math.Abs(s1_ydiff * s2_xdiff));
        }

        //[MethodImpl(MethodImplOptions.AggressiveInlining)]
        //        public static ST compare_lexicographically_xy(FT px, FT py,
        //                               FT qx, FT qy)
        //    {
        //            ST c = compare(px, qx);
        //        if (is_indeterminate(c)) return indeterminate<Cmp>();
        //        return (c != 0) ? c : compare(py, qy);
        //    }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ST orientation(FT px, FT py, FT qx, FT qy, FT rx, FT ry)
        {
            return sign_of_determinant(qx - px, qy - py, rx - px, ry - py);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool counterclockwise_in_between_2(in Direction2 p, in Direction2 q, in Direction2 r)
        {
            if (q < p)
                return (p < r) || (r <= q);
            else
                return (p < r) && (r <= q);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static OrientationEnum orientation(in Vector2 u, in Vector2 v)
        {
            return (OrientationEnum)sign_of_determinant(u.hx(), u.hy(), v.hx(), v.hy());
        }

        public static OrientationEnum orientation(Point2[] p)
        {
            var phx = p[0].hx();
            var phy = p[0].hy();
            var phw = p[0].hw();
            var qhx = p[1].hx();
            var qhy = p[1].hy();
            var qhw = p[1].hw();
            var rhx = p[2].hx();
            var rhy = p[2].hy();
            var rhw = p[2].hw();

            // | A B |
            // | C D |

            RT A = phx * rhw - phw * rhx;
            RT B = phy * rhw - phw * rhy;
            RT C = qhx * rhw - qhw * rhx;
            RT D = qhy * rhw - qhw * rhy;

            return (OrientationEnum)compare(A * D, B * C);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static OrientationEnum orientation(in Point2 p, in Point2 q, in Point2 r)
        {
            var phx = p.hx();
            var phy = p.hy();
            var phw = p.hw();
            var qhx = q.hx();
            var qhy = q.hy();
            var qhw = q.hw();
            var rhx = r.hx();
            var rhy = r.hy();
            var rhw = r.hw();

            // | A B |
            // | C D |

            RT A = phx * rhw - phw * rhx;
            RT B = phy * rhw - phw * rhy;
            RT C = qhx * rhw - qhw * rhx;
            RT D = qhy * rhw - qhw * rhy;

            return (OrientationEnum)compare(A * D, B * C);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ST orientation(FT ux, FT uy, FT vx, FT vy)
        {
            return sign_of_determinant(ux, uy, vx, vy);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ST angle(FT ux, FT uy, FT vx, FT vy)
        {
            return (sign(ux * vx + uy * vy));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ST angle(
                                        FT px, FT py,
                                        FT qx, FT qy,
                                        FT rx, FT ry)
        {
            return (sign((px - qx) * (rx - qx) + (py - qy) * (ry - qy)));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ST angle(
                                        FT px, FT py,
                                        FT qx, FT qy,
                                        FT rx, FT ry,
                                        FT sx, FT sy)
        {
            return (sign((px - qx) * (rx - sx) + (py - qy) * (ry - sy)));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ST angle(Point2 p, Point2 q, Point2 r, Point2 s) =>

            angle(p.X, p.Y, q.X, q.Y, q.X, r.Y, s.X, s.Y);

        public static bool
        collinear_are_ordered_along_line(FT px, FT py, FT qx, FT qy, FT rx, FT ry)
        {
            if (px < qx) return !(rx < qx);
            if (qx < px) return !(qx < rx);
            if (py < qy) return !(ry < qy);
            if (qy < py) return !(qy < ry);
            return true; // p==q
        }

        public static bool collinear_are_strictly_ordered_along_line(FT px, FT py, FT qx, FT qy, FT rx, FT ry)
        {
            if (px < qx) return (qx < rx);
            if (qx < px) return (rx < qx);
            if (py < qy) return (qy < ry);
            if (qy < py) return (ry < qy);
            return false;
        }

        // //CGAL_KERNEL_LARGE_INLINE
        public static ST side_of_oriented_circle(FT px, FT py,
                                  FT qx, FT qy,
                                  FT rx, FT ry,
                                  FT tx, FT ty)
        {
            //  sign_of_determinant(px, py, px*px + py*py, 1,
            //                         qx, qy, qx*qx + qy*qy, 1,
            //                         rx, ry, rx*rx + ry*ry, 1,
            //                         tx, ty, tx*tx + ty*ty, 1);
            // We first translate so that p is the new origin.
            FT qpx = qx - px;
            FT qpy = qy - py;
            FT rpx = rx - px;
            FT rpy = ry - py;
            FT tpx = tx - px;
            FT tpy = ty - py;
            // The usual 3x3 formula can be simplified a little bit to a 2x2.
            //         - sign_of_determinant(qpx, qpy, square(qpx) + square(qpy),
            //                                  rpx, rpy, square(rpx) + square(rpy),
            //                                  tpx, tpy, square(tpx) + square(tpy)));
            return sign_of_determinant(qpx * tpy - qpy * tpx, tpx * (tx - qx) + tpy * (ty - qy),
                                            qpx * rpy - qpy * rpx, rpx * (rx - qx) + rpy * (ry - qy));
        }

        ////CGAL_KERNEL_LARGE_INLINE

        public static ST? side_of_bounded_circle(FT px, FT py,
                                 FT qx, FT qy,
                                 FT rx, FT ry,
                                 FT tx, FT ty)
        {
            return (side_of_oriented_circle(px, py, qx, qy, rx, ry, tx, ty)
                                           * orientation(px, py, qx, qy, rx, ry));
        }

        public static ST? side_of_bounded_circle(FT px, FT py,
                                 FT qx, FT qy,
                                 FT tx, FT ty)
        {
            // Returns whether T lies inside or outside the circle which diameter is PQ.
            return (compare((tx - px) * (qx - tx), (ty - py) * (ty - qy)));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ST cmp_dist_to_point(FT px, FT py,
                            FT qx, FT qy,
                            FT rx, FT ry)
        {
            return compare(squared_distance(px, py, qx, qy),
                                    squared_distance(px, py, rx, ry));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ST compare_distance_2(Point2 p, Point2 q, Point2 r)
        {
            return compare(squared_distance(p, q),
                                    squared_distance(p, r));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool has_larger_dist_to_point(FT px, FT py, FT qx, FT qy, FT rx, FT ry)
        {
            return cmp_dist_to_point(px, py, qx, qy, rx, ry) == 1;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool has_smaller_dist_to_point(FT px, FT py, FT qx, FT qy, FT rx, FT ry)
        {
            return cmp_dist_to_point(px, py, qx, qy, rx, ry) == -1;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ST cmp_signed_dist_to_direction(FT la, FT lb,
                                       FT px, FT py,
                                       FT qx, FT qy)
        {
            return compare(scaled_distance_to_direction(la, lb, px, py), scaled_distance_to_direction(la, lb, qx, qy));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool has_larger_signed_dist_to_direction(FT la, FT lb,
                                              FT px, FT py,
                                              FT qx, FT qy)
        {
            return cmp_signed_dist_to_direction(la, lb, px, py, qx, qy) == 1;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool has_smaller_signed_dist_to_direction(FT la, FT lb,
                                               FT px, FT py,
                                               FT qx, FT qy)
        {
            return cmp_signed_dist_to_direction(la, lb, px, py, qx, qy) == -1;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ST cmp_signed_dist_to_line(FT px, FT py,
                                  FT qx, FT qy,
                                  FT rx, FT ry,
                                  FT sx, FT sy)
        {
            return compare(scaled_distance_to_line(px, py, qx, qy, rx, ry),
                                    scaled_distance_to_line(px, py, qx, qy, sx, sy));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool has_larger_signed_dist_to_line(FT px, FT py,
                                         FT qx, FT qy,
                                         FT rx, FT ry,
                                         FT sx, FT sy)
        {
            return cmp_signed_dist_to_line(px, py, qx, qy, rx, ry, sx, sy) == 1;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool has_smaller_signed_dist_to_line(FT px, FT py,
                                          FT qx, FT qy,
                                          FT rx, FT ry,
                                          FT sx, FT sy)
        {
            return cmp_signed_dist_to_line(px, py, qx, qy, rx, ry, sx, sy) == -1;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ST side_of_oriented_line(FT a, FT b, FT c,
                                FT x, FT y)
        {
            return sign(a * x + b * y + c);
        }

        public static ST compare_power_distance(FT px, FT py, FT pwt,
                                 FT qx, FT qy, FT qwt,
                                 FT rx, FT ry)
        {
            // returns -1 if r is closer to p w.r.t. the power metric
            FT d1 = square(rx - px) + square(ry - py) - pwt;
            FT d2 = square(rx - qx) + square(ry - qy) - qwt;
            return compare(d1, d2);
        }

        public static ST power_side_of_bounded_power_circle(
                                             FT px, FT py, FT pw,
                                             FT qx, FT qy, FT qw,
                                             FT tx, FT ty, FT tw)
        {
            FT dpx = px - qx;
            FT dpy = py - qy;
            FT dtx = tx - qx;
            FT dty = ty - qy;
            FT dpz = square(dpx) + square(dpy);

            return
              (sign(-(square(dtx) + square(dty) - tw + qw) * dpz
                             + (dpz - pw + qw) * (dpx * dtx + dpy * dty)));
        }

        public static OrientedSideEnum power_side_of_oriented_power_circle(FT px, FT py, FT pwt,
                                              FT qx, FT qy, FT qwt,
                                              FT rx, FT ry, FT rwt,
                                              FT tx, FT ty, FT twt)
        {
            // Note: maybe this can be further optimized like the usual in_circle() ?

            // We translate the 4 points so that T becomes the origin.
            FT dpx = px - tx;
            FT dpy = py - ty;
            FT dpz = square(dpx) + square(dpy) - pwt + twt;
            FT dqx = qx - tx;
            FT dqy = qy - ty;
            FT dqz = square(dqx) + square(dqy) - qwt + twt;
            FT drx = rx - tx;
            FT dry = ry - ty;
            FT drz = square(drx) + square(dry) - rwt + twt;

            return (OrientedSideEnum)sign_of_determinant(dpx, dpy, dpz,
                                       dqx, dqy, dqz,
                                       drx, dry, drz);
        }

        public static OrientedSideEnum power_side_of_oriented_power_circle(FT px, FT py, FT pwt,
                                          FT qx, FT qy, FT qwt,
                                          FT tx, FT ty, FT twt)
        {
            // Same translation as above.
            FT dpx = px - tx;
            FT dpy = py - ty;
            FT dpz = square(dpx) + square(dpy) - pwt + twt;
            FT dqx = qx - tx;
            FT dqy = qy - ty;
            FT dqz = square(dqx) + square(dqy) - qwt + twt;

            // We do an orthogonal projection on the (x) axis, if possible.
            ST cmpx = compare(px, qx);
            if (cmpx != 0)
                return (OrientedSideEnum)(cmpx * sign_of_determinant(dpx, dpz, dqx, dqz));

            // If not possible, then on the (y) axis.
            var cmpy = compare(py, qy);
            return (OrientedSideEnum)(cmpy * sign_of_determinant(dpy, dpz, dqy, dqz));
        }

        public static OrientedSideEnum circumcenter_oriented_side_of_oriented_segment(FT ax, FT ay,
                                                         FT bx, FT by,
                                                         FT p0x, FT p0y,
                                                         FT p1x, FT p1y,
                                                         FT p2x, FT p2y)
        {
            FT dX = bx - ax;
            FT dY = by - ay;
            FT R0 = p0x * p0x + p0y * p0y;
            FT R1 = p1x * p1x + p1y * p1y;
            FT R2 = p2x * p2x + p2y * p2y;
            FT denominator = (p1x - p0x) * (p2y - p0y) +
                                   (p0x - p2x) * (p1y - p0y);
            FT det = 2 * denominator * (ax * dY - ay * dX)
                             - (R2 - R1) * (p0x * dX + p0y * dY)
                             - (R0 - R2) * (p1x * dX + p1y * dY)
                             - (R1 - R0) * (p2x * dX + p2y * dY);
            return (OrientedSideEnum)sign(det);
        }
    }

    public static partial class Mathex
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void midpointC2(FT px, FT py,
                    FT qx, FT qy,
                   out FT x, out FT y)
        {
            x = (px + qx) / 2;
            y = (py + qy) / 2;
        }

        //CGAL_KERNEL_LARGE_INLINE
        public static void circumcenter_translateC2(FT dqx, FT dqy, FT drx, FT dry, out FT dcx, out FT dcy)
        {
            // Given 3 points P, Q, R, this function takes as input:
            // qx-px, qy-py, rx-px, ry-py.  And returns cx-px, cy-py,
            // where (cx, cy) are the coordinates of the circumcenter C.

            // What we do is intersect the bisectors.
            FT r2 = square(drx) + square(dry);
            FT q2 = square(dqx) + square(dqy);
            FT den = 2 * determinant(dqx, dqy, drx, dry);

            // The 3 points aren't collinear.
            // Hopefully, this is already checked at the upper level.
            CGAL_kernel_assertion(!is_zero(den));

            // One possible optimization here is to precompute 1/den, to avoid one
            // division.  However, we lose precision, and it's maybe not worth it (?).
            dcx = determinant(dry, dqy, r2, q2) / den;
            dcy = -determinant(drx, dqx, r2, q2) / den;
        }

        // CGAL_KERNEL_MEDIUM_INLINE
        public static void circumcenterC2(FT px, FT py, FT qx, FT qy, FT rx, FT ry, out FT x, out FT y)
        {
            circumcenter_translateC2(qx - px, qy - py, rx - px, ry - py, out x, out y);
            x += px;
            y += py;
        }

        public static void barycenterC2(FT p1x, FT p1y, FT w1,
                    FT p2x, FT p2y,
                    out FT x, out FT y)
        {
            FT w2 = 1 - w1;
            x = w1 * p1x + w2 * p2x;
            y = w1 * p1y + w2 * p2y;
        }

        public static void barycenterC2(FT p1x, FT p1y, FT w1,
                    FT p2x, FT p2y, FT w2,
                    out FT x, out FT y)
        {
            FT sum = w1 + w2;
            CGAL_kernel_assertion(sum != 0);
            x = (w1 * p1x + w2 * p2x) / sum;
            y = (w1 * p1y + w2 * p2y) / sum;
        }

        public static void barycenterC2(FT p1x, FT p1y, FT w1,
                    FT p2x, FT p2y, FT w2,
                    FT p3x, FT p3y,
                    out FT x, out FT y)
        {
            FT w3 = 1 - w1 - w2;
            x = w1 * p1x + w2 * p2x + w3 * p3x;
            y = w1 * p1y + w2 * p2y + w3 * p3y;
        }

        public static void barycenterC2(FT p1x, FT p1y, FT w1,
                    FT p2x, FT p2y, FT w2,
                    FT p3x, FT p3y, FT w3,
                    out FT x, out FT y)
        {
            FT sum = w1 + w2 + w3;
            CGAL_kernel_assertion(sum != 0);
            x = (w1 * p1x + w2 * p2x + w3 * p3x) / sum;
            y = (w1 * p1y + w2 * p2y + w3 * p3y) / sum;
        }

        public static void barycenterC2(FT p1x, FT p1y, FT w1,
                    FT p2x, FT p2y, FT w2,
                    FT p3x, FT p3y, FT w3,
                    FT p4x, FT p4y,
                    out FT x, out FT y)
        {
            FT w4 = 1 - w1 - w2 - w3;
            x = w1 * p1x + w2 * p2x + w3 * p3x + w4 * p4x;
            y = w1 * p1y + w2 * p2y + w3 * p3y + w4 * p4y;
        }

        public static void barycenterC2(FT p1x, FT p1y, FT w1,
                    FT p2x, FT p2y, FT w2,
                    FT p3x, FT p3y, FT w3,
                    FT p4x, FT p4y, FT w4,
                    out FT x, out FT y)
        {
            FT sum = w1 + w2 + w3 + w4;
            CGAL_kernel_assertion(sum != 0);
            x = (w1 * p1x + w2 * p2x + w3 * p3x + w4 * p4x) / sum;
            y = (w1 * p1y + w2 * p2y + w3 * p3y + w4 * p4y) / sum;
        }

        // CGAL_KERNEL_MEDIUM_INLINE
        public static void centroidC2(FT px, FT py,
                   FT qx, FT qy,
                   FT rx, FT ry,
                   out FT x, out FT y)
        {
            x = (px + qx + rx) / 3;
            y = (py + qy + ry) / 3;
        }

        // CGAL_KERNEL_MEDIUM_INLINE
        public static void centroidC2(FT px, FT py,
                   FT qx, FT qy,
                   FT rx, FT ry,
                   FT sx, FT sy,
                   out FT x, out FT y)
        {
            x = (px + qx + rx + sx) / 4;
            y = (py + qy + ry + sy) / 4;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void
        line_from_pointsC2(FT px, FT py,
                          FT qx, FT qy,
                          out FT a, out FT b, out FT c)
        {
            // The horizontal and vertical line get a special treatment
            // in order to make the intersection code robust for doubles
            if (py == qy)
            {
                a = 0;
                if (qx > px)
                {
                    b = 1;
                    c = -py;
                }
                else if (qx == px)
                {
                    b = 0;
                    c = 0;
                }
                else
                {
                    b = -1;
                    c = py;
                }
            }
            else if (qx == px)
            {
                b = 0;
                if (qy > py)
                {
                    a = -1;
                    c = px;
                }
                else if (qy == py)
                {
                    a = 0;
                    c = 0;
                }
                else
                {
                    a = 1;
                    c = -px;
                }
            }
            else
            {
                a = py - qy;
                b = qx - px;
                c = -px * a - py * b;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void
        line_from_point_directionC2(FT px, FT py,
                                   FT dx, FT dy,
                                   out FT a, out FT b, out FT c)
        {
            a = -dy;
            b = dx;
            c = px * dy - py * dx;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT approximate_sqrt(FT v) => Math.Sqrt(v);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT sqrt(FT v) => Math.Sqrt(v);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void
        bisector_of_pointsC2(FT px, FT py, FT dx, FT dy, FT qx, FT qy, out FT a, out FT b, out FT c)
        {
            a = 2 * (px - qx);
            b = 2 * (py - qy);
            c = square(qx) + square(qy) -
                square(px) - square(py);
        }

        //[MethodImpl(MethodImplOptions.AggressiveInlining)]
        //public static void
        //bisector_of_pointsC2(FT px, FT py, FT dx, FT dy, FT qx, FT qy, out FT a, out FT b, out FT c)
        //{
        //    a = 2 * (px - qx);
        //    b = 2 * (py - qy);
        //    c = square(qx) + square(qy) -
        //        square(px) - square(py);
        //}

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void bisector_of_linesC2(FT pa, FT pb, FT pc,
                          FT qa, FT qb, FT qc,
                          out FT a, out FT b, out FT c)
        {
            // We normalize the equations of the 2 lines, and we then add them.
            FT n1 = approximate_sqrt((square(pa) + square(pb)));
            FT n2 = approximate_sqrt((square(qa) + square(qb)));
            a = n2 * pa + n1 * qa;
            b = n2 * pb + n1 * qb;
            c = n2 * pc + n1 * qc;

            // Care must be taken for the case when this produces a degenerate line.
            if (a == 0 && b == 0)
            {
                a = n2 * pa - n1 * qa;
                b = n2 * pb - n1 * qb;
                c = n2 * pc - n1 * qc;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT line_y_at_xC2(FT a, FT b, FT c, FT x)
        {
            return (-a * x - c) / b;
        }

        // Silence a warning for MSVC 2017
        // > include\cgal\constructions\kernel_ftc2.h(287) :
        // >   warning C4723: potential divide by 0

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void line_get_pointC2(FT a, FT b, FT c, FT i, out FT x, out FT y)
        {
            if (is_zero(b))
            {
                x = -c / a;
                y = 1 - i * a;
            }
            else
            {
                x = 1 + i * b;
                y = -(a + c) / b - i * a;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void perpendicular_through_pointC2(FT la, FT lb,
                                    FT px, FT py,
                                    out FT a, out FT b, out FT c)
        {
            a = -lb;
            b = la;
            c = lb * px - la * py;
        }

        // CGAL_KERNEL_MEDIUM_INLINE
        public static void line_project_pointC2(FT la, FT lb, FT lc, FT px, FT py, out FT x, out FT y)
        {
            if ((is_zero(la))) // horizontal line
            {
                x = px;
                y = -lc / lb;
            }
            else if ((is_zero(lb))) // vertical line
            {
                x = -lc / la;
                y = py;
            }
            else
            {
                FT a2 = square(la);
                FT b2 = square(lb);
                FT d = a2 + b2;
                x = (b2 * px - la * lb * py - la * lc) / d;
                y = (-la * lb * px + a2 * py - lb * lc) / d;
            }
        }

        // CGAL_KERNEL_MEDIUM_INLINE
        public static FT squared_radiusC2(FT px, FT py, FT qx, FT qy, FT rx, FT ry, out FT x, out FT y)
        {
            circumcenter_translateC2(qx - px, qy - py, rx - px, ry - py, out x, out y);
            FT r2 = square(x) + square(y);
            x += px;
            y += py;
            return r2;
        }

        // CGAL_KERNEL_MEDIUM_INLINE
        public static FT squared_radiusC2(FT px, FT py,
               FT qx, FT qy,
               FT rx, FT ry)
        {
            FT x, y;
            circumcenter_translateC2(qx - px, qy - py, rx - px, ry - py, out x, out y);
            return square(x) + square(y);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT squared_distanceC2(FT px, FT py, FT qx, FT qy)
        {
            return square(px - qx) + square(py - qy);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT squared_radiusC2(FT px, FT py,
                       FT qx, FT qy)
        {
            return squared_distanceC2(px, py, qx, qy) / 4f;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT scaled_distance_to_lineC2(FT la, FT lb, FT lc,
                                 FT px, FT py)
        {
            // for comparisons, use distance_to_directionsC2 instead
            // since lc is irrelevant
            return la * px + lb * py + lc;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static FT scaled_distance_to_directionC2(FT la, FT lb, FT px, FT py)
        {
            // scalar product with direction
            return la * px + lb * py;
        }

        // CGAL_KERNEL_MEDIUM_INLINE
        public static FT scaled_distance_to_lineC2(FT px, FT py, FT qx, FT qy, FT rx, FT ry)
        {
            return determinant(px - rx, py - ry, qx - rx, qy - ry);
        }

        public static void weighted_circumcenter_translateC2(RT dqx, RT dqy, RT dqw,
                                        RT drx, RT dry, RT drw,
                                       out RT dcx, out RT dcy)
        {
            // Given 3 points P, Q, R, this function takes as input:
            // qx-px, qy-py,qw-pw,  rx-px, ry-py, rw-pw.  And returns cx-px, cy-py,
            // where (cx, cy) are the coordinates of the circumcenter C.

            // What we do is intersect the radical axis
            RT r2 = square(drx) + square(dry) - drw;
            RT q2 = square(dqx) + square(dqy) - dqw;

            RT den = 2 * determinant(dqx, dqy, drx, dry);

            // The 3 points aren't collinear.
            // Hopefully, this is already checked at the upper level.
            CGAL_assertion(den != 0);

            // One possible optimization here is to precompute 1/den, to avoid one
            // division.  However, we lose precision, and it's maybe not worth it (?).
            dcx = determinant(dry, dqy, r2, q2) / den;
            dcy = -determinant(drx, dqx, r2, q2) / den;
        }

        //template < class RT >

        public static void weighted_circumcenterC2(RT px, RT py, We pw,
                               RT qx, RT qy, We qw,
                               RT rx, RT ry, We rw,
                               out RT x, out RT y)
        {
            RT dqw = (RT)(qw - pw);
            RT drw = (RT)(rw - pw);

            weighted_circumcenter_translateC2(qx - px, qy - py, dqw, rx - px, ry - py, drw, out x, out y);
            x += px;
            y += py;
        }

        public static FT power_productC2(FT px, FT py, FT pw, FT qx, FT qy, FT qw)
        {
            // computes the power product of two weighted points
            FT qpx = qx - px;
            FT qpy = qy - py;
            FT qp2 = square(qpx) + square(qpy);
            return qp2 - pw - qw;
        }

        public static void radical_axisC2(RT px, RT py, We pw,
                     RT qx, RT qy, We qw,
                     out RT a, out RT b, out RT c)
        {
            a = (RT)(2) * (px - qx);
            b = (RT)(2) * (py - qy);
            c = -square(px) - square(py)
                + square(qx) + square(qy)
                + (RT)(pw) - (RT)(qw);
        }

        // CGAL_KERNEL_MEDIUM_INLINE
        public static FT
        squared_radius_orthogonal_circleC2(FT px, FT py, FT pw, FT qx, FT qy, FT qw, FT rx, FT ry, FT rw)
        {
            FT FT4 = (4);
            FT dpx = px - rx;
            FT dpy = py - ry;
            FT dqx = qx - rx;
            FT dqy = qy - ry;
            FT dpp = square(dpx) + square(dpy) - pw + rw;
            FT dqq = square(dqx) + square(dqy) - qw + rw;

            FT det0 = determinant(dpx, dpy, dqx, dqy);
            FT det1 = determinant(dpp, dpy, dqq, dqy);
            FT det2 = determinant(dpx, dpp, dqx, dqq);

            return (square(det1) + square(det2)) /
                                            (FT4 * square(det0)) - rw;
        }

        // CGAL_KERNEL_MEDIUM_INLINE
        public static FT squared_radius_smallest_orthogonal_circleC2(FT px, FT py, FT pw, FT qx, FT qy, FT qw)
        {
            FT FT4 = (4);
            FT dpz = square(px - qx) + square(py - qy);
            return (square(dpz - pw + qw) / (FT4 * dpz) - qw);
        }

        public static bool collinear_are_ordered_along_line(in Point2 p, in Point2 q, in Point2 r)
        {
            return collinear_are_ordered_along_line(p.X, p.Y, q.X, q.Y, r.X, r.Y);
        }

        public static bool collinear(in Point2 p, in Point2 q, in Point2 r)
        {
            return (orientation(p, q, r) == OrientationEnum.COLLINEAR);
        }

        public static bool are_ordered_along_line(in Point2 p, in Point2 q, in Point2 r) => collinear(p, q, r) && collinear_are_ordered_along_line(p, q, r);
    }
}