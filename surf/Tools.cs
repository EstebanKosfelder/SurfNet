using System.Text.RegularExpressions;
using System;
using System.Numerics;
using System.Security.Cryptography.X509Certificates;
using System.Diagnostics;
using TriangleNet.Topology.DCEL;

namespace SurfNet
{
    public static class Tools
    {
//        public static readonly Vector2 ORIGIN = new Vector2(0, 0);

//        public const double CORE_ZERO = 0.0;

//        public const int SMALLER = -1;
//        public const int EQUAL = 0;
//        public const int LARGER = 1;

//        public const int NEGATIVE = -1;
//        public const int ZERO = 0;
//        public const int POSITIVE = 1;

//        public const int LEFT_TURN = -1;
//        public const int COLLINEAR = 0;
//        public const int RIGHT_TURN = 1;

//        public const int CONVEX = LEFT_TURN;
//        public const int STRAIGHT = COLLINEAR;
//        public const int REFLEX = RIGHT_TURN;

//        public const int COUNTERCLOCKWISE = 1;
//        public const int CLOCKWISE = -1;

//        public const int ON_NEGATIVE_SIDE = -1;
//        public const int ON_ORIENTED_BOUNDARY = 0;
//        public const int ON_POSITIVE_SIDE = 1;

//        public const double EPS = 1e-9;
//        public const double FASTER_EDGE_WINS_IN_COLLINEAR_CASES = 1.0;



      
//        public static Polynomial1D differentiate(Polynomial1D p) => p.differentiate();

//        public static double evaluate(Polynomial1D p, double at) => p.evaluate(at);

//        public static (bool has_real_roots, bool is_square) solve_quadratic(Polynomial1D f, out double x0, out double x1) => f.solve_quadratic(out x0, out x1);

//        public static Polynomial1D compute_determinant(Polynomial1D x0, Polynomial1D y0,
//                             Polynomial1D x1, Polynomial1D y1,
//                             Polynomial1D x2, Polynomial1D y2)
//        {
//            return (
//              (x0 * y1
//              + x1 * y2
//              + x2 * y0
//              )
//              -
//              (y0 * x1
//              + y1 * x2
//              + y2 * x0
//              )
//            );
//        }

//        public static void swap<T>(ref T a, ref T b)
//        { (b, a) = (a, b); }

//        public static (int, int, int) indirect_sort_3<T>(T[] t) where T : IComparable
//        {
//            int i0 = 0;
//            int i1 = 1;
//            int i2 = 2;
//            if (t[i0].CompareTo(t[i1]) > 0) swap(ref i0, ref i1);
//            if (t[i1].CompareTo(t[i2]) > 0) swap(ref i1, ref i2);
//            if (t[i0].CompareTo(t[i1]) > 0) swap(ref i0, ref i1);
//            return (i0, i1, i2);
//        }

//        public static double dot(Vector2 a, Vector2 b) => a.X * b.X + a.Y * b.Y;

//        public static double squared_distance(Point2 a, Point2 b) => (a - b).squared_distance();

//        public static int sign(double value) => Math.Abs(value) < EPS ? 0 : Math.Sign(value);


//        public static double sqrt(double x) => Math.Sqrt(x);

//        public static int orientation(in Vector2 p, in Vector2 q) => orientation(p.X, p.Y, q.X, q.Y);

//        public static int orientation(in Point2 p, in Point2 q) => orientation(p.X, p.Y, q.X, q.Y);

//        public static int orientation(in Point2 p, in Point2 q, in Point2 r) => orientation(p.X, p.Y, q.X, q.Y, r.X, r.Y);

//        internal static int CompareToEps(this double a, double b, double eps = EPS)
//        {
//            return a.AreNear(b, EPS) ? 0 : a.CompareTo(b);
//        }

        
//        internal static bool AreNear(this Point2 a, Point2 b, double eps = EPS)
//        {
//            return Mathex.AreNear(a.X, b.X, eps) && Mathex.AreNear(a.Y, b.Y, eps);
//        }

//        // Not documented
//        public static bool
//                are_sorted(double a, double b, double c)
//        {
//            return a <= b && b <= c;
//        }

//        /*
//        // Not documented
//        public static
//            bool
//                are_sorted(double a, double b, double c, Compare cmp)
//        {
//            return !cmp(b, a) && !cmp(c, b);
//        }
//*/

//        // Not documented
//        public static bool
//                are_strictly_sorted(double a, double b, double c)
//        {
//            return a < b && b < c;
//        }

//        /*
//                // Not documented
//                public static
//                    bool
//                        are_strictly_sorted(double a, double b, double c, Compare cmp)
//                {
//                    return cmp(a, b) && cmp(b, c);
//                }
//        */

//        // Not documented
//        // Checks that b is in the interval [min(a, c) , max(a, c)].
//        public static bool
//                are_ordered(double a, double b, double c)
//        {
//            double min = Math.Min(a, c);
//            double max = Math.Max(a, c);
//            return min <= b && b <= max;
//        }

//        /*
//        // Not documented
//        // Checks that b is in the interval [min(a, c) , max(a, c)].
//        public static
//            bool
//                are_ordered(double a, double b, double c, Compare cmp)
//        {
//            double min = Math.Min(a, c, cmp);
//            double max = Math.Max(a, c, cmp);
//            return !cmp(b, min) && !cmp(max, b);
//        }
//*/

//        // Not documented
//        // Checks that b is in the interval ]min(a, c) , max(a, c)[.
//        public static bool
//                are_strictly_ordered(double a, double b, double c)
//        {
//            double min = Math.Min(a, c);
//            double max = Math.Max(a, c);
//            return min < b && b < max;
//        }

//        /*
//        // Not documented
//        // Checks that b is in the interval ]min(a, c) , max(a, c)[.
//        public static
//            bool
//                are_strictly_ordered(double a, double b, double c, Compare cmp)
//        {
//            double min = (std.min)(a, c, cmp);
//            double max = (std.max)(a, c, cmp);
//            return cmp(min, b) && cmp(b, max);
//        }
//*/

//        public static double ToDegree(double a)
//        {
//            return a * 180 / Math.PI;
//        }

//        public static void Swap(ref double a, ref double b)
//        {
//            double c = a;
//            a = b;
//            b = c;
//        }

//        public static void Swap<T>(ref T a, ref T b)
//        {
//            T c = a;
//            a = b;
//            b = c;
//        }

//        public static double square(double a)
//        {
//            return a * a;
//        }

//        public static void
//                midpointC2(double px, double py,
//                     double qx, double qy,
//                     out double x, out double y)
//        {
//            x = (px + qx) / 2;
//            y = (py + qy) / 2;
//        }

//        public static void
//                circumcedoubleer_translateC2(double dqx, double dqy,
//                                       double drx, double dry,
//                                       out double dcx, out double dcy)
//        {
//            // Given 3 points P, Q, R, this function takes as input:
//            // qx-px, qy-py, rx-px, ry-py.  And returns cx-px, cy-py,
//            // where (cx, cy) are the coordinates of the circumcedoubleer C.

//            // What we do is idoubleersect the bisectors.
//            double r2 = square(drx) + square(dry);
//            double q2 = square(dqx) + square(dqy);
//            double den = 2 * determinant(dqx, dqy, drx, dry);

//            // The 3 points aren't collinear.
//            // Hopefully, this is already checked at the upper level.
//            if (isZero(den))
//                throw new Exception();

//            // One possible optimization here is to precompute 1/den, to avoid one
//            // division.  However, we loose precision, and it's maybe not worth it (?).
//            dcx = determinant(dry, dqy, r2, q2) / den;
//            dcy = -determinant(drx, dqx, r2, q2) / den;
//        }

//        public static void
//                circumcedoubleerC2(double px, double py,
//                             double qx, double qy,
//                             double rx, double ry,
//                             out double x, out double y)
//        {
//            circumcedoubleer_translateC2(qx - px, qy - py, rx - px, ry - py, out x, out y);
//            x += px;
//            y += py;
//        }

//        public static void
//                barycedoubleerC2(double p1x, double p1y, double w1,
//                           double p2x, double p2y,
//                           out double x, out double y)
//        {
//            double w2 = 1 - w1;
//            x = w1 * p1x + w2 * p2x;
//            y = w1 * p1y + w2 * p2y;
//        }

//        public static void
//                barycedoubleerC2(double p1x, double p1y, double w1,
//                           double p2x, double p2y, double w2,
//                           out double x, out double y)
//        {
//            double sum = w1 + w2;
//            if (isZero(sum))
//                throw new Exception();
//            x = (w1 * p1x + w2 * p2x) / sum;
//            y = (w1 * p1y + w2 * p2y) / sum;
//        }

//        public static void
//                barycedoubleerC2(double p1x, double p1y, double w1,
//                           double p2x, double p2y, double w2,
//                           double p3x, double p3y,
//                           out double x, out double y)
//        {
//            double w3 = 1 - w1 - w2;
//            x = w1 * p1x + w2 * p2x + w3 * p3x;
//            y = w1 * p1y + w2 * p2y + w3 * p3y;
//        }

//        public static void
//                barycedoubleerC2(double p1x, double p1y, double w1,
//                           double p2x, double p2y, double w2,
//                           double p3x, double p3y, double w3,
//                           out double x, out double y)
//        {
//            double sum = w1 + w2 + w3;
//            if (isZero(sum))
//                throw new Exception();
//            x = (w1 * p1x + w2 * p2x + w3 * p3x) / sum;
//            y = (w1 * p1y + w2 * p2y + w3 * p3y) / sum;
//        }

//        public static void
//                barycedoubleerC2(double p1x, double p1y, double w1,
//                           double p2x, double p2y, double w2,
//                           double p3x, double p3y, double w3,
//                           double p4x, double p4y,
//                           out double x, out double y)
//        {
//            double w4 = 1 - w1 - w2 - w3;
//            x = w1 * p1x + w2 * p2x + w3 * p3x + w4 * p4x;
//            y = w1 * p1y + w2 * p2y + w3 * p3y + w4 * p4y;
//        }

//        public static void
//                barycedoubleerC2(double p1x, double p1y, double w1,
//                           double p2x, double p2y, double w2,
//                           double p3x, double p3y, double w3,
//                           double p4x, double p4y, double w4,
//                           out double x, out double y)
//        {
//            double sum = w1 + w2 + w3 + w4;
//            if (isZero(sum))
//                throw new Exception();
//            x = (w1 * p1x + w2 * p2x + w3 * p3x + w4 * p4x) / sum;
//            y = (w1 * p1y + w2 * p2y + w3 * p3y + w4 * p4y) / sum;
//        }

//        public static void
//                cedoubleroidC2(double px, double py,
//                         double qx, double qy,
//                         double rx, double ry,
//                         out double x, out double y)
//        {
//            x = (px + qx + rx) / 3;
//            y = (py + qy + ry) / 3;
//        }

//        public static void
//                cedoubleroidC2(double px, double py,
//                         double qx, double qy,
//                         double rx, double ry,
//                         double sx, double sy,
//                         out double x, out double y)
//        {
//            x = (px + qx + rx + sx) / 4;
//            y = (py + qy + ry + sy) / 4;
//        }

//        public static void
//                line_from_pointsC2(double px, double py,
//                             double qx, double qy,
//                             out double a, out double b, out double c)
//        {
//            // The horizodoubleal and vertical line get a special treatmedouble
//            // in order to make the idoubleersection code robust for doubles
//            if (py == qy)
//            {
//                a = 0;
//                if (qx > px)
//                {
//                    b = 1;
//                    c = -py;
//                }
//                else if (qx == px)
//                {
//                    b = 0;
//                    c = 0;
//                }
//                else
//                {
//                    b = -1;
//                    c = py;
//                }
//            }
//            else if (qx == px)
//            {
//                b = 0;
//                if (qy > py)
//                {
//                    a = -1;
//                    c = px;
//                }
//                else if (qy == py)
//                {
//                    a = 0;
//                    c = 0;
//                }
//                else
//                {
//                    a = 1;
//                    c = -px;
//                }
//            }
//            else
//            {
//                a = py - qy;
//                b = qx - px;
//                c = -px * a - py * b;
//            }
//        }

//        public static void
//                line_from_point_directionC2(double px, double py,
//                                      double dx, double dy,
//                                      out double a, out double b, out double c)
//        {
//            a = -dy;
//            b = dx;
//            c = px * dy - py * dx;
//        }

//        public static void
//                bisector_of_pointsC2(double px, double py,
//                               double qx, double qy,
//                               out double a, out double b, out double c)
//        {
//            a = 2 * (px - qx);
//            b = 2 * (py - qy);
//            c = square(qx) + square(qy) -
//            square(px) - square(py);
//        }

//        public static void
//                bisector_of_linesC2(double pa, double pb, double pc,
//                              double qa, double qb, double qc,
//                              out double a, out double b, out double c)
//        {
//            // We normalize the equations of the 2 lines, and we then add them.
//            double n1 = Math.Sqrt(square(pa) + square(pb));
//            double n2 = Math.Sqrt(square(qa) + square(qb));
//            a = n2 * pa + n1 * qa;
//            b = n2 * pb + n1 * qb;
//            c = n2 * pc + n1 * qc;

//            // Care must be taken for the case when this produces a degenerate line.
//            if (a == 0 && b == 0)
//            {
//                a = n2 * pa - n1 * qa;
//                b = n2 * pb - n1 * qb;
//                c = n2 * pc - n1 * qc;
//            }
//        }

//        public static double
//                line_y_at_xC2(double a, double b, double c, double x)
//        {
//            return (-a * x - c) / b;
//        }

//        public static void
//                line_get_pointC2(double a, double b, double c, double i,
//                           out double x, out double y)
//        {
//            if (isZero(b))
//            {
//                x = (-b - c) / a + i * b;
//                y = 1 - i * a;
//            }
//            else
//            {
//                x = 1 + i * b;
//                y = -(a + c) / b - i * a;
//            }
//        }

//        public static void
//                perpendicular_through_pointC2(double la, double lb,
//                                        double px, double py,
//                                        out double a, out double b, out double c)
//        {
//            a = -lb;
//            b = la;
//            c = lb * px - la * py;
//        }

//        public static Vector2 line_project_point(Line2 l, Vector2 p)
//        {
//            double x, y;
//            // New version, with more multiplications, but less divisions and tests.
//            // Let's compare the results of the 2, benchmark them, as well as check
//            // the precision with the idoubleervals.
//            double a2 = l.A * l.A;
//            double b2 = l.B * l.B;
//            double d = a2 + b2;
//            x = (l.A * (l.B * p.Y - l.C) - p.X * b2) / d;
//            y = (l.B * (l.C - l.A * p.X) + p.Y * a2) / d;

//            return new Vector2(x, y);
//        }

//        public static double
//                squared_radiusC2(double px, double py,
//                           double qx, double qy,
//                           double rx, double ry,
//                           out double x, out double y)
//        {
//            circumcedoubleer_translateC2(qx - px, qy - py, rx - px, ry - py, out x, out y);
//            double r2 = square(x) + square(y);
//            x += px;
//            y += py;
//            return r2;
//        }

//        public static double squared_radiusC2(double px, double py,
//                                  double qx, double qy,
//                                  double rx, double ry)
//        {
//            double x, y;
//            circumcedoubleer_translateC2(qx - px, qy - py, rx - px, ry - py, out x, out y);
//            return square(x) + square(y);
//        }

//        public static double
//                squared_distanceC2(double px, double py,
//                             double qx, double qy)
//        {
//            return square(px - qx) + square(py - qy);
//        }

//        public static double
//                squared_radiusC2(double px, double py,
//                           double qx, double qy)
//        {
//            return squared_distanceC2(px, py, qx, qy) / 4;
//        }

//        public static double
//                scaled_distance_to_lineC2(double la, double lb, double lc,
//                                    double px, double py)
//        {
//            // for comparisons, use distance_to_directionsC2 instead
//            // since lc is irrelevadouble
//            return la * px + lb * py + lc;
//        }

//        public static double
//                scaled_distance_to_directionC2(double la, double lb,
//                                         double px, double py)
//        {
//            // scalar product with direction
//            return la * px + lb * py;
//        }

//        public static double
//                scaled_distance_to_lineC2(double px, double py,
//                                    double qx, double qy,
//                                    double rx, double ry)
//        {
//            return determinant(px - rx, py - ry, qx - rx, qy - ry);
//        }

//        public static void
//                rational_rotation_approximation(double dirx,     // dir.x()
//                                          double diry,     // dir.y()
//                                          ref double sin_num,  // return
//                                          ref double cos_num,  // return
//                                          ref double denom,    // return
//                                          double eps_num,  // quality_bound
//                                          double eps_den)
//        {
//            double n = eps_num;
//            double d = eps_den;

//            if (isZero(dirx) && isZero(diry))
//                throw new Exception();
//            if (n <= 0.0d)
//                throw new Exception();
//            if (d <= 0.0d)
//                throw new Exception();

//            double sin = sin_num;
//            double cos = cos_num;
//            double den = denom;
//            double dx = Math.Abs(dirx);
//            double dy = Math.Abs(diry);
//            double sq_hypotenuse = dx * dx + dy * dy;
//            double common_part;
//            double diff_part;
//            double rhs;
//            bool lower_ok;
//            bool upper_ok;

//            if (dy > dx)
//            {
//                Swap(ref dx, ref dy);
//            }
//            // approximate sin = dy / sqrt(sq_hypotenuse)
//            // if ( dy / sqrt(sq_hypotenuse) < n/d )
//            if (dy * dy * d * d < sq_hypotenuse * n * n)
//            {
//                cos = 1.0d;
//                sin = 0.0d;
//                den = 1.0d;
//            }
//            else
//            {
//                double p;
//                double q;
//                double p0 = 0.0d;
//                double q0 = 1.0d;
//                double p1 = 1.0d;
//                double q1 = 1.0d;

//                for (; ; )
//                {
//                    p = p0 + p1;
//                    q = q0 + q1;
//                    sin = (2.0d) * p * q;
//                    den = p * p + q * q;

//                    // sanity check for approximation
//                    //        sin/den < dy/sqrt(hypotenuse) + n/d
//                    //    &&  sin/den > dy/sqrt(hypotenuse) - n/d
//                    // ===    sin/den - n/d  <   dy/sqrt(sq_hypotenuse)
//                    //    &&  sin/den + n/d  >   dy/sqrt(sq_hypotenuse)
//                    // ===    (sin^2 d^2 + n^2 den^2)sq_hypotenuse - 2... < dy^2 d^2 den^2
//                    //    &&  (sin^2 d^2 + n^2 den^2)sq_hypotenuse + 2... > dy^2 d^2 den^2

//                    common_part = (sin * sin * d * d + n * n * den * den) * sq_hypotenuse;
//                    diff_part = 2.0d * n * sin * d * den * sq_hypotenuse;
//                    rhs = dy * dy * d * d * den * den;

//                    upper_ok = (common_part - diff_part < rhs);
//                    lower_ok = (common_part + diff_part > rhs);

//                    if (lower_ok && upper_ok)
//                    {
//                        // if ( (p*p)%2 + (q*q)%2 > double1)
//                        // {
//                        //     sin = p*q;
//                        //     cos = (q*q - p*p)/2;    // exact division
//                        //     den = (p*p + q*q)/2;    // exact division
//                        // }
//                        // else
//                        // {
//                        cos = q * q - p * p;
//                        // }

//                        break;
//                    }
//                    else
//                    {
//                        // if ( dy/sqrt(sq_hypotenuse) < sin/den )
//                        if (dy * dy * den * den < sin * sin * sq_hypotenuse)
//                        {
//                            p1 = p;
//                            q1 = q;
//                        }
//                        else
//                        {
//                            p0 = p;
//                            q0 = q;
//                        }
//                    }
//                } // for(;;)
//            }
//            dx = dirx;
//            dy = diry;

//            if (Math.Abs(dy) > Math.Abs(dx))
//            {
//                Swap(ref sin, ref cos);
//            }

//            if (dx < 0.0d)
//            {
//                cos = -cos;
//            }

//            if (dy < 0.0d)
//            {
//                sin = -sin;
//            }

//            sin_num = sin;
//            cos_num = cos;
//            denom = den;
//        }

//        public static void
//                rational_rotation_approximation(double angle,
//                                          ref double sin_num,  // return
//                                          ref double cos_num,  // return
//                                          ref double denom,    // return
//                                          double eps_num,  // quality_bound
//                                          double eps_den)
//        {
//            double n = eps_num;
//            double d = eps_den;
//            if (n <= 0.0d)
//                throw new Exception();
//            if (d <= 0.0d)
//                throw new Exception();
//            double isin = sin_num;
//            double icos = cos_num;
//            double iden = denom;
//            double dsin = Math.Sin(angle);
//            double dcos = Math.Cos(angle);
//            double dn = n;
//            double dd = d;
//            double eps = dn / dd;
//            dsin = Math.Abs(dsin);
//            dcos = Math.Abs(dcos);
//            double common_part;
//            double diff_part;
//            double os;
//            bool lower_ok;
//            bool upper_ok;
//            bool swapped = false;

//            if (dsin > dcos)
//            {
//                swapped = true;
//                Swap(ref dsin, ref dcos);
//            }
//            if (dsin < eps)
//            {
//                icos = 1;
//                isin = 0;
//                iden = 1;
//            }
//            else
//            {
//                double p;
//                double q;
//                double p0 = 0.0d;
//                double q0 = 1.0d;
//                double p1 = 1.0d;
//                double q1 = 1.0d;

//                for (; ; )
//                {
//                    p = p0 + p1;
//                    q = q0 + q1;
//                    isin = 2.0d * p * q;
//                    iden = p * p + q * q;

//                    // XXX sanity check for approximation
//                    //        sin/den < dsin + n/d
//                    //    &&  sin/den > dsin - n/d
//                    //        sin < dsin * den + n/d * den
//                    //    &&  sin > dsin * den - n/d * den
//                    os = isin;
//                    diff_part = eps * iden;
//                    common_part = dsin * iden;

//                    upper_ok = (common_part - diff_part < os);
//                    lower_ok = (os < common_part + diff_part);

//                    if (lower_ok && upper_ok)
//                    {
//                        // if ( (p*p)%2 + (q*q)%2 > double1)
//                        // {
//                        //     isin = p*q;
//                        //     icos = (q*q - p*p)/2;    // exact division
//                        //     iden = (p*p + q*q)/2;    // exact division
//                        // }
//                        // else
//                        // {
//                        icos = q * q - p * p;
//                        // }

//                        break;
//                    }
//                    else
//                    {
//                        // XXX if ( dsin < sin/den )
//                        if (dsin * iden < isin)
//                        {
//                            p1 = p;
//                            q1 = q;
//                        }
//                        else
//                        {
//                            p0 = p;
//                            q0 = q;
//                        }
//                    }
//                } // for(;;)
//            }

//            if (swapped)
//            {
//                Swap(ref isin, ref icos);
//            }

//            dsin = Math.Sin(angle);
//            dcos = Math.Cos(angle);
//            if (dcos < 0.0d)
//            {
//                icos = -icos;
//            }
//            if (dsin < 0.0d)
//            {
//                isin = -isin;
//            }

//            sin_num = isin;
//            cos_num = icos;
//            denom = iden;
//        }

//        public static bool IsFinite(double x)
//        {
//            return !double.IsInfinity(x) && !double.IsNaN(x);
//        }

//        public static bool isZero(double a)
//        {
//            return Math.Abs(a) < 1e-9;
//        }

//        public static bool IsNear(Vector2 a, Vector2 b)
//        {
//            Vector2 d = (a - b);
//            return (isZero(d.X) && isZero(d.Y));
//        }

//        public static int compare(double a, double b)
//        {
//            return a.CompareTo(b);
//        }

//        public static double
//                    determinant(
//            double a00, double a01,
//            double a10, double a11)
//        {
//            // First compute the det2x2
//            double m01 = a00 * a11 - a10 * a01;
//            return m01;
//        }

//        public static double
//                    determinant(
//            double a00, double a01, double a02,
//            double a10, double a11, double a12,
//            double a20, double a21, double a22)
//        {
//            // First compute the det2x2
//            double m01 = a00 * a11 - a10 * a01;
//            double m02 = a00 * a21 - a20 * a01;
//            double m12 = a10 * a21 - a20 * a11;
//            // Now compute the minors of rank 3
//            double m012 = m01 * a22 - m02 * a12 + m12 * a02;
//            return m012;
//        }

//        public static double
//                    determinant(
//            double a00, double a01, double a02, double a03,
//            double a10, double a11, double a12, double a13,
//            double a20, double a21, double a22, double a23,
//            double a30, double a31, double a32, double a33)
//        {
//            // First compute the det2x2
//            double m01 = a10 * a01 - a00 * a11;
//            double m02 = a20 * a01 - a00 * a21;
//            double m03 = a30 * a01 - a00 * a31;
//            double m12 = a20 * a11 - a10 * a21;
//            double m13 = a30 * a11 - a10 * a31;
//            double m23 = a30 * a21 - a20 * a31;
//            // Now compute the minors of rank 3
//            double m012 = m12 * a02 - m02 * a12 + m01 * a22;
//            double m013 = m13 * a02 - m03 * a12 + m01 * a32;
//            double m023 = m23 * a02 - m03 * a22 + m02 * a32;
//            double m123 = m23 * a12 - m13 * a22 + m12 * a32;
//            // Now compute the minors of rank 4
//            double m0123 = m123 * a03 - m023 * a13 + m013 * a23 - m012 * a33;
//            return m0123;
//        }

//        public static double
//                    determinant(
//            double a00, double a01, double a02, double a03, double a04,
//            double a10, double a11, double a12, double a13, double a14,
//            double a20, double a21, double a22, double a23, double a24,
//            double a30, double a31, double a32, double a33, double a34,
//            double a40, double a41, double a42, double a43, double a44)
//        {
//            // First compute the det2x2
//            double m01 = a10 * a01 - a00 * a11;
//            double m02 = a20 * a01 - a00 * a21;
//            double m03 = a30 * a01 - a00 * a31;
//            double m04 = a40 * a01 - a00 * a41;
//            double m12 = a20 * a11 - a10 * a21;
//            double m13 = a30 * a11 - a10 * a31;
//            double m14 = a40 * a11 - a10 * a41;
//            double m23 = a30 * a21 - a20 * a31;
//            double m24 = a40 * a21 - a20 * a41;
//            double m34 = a40 * a31 - a30 * a41;
//            // Now compute the minors of rank 3
//            double m012 = m12 * a02 - m02 * a12 + m01 * a22;
//            double m013 = m13 * a02 - m03 * a12 + m01 * a32;
//            double m014 = m14 * a02 - m04 * a12 + m01 * a42;
//            double m023 = m23 * a02 - m03 * a22 + m02 * a32;
//            double m024 = m24 * a02 - m04 * a22 + m02 * a42;
//            double m034 = m34 * a02 - m04 * a32 + m03 * a42;
//            double m123 = m23 * a12 - m13 * a22 + m12 * a32;
//            double m124 = m24 * a12 - m14 * a22 + m12 * a42;
//            double m134 = m34 * a12 - m14 * a32 + m13 * a42;
//            double m234 = m34 * a22 - m24 * a32 + m23 * a42;
//            // Now compute the minors of rank 4
//            double m0123 = m123 * a03 - m023 * a13 + m013 * a23 - m012 * a33;
//            double m0124 = m124 * a03 - m024 * a13 + m014 * a23 - m012 * a43;
//            double m0134 = m134 * a03 - m034 * a13 + m014 * a33 - m013 * a43;
//            double m0234 = m234 * a03 - m034 * a23 + m024 * a33 - m023 * a43;
//            double m1234 = m234 * a13 - m134 * a23 + m124 * a33 - m123 * a43;
//            // Now compute the minors of rank 5
//            double m01234 = m1234 * a04 - m0234 * a14 + m0134 * a24 - m0124 * a34 + m0123 * a44;
//            return m01234;
//        }

//        public static double
//                    determinant(
//            double a00, double a01, double a02, double a03, double a04,
//            double a05,
//            double a10, double a11, double a12, double a13, double a14,
//            double a15,
//            double a20, double a21, double a22, double a23, double a24,
//            double a25,
//            double a30, double a31, double a32, double a33, double a34,
//            double a35,
//            double a40, double a41, double a42, double a43, double a44,
//            double a45,
//            double a50, double a51, double a52, double a53, double a54,
//            double a55)
//        {
//            // First compute the det2x2
//            double m01 = a00 * a11 - a10 * a01;
//            double m02 = a00 * a21 - a20 * a01;
//            double m03 = a00 * a31 - a30 * a01;
//            double m04 = a00 * a41 - a40 * a01;
//            double m05 = a00 * a51 - a50 * a01;
//            double m12 = a10 * a21 - a20 * a11;
//            double m13 = a10 * a31 - a30 * a11;
//            double m14 = a10 * a41 - a40 * a11;
//            double m15 = a10 * a51 - a50 * a11;
//            double m23 = a20 * a31 - a30 * a21;
//            double m24 = a20 * a41 - a40 * a21;
//            double m25 = a20 * a51 - a50 * a21;
//            double m34 = a30 * a41 - a40 * a31;
//            double m35 = a30 * a51 - a50 * a31;
//            double m45 = a40 * a51 - a50 * a41;
//            // Now compute the minors of rank 3
//            double m012 = m01 * a22 - m02 * a12 + m12 * a02;
//            double m013 = m01 * a32 - m03 * a12 + m13 * a02;
//            double m014 = m01 * a42 - m04 * a12 + m14 * a02;
//            double m015 = m01 * a52 - m05 * a12 + m15 * a02;
//            double m023 = m02 * a32 - m03 * a22 + m23 * a02;
//            double m024 = m02 * a42 - m04 * a22 + m24 * a02;
//            double m025 = m02 * a52 - m05 * a22 + m25 * a02;
//            double m034 = m03 * a42 - m04 * a32 + m34 * a02;
//            double m035 = m03 * a52 - m05 * a32 + m35 * a02;
//            double m045 = m04 * a52 - m05 * a42 + m45 * a02;
//            double m123 = m12 * a32 - m13 * a22 + m23 * a12;
//            double m124 = m12 * a42 - m14 * a22 + m24 * a12;
//            double m125 = m12 * a52 - m15 * a22 + m25 * a12;
//            double m134 = m13 * a42 - m14 * a32 + m34 * a12;
//            double m135 = m13 * a52 - m15 * a32 + m35 * a12;
//            double m145 = m14 * a52 - m15 * a42 + m45 * a12;
//            double m234 = m23 * a42 - m24 * a32 + m34 * a22;
//            double m235 = m23 * a52 - m25 * a32 + m35 * a22;
//            double m245 = m24 * a52 - m25 * a42 + m45 * a22;
//            double m345 = m34 * a52 - m35 * a42 + m45 * a32;
//            // Now compute the minors of rank 4
//            double m0123 = m012 * a33 - m013 * a23 + m023 * a13 - m123 * a03;
//            double m0124 = m012 * a43 - m014 * a23 + m024 * a13 - m124 * a03;
//            double m0125 = m012 * a53 - m015 * a23 + m025 * a13 - m125 * a03;
//            double m0134 = m013 * a43 - m014 * a33 + m034 * a13 - m134 * a03;
//            double m0135 = m013 * a53 - m015 * a33 + m035 * a13 - m135 * a03;
//            double m0145 = m014 * a53 - m015 * a43 + m045 * a13 - m145 * a03;
//            double m0234 = m023 * a43 - m024 * a33 + m034 * a23 - m234 * a03;
//            double m0235 = m023 * a53 - m025 * a33 + m035 * a23 - m235 * a03;
//            double m0245 = m024 * a53 - m025 * a43 + m045 * a23 - m245 * a03;
//            double m0345 = m034 * a53 - m035 * a43 + m045 * a33 - m345 * a03;
//            double m1234 = m123 * a43 - m124 * a33 + m134 * a23 - m234 * a13;
//            double m1235 = m123 * a53 - m125 * a33 + m135 * a23 - m235 * a13;
//            double m1245 = m124 * a53 - m125 * a43 + m145 * a23 - m245 * a13;
//            double m1345 = m134 * a53 - m135 * a43 + m145 * a33 - m345 * a13;
//            double m2345 = m234 * a53 - m235 * a43 + m245 * a33 - m345 * a23;
//            // Now compute the minors of rank 5
//            double m01234 = m0123 * a44 - m0124 * a34 + m0134 * a24 - m0234 * a14 + m1234 * a04;
//            double m01235 = m0123 * a54 - m0125 * a34 + m0135 * a24 - m0235 * a14 + m1235 * a04;
//            double m01245 = m0124 * a54 - m0125 * a44 + m0145 * a24 - m0245 * a14 + m1245 * a04;
//            double m01345 = m0134 * a54 - m0135 * a44 + m0145 * a34 - m0345 * a14 + m1345 * a04;
//            double m02345 = m0234 * a54 - m0235 * a44 + m0245 * a34 - m0345 * a24 + m2345 * a04;
//            double m12345 = m1234 * a54 - m1235 * a44 + m1245 * a34 - m1345 * a24 + m2345 * a14;
//            // Now compute the minors of rank 6
//            double m012345 = m01234 * a55 - m01235 * a45 + m01245 * a35 - m01345 * a25
//                             + m02345 * a15 - m12345 * a05;
//            return m012345;
//        }

//        public static int
//                    sign_of_determinant(double a00, double a01,
//                               double a10, double a11)
//        {
//            return (a00 * a11).CompareTo(a10 * a01);
//        }

//        public static int
//                    sign_of_determinant(double a00, double a01, double a02,
//                               double a10, double a11, double a12,
//                               double a20, double a21, double a22)
//        {
//            return Math.Sign(determinant(a00, a01, a02,
//                a10, a11, a12,
//                a20, a21, a22));
//        }

//        public static int
//                    sign_of_determinant(
//            double a00, double a01, double a02, double a03,
//            double a10, double a11, double a12, double a13,
//            double a20, double a21, double a22, double a23,
//            double a30, double a31, double a32, double a33)
//        {
//            return Math.Sign(determinant(a00, a01, a02, a03,
//                a10, a11, a12, a13,
//                a20, a21, a22, a23,
//                a30, a31, a32, a33));
//        }

//        public static int
//                    sign_of_determinant(
//            double a00, double a01, double a02, double a03, double a04,
//            double a10, double a11, double a12, double a13, double a14,
//            double a20, double a21, double a22, double a23, double a24,
//            double a30, double a31, double a32, double a33, double a34,
//            double a40, double a41, double a42, double a43, double a44)
//        {
//            return Math.Sign(determinant(a00, a01, a02, a03, a04,
//                a10, a11, a12, a13, a14,
//                a20, a21, a22, a23, a24,
//                a30, a31, a32, a33, a34,
//                a40, a41, a42, a43, a44));
//        }

//        public static int
//                    sign_of_determinant(
//            double a00, double a01, double a02, double a03, double a04,
//            double a05,
//            double a10, double a11, double a12, double a13, double a14,
//            double a15,
//            double a20, double a21, double a22, double a23, double a24,
//            double a25,
//            double a30, double a31, double a32, double a33, double a34,
//            double a35,
//            double a40, double a41, double a42, double a43, double a44,
//            double a45,
//            double a50, double a51, double a52, double a53, double a54,
//            double a55)
//        {
//            return Math.Sign(determinant(a00, a01, a02, a03, a04, a05,
//                a10, a11, a12, a13, a14, a15,
//                a20, a21, a22, a23, a24, a25,
//                a30, a31, a32, a33, a34, a35,
//                a40, a41, a42, a43, a44, a45,
//                a50, a51, a52, a53, a54, a55));
//        }

//        public static bool parallelC2(double l1a, double l1b,
//                                       double l2a, double l2b)
//        {
//            return sign_of_determinant(l1a, l1b, l2a, l2b) == (int)ZERO;
//        }

//        public static bool parallelC2(double s1sx, double s1sy,
//                                       double s1tx, double s1ty,
//                                       double s2sx, double s2sy,
//                                       double s2tx, double s2ty)
//        {
//            return sign_of_determinant(s1tx - s1sx, s1ty - s1sy,
//                s2tx - s2sx, s2ty - s2sy) == (int)ZERO;
//        }

//        public static bool equal_lineC2(double l1a, double l1b, double l1c,
//                                         double l2a, double l2b, double l2c)
//        {
//            if (sign_of_determinant(l1a, l1b, l2a, l2b) != (int)ZERO)
//                return false; // Not parallel.

//            int s1a = Math.Sign(l1a);
//            if (s1a != (int)ZERO)
//                return s1a == Math.Sign(l2a)
//                && sign_of_determinant(l1a, l1c, l2a, l2c) == (int)ZERO;
//            return Math.Sign(l1b) == Math.Sign(l2b)
//            && sign_of_determinant(l1b, l1c, l2b, l2c) == (int)ZERO;
//        }

//        public static int
//                compare_xC2(double px,
//                      double la, double lb, double lc,
//                      double ha, double hb, double hc)
//        {
//            // The abscissa of the intersection point is num/den.
//            double num = determinant(lb, lc, hb, hc);
//            double den = determinant(la, lb, ha, hb);
//            int s = Math.Sign(den);
//            if (s == 0)
//                throw new Exception();
//            return s * compare(px * den, num);
//        }

//        public static int
//                compare_xC2(double la, double lb, double lc,
//                      double h1a, double h1b, double h1c,
//                      double h2a, double h2b, double h2c)
//        {
//            /*
//  double num1 = determinant( lb, lc, h1b, h1c);
//  double den1 = determinant( la, lb, h1a, h1b);
//  double num2 = determinant( lb, lc, h2b, h2c);
//  double den2 = determinant( la, lb, h2a, h2b);
//  Sign s = Sign (Math.Sign(den1) * Math.Sign(den2));
//  CGAL_kernel_assertion( s != Sign.ZERO );
//  return s * sign_of_determinant(num1, num2, den1, den2);
//  */
//            double num1 = determinant(la, lc, h1a, h1c);
//            double num2 = determinant(la, lc, h2a, h2c);
//            double num = determinant(h1a, h1c, h2a, h2c) * lb
//                         + determinant(num1, num2, h1b, h2b);
//            double den1 = determinant(la, lb, h1a, h1b);
//            double den2 = determinant(la, lb, h2a, h2b);
//            return Math.Sign(lb) *
//            Math.Sign(num) *
//            Math.Sign(den1) *
//            Math.Sign(den2);
//        }

//        public static int
//                    compare_xC2(double l1a, double l1b, double l1c,
//                       double h1a, double h1b, double h1c,
//                       double l2a, double l2b, double l2c,
//                       double h2a, double h2b, double h2c)
//        {
//            double num1 = determinant(l1b, l1c, h1b, h1c);
//            double den1 = determinant(l1a, l1b, h1a, h1b);
//            double num2 = determinant(l2b, l2c, h2b, h2c);
//            double den2 = determinant(l2a, l2b, h2a, h2b);
//            int s = Math.Sign(den1) * Math.Sign(den2);
//            if (s == 0)
//                throw new Exception();

//            return s * sign_of_determinant(num1, num2, den1, den2);
//        }

//        public static int
//                    compare_y_at_xC2(double px, double py,
//                            double la, double lb, double lc)
//        {
//            int s = Math.Sign(lb);
//            if (s == 0)
//                throw new Exception();
//            return s * Math.Sign(la * px + lb * py + lc);
//        }

//        public static int
//                    compare_y_at_xC2(double px,
//                            double l1a, double l1b, double l1c,
//                            double l2a, double l2b, double l2c)
//        {
//            int s = Math.Sign(l1b) * Math.Sign(l2b);
//            if (s == 0)
//                throw new Exception();
//            return s * sign_of_determinant(l2a * px + l2c, l2b,
//                l1a * px + l1c, l1b);
//        }

//        public static int
//                    compare_y_at_xC2(double l1a, double l1b, double l1c,
//                            double l2a, double l2b, double l2c,
//                            double ha, double hb, double hc)
//        {
//            int s = Math.Sign(hb) *
//                    sign_of_determinant(l1a, l1b, l2a, l2b);
//            if (s == 0)
//                throw new Exception();
//            return s * sign_of_determinant(l1a, l1b, l1c,
//                l2a, l2b, l2c,
//                ha, hb, hc);
//        }

//        public static int
//                    compare_y_at_xC2(double l1a, double l1b, double l1c,
//                            double l2a, double l2b, double l2c,
//                            double h1a, double h1b, double h1c,
//                            double h2a, double h2b, double h2c)
//        {
//            // The abscissa of the intersection point is num/den.
//            double num = determinant(l1b, l1c, l2b, l2c);
//            double den = determinant(l1a, l1b, l2a, l2b);
//            int s = Math.Sign(h1b) *
//                    Math.Sign(h2b) *
//                    Math.Sign(den);
//            if (s == 0)
//                throw new Exception();
//            return s * sign_of_determinant(h2a * num + h2c * den, h2b,
//                h1a * num + h1c * den, h1b);
//        }

//        public static int
//                    compare_y_at_xC2(double px, double py,
//                            double ssx, double ssy,
//                            double stx, double sty)
//        {
//            // compares the y-coordinates of p and the vertical projection of p on s.
//            // Precondition : p is in the x-range of s.

//            if (!are_ordered(ssx, px, stx))
//                throw new Exception();
//            if (ssx < stx)
//                return orientation(px, py, ssx, ssy, stx, sty);
//            else if (ssx > stx)
//                return orientation(px, py, stx, sty, ssx, ssy);
//            else
//            {
//                if (py < Math.Min(sty, ssy))
//                    return (int)SMALLER;

//                if (py > Math.Max(sty, ssy))
//                    return (int)LARGER;
//                return (int)EQUAL;
//            }
//        }

//        public static int
//                    compare_y_at_x_segment_C2(double px,
//                                     double s1sx, double s1sy,
//                                     double s1tx, double s1ty,
//                                     double s2sx, double s2sy,
//                                     double s2tx, double s2ty)
//        {
//            // compares the y-coordinates of the vertical projections of p on s1 and s2
//            // Precondition : p is in the x-range of s1 and s2.
//            // - if one or two segments are vertical :
//            //   - if the segments intersect, return EQUAL
//            //   - if not, return the obvious SMALLER/LARGER.
//            /*
//            CGAL_kernel_precondition(are_ordered(s1sx, px, s1tx));
//            CGAL_kernel_precondition(are_ordered(s2sx, px, s2tx));
//*/
//            if (s1sx != s1tx && s2sx != s2tx)
//            {
//                double s1stx = s1sx - s1tx;
//                double s2stx = s2sx - s2tx;

//                return compare(s1sx, s1tx) *
//                compare(s2sx, s2tx) *
//                compare(-(s1sx - px) * (s1sy - s1ty) * s2stx,
//                    (s2sy - s1sy) * s2stx * s1stx
//                    - (s2sx - px) * (s2sy - s2ty) * s1stx);
//            }
//            else
//            {
//                if (s1sx == s1tx)
//                { // s1 is vertical
//                    int c1, c2;
//                    c1 = compare_y_at_xC2(px, s1sy, s2sx, s2sy, s2tx, s2ty);
//                    c2 = compare_y_at_xC2(px, s1ty, s2sx, s2sy, s2tx, s2ty);
//                    if (c1 == c2)
//                        return c1;
//                    return (int)EQUAL;
//                }
//                // s2 is vertical
//                int c3, c4;
//                c3 = compare_y_at_xC2(px, s2sy, s1sx, s1sy, s1tx, s1ty);
//                c4 = compare_y_at_xC2(px, s2ty, s1sx, s1sy, s1tx, s1ty);
//                if (c3 == c4)
//                    return -c3;
//                return (int)EQUAL;
//            }
//        }

//        public static bool
//                equal_directionC2(double dx1, double dy1,
//                            double dx2, double dy2)
//        {
//            return Math.Sign(dx1) == Math.Sign(dx2)
//            && Math.Sign(dy1) == Math.Sign(dy2)
//            && sign_of_determinant(dx1, dy1, dx2, dy2) == (int)ZERO;
//        }

//        public static int
//                compare_angle_with_x_axisC2(double dx1, double dy1,
//                                      double dx2, double dy2)
//        {
//            // angles are in [-pi,pi], and the angle between Ox and d1 is compared
//            // with the angle between Ox and d2
//            int quadrant_1 = (dx1 >= 0) ? (dy1 >= 0 ? 1 : 4)
//                : (dy1 >= 0 ? 2 : 3);
//            int quadrant_2 = (dx2 >= 0) ? (dy2 >= 0 ? 1 : 4)
//                : (dy2 >= 0 ? 2 : 3);
//            // We can't use compare(quadrant_1,quadrant_2) because in case
//            // of tie, we need additional computation
//            if (quadrant_1 > quadrant_2)
//                return (int)LARGER;
//            else if (quadrant_1 < quadrant_2)
//                return (int)SMALLER;
//            return -sign_of_determinant(dx1, dy1, dx2, dy2);
//        }

//        public static int
//                compare_slopesC2(double l1a, double l1b, double l2a, double l2b)
//        {
//            if ((int)ZERO == (l1a))  // l1 is horizontal
//                return (int)ZERO == (l2b) ? (int)(SMALLER)
//                    : Math.Sign(l2a) * Math.Sign(l2b);
//            if ((int)ZERO == (l2a)) // l2 is horizontal
//                return (int)ZERO == (l1b) ? (int)(LARGER)
//                    : -Math.Sign(l1a) * Math.Sign(l1b);
//            if ((int)ZERO == (l1b))
//                return (int)ZERO == (l2b) ? (int)EQUAL : (int)LARGER;
//            if ((int)ZERO == (l2b))
//                return (int)SMALLER;
//            int l1_sign = -Math.Sign(l1a) * Math.Sign(l1b);
//            int l2_sign = -Math.Sign(l2a) * Math.Sign(l2b);

//            if (l1_sign < l2_sign)
//                return (int)SMALLER;
//            if (l1_sign > l2_sign)
//                return (int)LARGER;

//            if (l1_sign > (int)ZERO)
//                return compare(Math.Abs(l1a * l2b),
//                    Math.Abs(l2a * l1b));

//            return compare(Math.Abs(l2a * l1b),
//                Math.Abs(l1a * l2b));
//        }

//        public static int
//                compare_slopesC2(double s1_src_x, double s1_src_y, double s1_tgt_x,
//                           double s1_tgt_y, double s2_src_x, double s2_src_y,
//                           double s2_tgt_x, double s2_tgt_y)
//        {
//            int cmp_x1;
//            int cmp_x2;

//            int cmp_y1 = compare(s1_src_y, s1_tgt_y);
//            if (cmp_y1 == (int)EQUAL)
//            { // horizontal
//                cmp_x2 = compare(s2_src_x, s2_tgt_x);

//                if (cmp_x2 == (int)EQUAL)
//                    return (int)SMALLER;
//                return -Math.Sign(s2_src_y - s2_tgt_y) * Math.Sign(s2_src_x - s2_tgt_x);
//            }

//            int cmp_y2 = compare(s2_src_y, s2_tgt_y);
//            if (cmp_y2 == (int)EQUAL)
//            {
//                cmp_x1 = compare(s1_src_x, s1_tgt_x);

//                if (cmp_x1 == (int)EQUAL)
//                    return (int)LARGER;
//                return Math.Sign(s1_src_y - s1_tgt_y) * Math.Sign(s1_src_x - s1_tgt_x);
//            }

//            cmp_x1 = compare(s1_src_x, s1_tgt_x);
//            cmp_x2 = compare(s2_src_x, s2_tgt_x);

//            if (cmp_x1 == (int)EQUAL)
//                return cmp_x2 == (int)EQUAL ? (int)EQUAL : (int)LARGER;

//            if (cmp_x2 == (int)EQUAL)
//                return (int)SMALLER;

//            double s1_xdiff = s1_src_x - s1_tgt_x;
//            double s1_ydiff = s1_src_y - s1_tgt_y;
//            double s2_xdiff = s2_src_x - s2_tgt_x;
//            double s2_ydiff = s2_src_y - s2_tgt_y;
//            int s1_sign = Math.Sign(s1_ydiff) * Math.Sign(s1_xdiff);
//            int s2_sign = Math.Sign(s2_ydiff) * Math.Sign(s2_xdiff);

//            if (s1_sign < s2_sign)
//                return (int)SMALLER;
//            if (s1_sign > s2_sign)
//                return (int)LARGER;

//            if (s1_sign > (int)ZERO)
//                return compare(Math.Abs(s1_ydiff * s2_xdiff),
//                    Math.Abs(s2_ydiff * s1_xdiff));

//            return compare(Math.Abs(s2_ydiff * s1_xdiff),
//                Math.Abs(s1_ydiff * s2_xdiff));
//        }

//        public static int
//                compare_lexicographically_xyC2(double px, double py,
//                                         double qx, double qy)
//        {
//            int c = compare(px, qx);
//            return (c != (int)EQUAL) ? c : compare(py, qy);
//        }

//        public static int orientation(Vector2 p, Vector2 q, Vector2 r)
//        {
//            return orientation(p.X, p.Y, q.X, q.Y, r.X, r.Y);
//        }

//        public static int
//                orientation(double px, double py,
//                      double qx, double qy,
//                      double rx, double ry)
//        {
//            return sign_of_determinant(qx - px, qy - py, rx - px, ry - py);
//        }

//        public static int
//                orientation(double ux, double uy, double vx, double vy)
//        {
//            return sign_of_determinant(ux, uy, vx, vy);
//        }

//        public static int
//                angleC2(double ux, double uy,
//                  double vx, double vy)
//        {
//            return (Math.Sign(ux * vx + uy * vy));
//        }

//        public static int
//                angleC2(double px, double py,
//                  double qx, double qy,
//                  double rx, double ry)
//        {
//            return (Math.Sign((px - qx) * (rx - qx) + (py - qy) * (ry - qy)));
//        }

//        public static int
//                angleC2(double px, double py,
//                  double qx, double qy,
//                  double rx, double ry,
//                  double sx, double sy)
//        {
//            return (Math.Sign((px - qx) * (rx - sx) + (py - qy) * (ry - sy)));
//        }

//        public static bool
//            collinear_are_ordered_along_lineC2(in Point2 p, in Point2 q, in Point2 r) =>
//            collinear_are_ordered_along_lineC2(p.X, p.Y, q.X, q.Y, r.X, r.Y);

//        public static bool
//                collinear_are_ordered_along_lineC2(double px, double py,
//                                             double qx, double qy,
//                                             double rx, double ry)
//        {
//            if (px < qx)
//                return !(rx < qx);
//            if (qx < px)
//                return !(qx < rx);
//            if (py < qy)
//                return !(ry < qy);
//            if (qy < py)
//                return !(qy < ry);
//            return true; // p==q
//        }

//        public static bool
//                collinear_are_strictly_ordered_along_lineC2(double px, double py,
//                                                      double qx, double qy,
//                                                      double rx, double ry)
//        {
//            if (px < qx)
//                return (qx < rx);
//            if (qx < px)
//                return (rx < qx);
//            if (py < qy)
//                return (qy < ry);
//            if (qy < py)
//                return (ry < qy);
//            return false;
//        }

//        public static int
//                side_of_oriented_circleC2(double px, double py,
//                                    double qx, double qy,
//                                    double rx, double ry,
//                                    double tx, double ty)
//        {
//            //  sign_of_determinant(px, py, px*px + py*py, 1,
//            //                         qx, qy, qx*qx + qy*qy, 1,
//            //                         rx, ry, rx*rx + ry*ry, 1,
//            //                         tx, ty, tx*tx + ty*ty, 1);
//            // We first translate so that p is the new origin.
//            double qpx = qx - px;
//            double qpy = qy - py;
//            double rpx = rx - px;
//            double rpy = ry - py;
//            double tpx = tx - px;
//            double tpy = ty - py;
//            // The usual 3x3 formula can be simplified a little bit to a 2x2.
//            //         - sign_of_determinant(qpx, qpy, square(qpx) + square(qpy),
//            //                                  rpx, rpy, square(rpx) + square(rpy),
//            //                                  tpx, tpy, square(tpx) + square(tpy)));
//            return sign_of_determinant(qpx * tpy - qpy * tpx, tpx * (tx - qx) + tpy * (ty - qy),
//                qpx * rpy - qpy * rpx, rpx * (rx - qx) + rpy * (ry - qy));
//        }

//        public static int
//                side_of_bounded_circleC2(double px, double py,
//                                   double qx, double qy,
//                                   double rx, double ry,
//                                   double tx, double ty)
//        {
//            return (side_of_oriented_circleC2(px, py, qx, qy, rx, ry, tx, ty)
//            * orientation(px, py, qx, qy, rx, ry));
//        }

//        public static int
//                side_of_bounded_circleC2(double px, double py,
//                                   double qx, double qy,
//                                   double tx, double ty)
//        {
//            // Returns whether T lies inside or outside the circle which diameter is PQ.
//            return
//                 compare((tx - px) * (qx - tx), (ty - py) * (ty - qy));
//        }

//        public static int
//                cmp_dist_to_pointC2(double px, double py,
//                              double qx, double qy,
//                              double rx, double ry)
//        {
//            return compare(squared_distanceC2(px, py, qx, qy),
//                squared_distanceC2(px, py, rx, ry));
//        }

//        public static bool
//                is_larger_dist_to_pointC2(double px, double py,
//                                    double qx, double qy,
//                                    double rx, double ry)
//        {
//            return cmp_dist_to_pointC2(px, py, qx, qy, rx, ry) == (int)LARGER;
//        }

//        public static bool
//                si_smaller_dist_to_pointC2(double px, double py,
//                                     double qx, double qy,
//                                     double rx, double ry)
//        {
//            return cmp_dist_to_pointC2(px, py, qx, qy, rx, ry) == (int)SMALLER;
//        }

//        public static int
//                cmp_signed_dist_to_directionC2(double la, double lb,
//                                         double px, double py,
//                                         double qx, double qy)
//        {
//            return compare(scaled_distance_to_directionC2(la, lb, px, py),
//                scaled_distance_to_directionC2(la, lb, qx, qy));
//        }

//        public static bool
//                is_larger_signed_dist_to_directionC2(double la, double lb,
//                                               double px, double py,
//                                               double qx, double qy)
//        {
//            return cmp_signed_dist_to_directionC2(la, lb, px, py, qx, qy) == (int)LARGER;
//        }

//        public static bool
//                is_smaller_signed_dist_to_directionC2(double la, double lb,
//                                                double px, double py,
//                                                double qx, double qy)
//        {
//            return cmp_signed_dist_to_directionC2(la, lb, px, py, qx, qy) == (int)SMALLER;
//        }

//        public static int
//                cmp_signed_dist_to_lineC2(double px, double py,
//                                    double qx, double qy,
//                                    double rx, double ry,
//                                    double sx, double sy)
//        {
//            return compare(scaled_distance_to_lineC2(px, py, qx, qy, rx, ry),
//                scaled_distance_to_lineC2(px, py, qx, qy, sx, sy));
//        }

//        public static bool
//                is_larger_signed_dist_to_lineC2(double px, double py,
//                                          double qx, double qy,
//                                          double rx, double ry,
//                                          double sx, double sy)
//        {
//            return cmp_signed_dist_to_lineC2(px, py, qx, qy, rx, ry, sx, sy) == (int)LARGER;
//        }

//        public static bool is_smaller_signed_dist_to_lineC2(double px, double py,
//                                                double qx, double qy,
//                                                double rx, double ry,
//                                                double sx, double sy)
//        {
//            return cmp_signed_dist_to_lineC2(px, py, qx, qy, rx, ry, sx, sy) == (int)SMALLER;
//        }

//        public static int
//                side_of_oriented_lineC2(double a, double b, double c,
//                                  double x, double y)
//        {
//            return Math.Sign(a * x + b * y + c);
//        }
    }
}