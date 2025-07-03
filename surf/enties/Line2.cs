using static SurfNet.Mathex;
using FT = double;

namespace SurfNet
{
    public struct Line2
    {

        public override string ToString()
        {
            return $" a:{A} b:{B} c:{C} n:{to_vector().Normal()} ";
        }
        public Line2(FT a, FT b, FT c) { A = a; B = b; C = c; }

        public Line2(Point2 start, Point2 end)
        {
            line_from_pointsC2(start.X, start.Y, end.X, end.Y, out A, out B, out C);
        }
        public Line2(Segment2 s) : this(s.Source, s.Target) { }

        public FT A;
        public FT B;
        public FT C;

        public static bool operator !=(Line2 a, Line2 b)
        {
            return !(a == b);
        }

        public static bool operator ==(Line2 a, Line2 b)
        {
            return equal_lineC2(a.a(), a.b(), a.c(), b.a(), b.b(), b.c());
        }

        public FT a() => A;

        public FT b() => B;

        public FT c() => C;
        public FT Distance(Point2 point) => A * point.X + B * point.Y + C;

      

        public override bool Equals(object? obj)
        {
            if (obj == null || !(obj is Line2 line))
            {
                return false;
            }

            return equal_lineC2(A, B, C, line.a(), line.b(), line.c());
        }

        public override int GetHashCode()
        {
            var hashCode = new HashCode();

            hashCode.Add(A);
            hashCode.Add(B);
            hashCode.Add(C);
            return hashCode.ToHashCode();
        }

        public Point2 line_project_point(FT px, FT py)
        {
            if (are_near(A, 0)) // horizontal line
            {
                return new Point2(px, -C / B);
            }
            else if (are_near(B, 0)) // vertical line
            {
                return new Point2(-C / A, py);
            }
            else
            {
                FT a2 = A * A;
                FT b2 = B * B;
                FT d = a2 + b2;
                return new Point2((b2 * px - A * B * py - A * C) / d, (-A * B * px + a2 * py - B * C) / d);
            }
        }

        public int SideOfOriented(Point2 point)
        {
            return side_of(this.Distance(point));
        }

        public Vector2 to_vector() => new Vector2(b(), -a());

        internal static Line2 perpendicular_through_point(FT la, FT lb, FT px, FT py)
        {
            return new Line2(-lb, la, lb * px - la * py);
        }
        public Line2 Perpendicular(Point2 p)
        {
         
            
            return perpendicular_through_point(A,B,p.X,p.Y);
        }
        public bool HasOnPositiveSide(Point2 p, double epsilon = Mathex.EPS)
        {
            double val = A * p.X + B * p.Y + C;
            return val > epsilon; // Valor positivo → lado positivo
        }

        private static bool equal_lineC2(FT l1a, FT l1b, FT l1c, FT l2a, FT l2b, FT l2c)
        {
            if (Mathex.sign_of_determinant(l1a, l1b, l2a, l2b) != 0)
                return false; // Not parallel.
            var s1a = Mathex.is_zero(l1a) ? 0 : sign(l1a);
            if (s1a != 0) return s1a == sign(l2a) && Mathex.sign_of_determinant(l1a, l1c, l2a, l2c) == 0;
            return (Mathex.is_zero(l1a) ? 0 : sign(l1b)) == (Mathex.is_zero(l2b) ? 0 : Mathex.sign(l2b)) && Mathex.sign_of_determinant(l1b, l1c, l2b, l2c) == 0;
        }

        public Line2 Translate(Vector2 vector)
        {
            // Ajustar el término independiente
            return new Line2(this.A, this.B, this.C - (A * vector.X + B * vector.Y));
        }

        public bool has_on_boundary(Point2 p)
        {
            var o = side_of_oriented_lineC2(A, B, C, p.X, p.Y);
            return o == OrientationEnum.ON_ORIENTED_BOUNDARY;
        }
        public bool has_on(Point2 p) => has_on_boundary(p);

        public bool is_vertical()
        {
            return this.B == 0.0;
        }

        public Point2 point()
        {
            if ( is_vertical() ) return new Point2(-C / A, 0.0);
            return new Point2(0, -(C / B));
        }

        public Line2 Opposite()
        {
            return new Line2(-A, -B, -C);
        }
    

    }
} //namespace SurfNet