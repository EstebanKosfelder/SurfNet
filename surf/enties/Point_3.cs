
namespace SurfNet
{
    public struct Point_3
    {
        public static readonly Point_3 NaN = new Point_3(double.NaN, double.NaN, double.NaN);
        public static readonly Point_3 ORIGIN = new Point_3(0.0, 0.0, 0.0);

        //public static explicit operator Point_2(Vector_2 v) => new Point_2(v.X, v.Y);
        //public static explicit operator Vector_2(Point_2 v) => new Vector_2(v.X, v.Y);

        public Point_3(double x, double y, double z)
        {
            this.X = x; this.Y = y; this.Z = z;
        }

        public bool IsNaN()=>X==double.NaN||Y== double.NaN ||Z== double.NaN;
        public static bool operator ==(Point_3 lhs, Point_3 rhs) { return lhs.Equals(rhs); }

        public static bool operator !=(Point_3 lhs, Point_3 rhs) { return !lhs.Equals(rhs); }


        public static Point_3 operator +(Point_3 a, Point_3 b)
        {
            return new Point_3() { X = a.X + b.X, Y = a.Y + b.Y, Z = a.Z + b.Z };
        }

        public static Point_3 operator -(Point_3 a, Point_3 b)
        {
            return new Point_3() { X = a.X - b.X, Y = a.Y - b.Y, Z = a.Z - b.Z };
        }

        public static Point_3 operator *(double a, Point_3 b)
        {
            return new Point_3() { X = a * b.X, Y = a * b.Y, Z = a * b.Z };
        }

        public static Point_3 operator *(Point_3 b, double a)
        {
            return new Point_3() { X = a * b.X, Y = a * b.Y, Z = a * b.Z };
        }

        public double X;
        public double Y;

        public double Z;

        public override bool Equals(object? obj)
        {
            return obj is Point_3 point &&
                   X == point.X &&
                   Y == point.Y &&
                   Z == point.Z;
        }


    }

    public struct Vector_3
    {
        public static readonly Point_3 NaN = new Point_3(double.NaN, double.NaN, double.NaN);
        public static readonly Point_3 ORIGIN = new Point_3(0.0, 0.0, 0.0);

        //public static explicit operator Point_2(Vector_2 v) => new Point_2(v.X, v.Y);
        //public static explicit operator Vector_2(Point_2 v) => new Vector_2(v.X, v.Y);

        public Vector_3(double x, double y, double z)
        {
            this.X = x; this.Y = y; this.Z = z;
        }


        public static bool operator ==(Vector_3 lhs, Vector_3 rhs) { return lhs.Equals(rhs); }

        public static bool operator !=(Vector_3 lhs, Vector_3 rhs) { return !lhs.Equals(rhs); }


        public static Vector_3 operator +(Vector_3 a, Vector_3 b)
        {
            return new Vector_3() { X = a.X + b.X, Y = a.Y + b.Y, Z = a.Z + b.Z };
        }

        public static Vector_3 operator -(Vector_3 a, Vector_3 b)
        {
            return new Vector_3() { X = a.X - b.X, Y = a.Y - b.Y, Z = a.Z - b.Z };
        }

        public static Vector_3 operator *(double a, Vector_3 b)
        {
            return new Vector_3() { X = a * b.X, Y = a * b.Y, Z = a * b.Z };
        }

        public static Vector_3 operator *(Vector_3 b, double a)
        {
            return new Vector_3() { X = a * b.X, Y = a * b.Y, Z = a * b.Z };
        }

        public double X;
        public double Y;

        public double Z;

        public override bool Equals(object? obj)
        {
            return obj is Vector_3 point &&
                   X == point.X &&
                   Y == point.Y &&
                   Z == point.Z;
        }


    }

}