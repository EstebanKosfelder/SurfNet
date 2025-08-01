namespace SurfNetPolyBool
{
	using SurfNet;
	public static class Point2Utils
    {
        public static bool PointAboveOrOnLine(Point2 point, Point2 left, Point2 right)
        {
            return (right.X - left.X) * (point.Y - left.Y) - (right.Y - left.Y) * (point.X - left.X) >= -Mathex.EPS;
        }

        public static bool PointBetween(Point2 point, Point2 left, Point2 right)
        {
            // p must be collinear with left->right
            // returns false if p == left, p == right, or left == right
            double dPyLy = point.Y - left.Y;
            double dRxLx = right.X - left.X;
            double dPxLx = point.X - left.X;
            double dRyLy = right.Y - left.Y;

            double dot = dPxLx * dRxLx + dPyLy * dRyLy;

            if (dot < Mathex.EPS)
            {
                return false;
            }

            double sqlen = dRxLx * dRxLx + dRyLy * dRyLy;
            if (dot - sqlen > -Mathex.EPS)
            {
                return false;
            }

            return true;
        }

        private static bool PointsSameX(Point2 point1, Point2 point2)
        {
            return Math.Abs(point1.X - point2.X) < Mathex.EPS;
        }

        private static bool PointsSameY(Point2 point1, Point2 point2)
        {
            return Math.Abs(point1.Y - point2.Y) < Mathex.EPS;
        }

        public static bool PointsSame(Point2 point1, Point2 point2)
        {
            return PointsSameX(point1, point2) && PointsSameY(point1, point2);
        }

        public static int PointsCompare(Point2 point1, Point2 point2)
        {
            if (PointsSameX(point1, point2))
            {
                return PointsSameY(point1, point2) ? 0 : (point1.Y < point2.Y ? -1 : 1);
            }
            return point1.X < point2.X ? -1 : 1;
        }

        public static bool PointsCollinear(Point2 pt1, Point2 pt2, Point2 pt3)
        {
            var dx1 = pt1.X - pt2.X;
            var dy1 = pt1.Y - pt2.Y;
            var dx2 = pt2.X - pt3.X;
            var dy2 = pt2.Y - pt3.Y;
            return Math.Abs(dx1 * dy2 - dx2 * dy1) < Mathex.EPS;
        }

        public static IntersectionPoint LinesIntersect(Point2 a0, Point2 a1, Point2 b0, Point2 b1)
        {
            double adx = a1.X - a0.X;
            double ady = a1.Y - a0.Y;
            double bdx = b1.X - b0.X;
            double bdy = b1.Y - b0.Y;

            double axb = adx * bdy - ady * bdx;

            if (Math.Abs(axb) < Mathex.EPS)
            {
                return null;
            }

            double dx = a0.X - b0.X;
            double dy = a0.Y - b0.Y;

            double a = (bdx * dy - bdy * dx) / axb;
            double b = (adx * dy - ady * dx) / axb;

            return new IntersectionPoint(CalcAlongUsingValue(a), CalcAlongUsingValue(b), new Point2(a0.X + a * adx, a0.Y + a * ady));
        }

        private static int CalcAlongUsingValue(double value)
        {
            if (value <= -Mathex.EPS)
            {
                return -2;
            }
            else if (value < Mathex.EPS)
            {
                return -1;
            }
            else if (value - 1 <= -Mathex.EPS)
            {
                return 0;
            }
            else if (value - 1 < Mathex.EPS)
            {
                return 1;
            }
            else
            {
                return 2;
            }
        }
    }
}