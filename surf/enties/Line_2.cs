﻿namespace SurfNet
{


    //public enum LineLineIntersectionResult
    //{
    //    NO_INTERSECTION = 0,
    //    POINT = 1,
    //    LINE = 2,
    //}
    //public class LineLineIntersector
    //{
    //    public LineLineIntersector(Line2 one, Line2 other)
    //    {
    //        this.one = one;
    //        this.other = other;
    //        this.point = Vector2.NaN;
    //        this.line = new Line2(0, 0, 0);
    //    }

    //    public Line2 one;
    //    public Line2 other;
    //    public Vector2 point;
    //    public Line2 line;

    //    public LineLineIntersectionResult intersection_type()
    //    {
    //        /// solve intersection

    //        // == 3 possible outcomes in 2D:
    //        //
    //        // 0. overlapping lines - always intersecting in a line
    //        // 1. crossing - point2
    //        // 2. parallel - no intersection

    //        (double a1, double b1, double c1) = (one.A, one.B, one.C);
    //        (double a2, double b2, double c2) = (other.C, other.B, other.C);
    //        var denom = a1 * b2 - a2 * b1;
    //        //        print(f'denom {denom}')
    //        // if denom == 0: // ^FIXME: use near_zero ?
    //        if (Mathex.IsZero(denom))
    //        {
    //            var x1 = a1 * c2 - a2 * c1;
    //            var x2 = b1 * c2 - b2 * c1;
    //            //            print(f'cross1 {x1}')
    //            //            print(f'cross2 {x2}')
    //            //// if (x1 == 0) and (x2 == 0): // ^FIXME: use near_zero ?
    //            if (Mathex.IsZero(x1) && Mathex.IsZero(x2))
    //            {
    //                // overlapping lines, always intersecting in this configuration
    //                this.line = this.one;
    //                return LineLineIntersectionResult.LINE;
    //            }
    //            else
    //            {
    //                // parallel lines, but not intersecting in this configuration
    //                return LineLineIntersectionResult.NO_INTERSECTION;
    //            }
    //        }
    //        else
    //        {
    //            // crossing lines
    //            var x = b1 * c2 - b2 * c1;
    //            var y = a2 * c1 - a1 * c2;
    //            var w = denom;
    //            var xw = x / w;
    //            var yw = y / w;
    //            var point = new Vector2(xw, yw);
    //            return LineLineIntersectionResult.POINT;
    //        }
    //    }
    //}
}