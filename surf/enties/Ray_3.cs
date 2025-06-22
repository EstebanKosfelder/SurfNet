
namespace SurfNet
{

    using static Mathex;
    using static DebugLog;
    public class Ray3
    {
        public Ray3(in Point_3 sp, in Point_3 secondp)
        {
            Source = sp;
            SecondPoint = secondp;
        }

        public Ray3(in Point_3 sp, in Vector_3 vector) : this(sp, new Point_3(sp.X + vector.X, sp.Y + vector.Y, sp.Z + vector.Z)) { }

   //     public Ray3(in Point_3 sp, in Direction2 d) : this(sp, d.to_vector()) { }

        public Point_3 Source;
        public Point_3 SecondPoint;

        public Point_3 source() => Source;
        public Point_3 second_point() => SecondPoint;

        public bool is_degenerate() => Source == SecondPoint;

        //        translate
        //        Point_2
        //   operator()( const Point_2& p, const Vector_2& v) const
        //    {
        //      Construct_point_2 construct_point_2;
        //      return construct_point_2(p.x() + v.x(), p.y() + v.y());
        //    }

        //    Point_2
        //    operator()( const Origin& , const Vector_2& v) const
        //    {
        //      Construct_point_2 construct_point_2;
        //      return construct_point_2(v.x(), v.y());
        //}

        public Ray3(Point_3 p, Vector_3 v) : this(p, new Point_3(p.X + v.X, p.Y + v.Y, p.Z + v.Z)) { }

        //public Ray3(Point_3 p, Direction2 d) : this(p, d.to_vector()) { }

        //public Ray3(Point_3 p, Line2 l) : this(p, l.to_vector()) { }

        public Ray3(in Ray3 r) : this(r.Source, r.SecondPoint) { }

        /*
        : RRay_2(typename R::Construct_ray_2()(Return_base_tag(), sp, secondp)) {}

      Ray_2(const Point_2 &sp, const Direction_2 &d)
        : RRay_2(typename R::Construct_ray_2()(Return_base_tag(), sp, d)) {}

      Ray_2(const Point_2 &sp, const Vector_2 &v)
        : RRay_2(typename R::Construct_ray_2()(Return_base_tag(), sp, v)) {}

      Ray_2(const Point_2 &sp, const Linel)
        : RRay_2(typename R::Construct_ray_2()(Return_base_tag(), sp, l)) {}
        */

        /*
  decltype(auto)
 */

        public Point_3 point(double i)
        {
            assert(i >= 0);

            if (i == (0)) return source();
            if (i == (1)) return second_point();

            return Source + (SecondPoint - Source) * i;
        }
        public Point_3 start() => source();

        public bool is_horizontal() => are_near(Source.Y, SecondPoint.Y);
        public bool is_vertical() => are_near(Source.Y, SecondPoint.Y);

        public Direction2 direction() => new Direction2(to_vector());

        public Vector2 to_vector() => new Vector2(SecondPoint.X - Source.X, SecondPoint.Y - Source.Y);

        // public  bool has_on(Point_2 p)
        //  {
        //return are_near( p, source()) || collinear2(source(), p, second_point()) &&
        //       Direction_2(construct_vector(source(), p)) == direction());
        //  }

        //public bool collinear_has_on(Point_3 p)
        //{
        //    Point_3 source = Source;
        //    Point_3 second = SecondPoint;
        //    switch ((CompareResultEnum)compare_x(source, second))
        //    {
        //        case CompareResultEnum.SMALLER:
        //            return (CompareResultEnum)compare_x(Source, p) != CompareResultEnum.LARGER;

        //        case CompareResultEnum.LARGER:
        //            return (CompareResultEnum)compare_x(p, Source) != CompareResultEnum.LARGER;

        //        default:
        //            switch ((CompareResultEnum)compare_y(source, second))
        //            {
        //                case CompareResultEnum.SMALLER:
        //                    return (CompareResultEnum)compare_y(source, p) != CompareResultEnum.LARGER;

        //                case CompareResultEnum.LARGER:
        //                    return (CompareResultEnum)compare_y(p, source) != CompareResultEnum.LARGER;

        //                default:
        //                    return true; // p == source
        //            }
        //    } // switch
        //}

        //Ray_2
        //opposite() const
        //{
        //  return Ray_2( source(), - direction() );
        //}

        //public Line2 supporting_line()
        //{
        //    return new Line2(source(), second_point());
        //}

    }

}