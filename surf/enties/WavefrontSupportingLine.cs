﻿namespace SurfNet
{

    public class WavefrontSupportingLine
    {
        // private  using Transformation = SurfNet.Aff_transformation_2<Kernel>;

        public Line2 l;

        public double weight;

        /// <summary>
        /// arbitrary length
        /// </summary>
        public Vector2 normal_direction;

        /// <summary>
        /// unit length
        /// </summary>
        public Vector2 normal_unit;

        /// <summary>
        /// weighted
        /// </summary>
        public Vector2 normal;

        public WavefrontSupportingLine(Point2 u, Point2 v, double p_weight = 0.0) : this(new Line2(u, v), p_weight)
        {
        }

        public WavefrontSupportingLine(Line2 p_l, double p_weight = 0.0)
        {
            this.l = p_l;
            this.weight = p_weight;
            normal_direction = l.to_vector().perpendicular(OrientationEnum.COUNTERCLOCKWISE);
            normal_unit = normal_direction / Mathex.sqrt(normal_direction.squared_length());
            normal = (normal_unit * weight);
        }

        public Line2 line_at_one()
        {

            // Traslada la línea por un vector
            var result = l.Translate(normal);
            return result;
        }



        public Line2 line_at(double now)
        {
            return l.Translate(normal * now);
        }
    };

    public class WavefrontSupportingLineList : List<WavefrontSupportingLine>
    { };
}