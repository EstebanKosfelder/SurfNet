﻿using System.Diagnostics;


namespace SurfNet
{
#if USE_Z
    [DebuggerDisplay("ID {ID} [{X}, {Y}, {Z}]")]
#else
    [DebuggerDisplay("ID {ID} [{X}, {Y}]")]
#endif
    public class DcelPoint : IComparable<DcelPoint>, IEquatable<DcelPoint>
    {
        public int id;
        public int label;

        public double x;
        public double y;
#if USE_Z
        public double z;
#endif

        /// <summary>
        /// Initializes a new instance of the <see cref="DcelPoint" /> class.
        /// </summary>
        public DcelPoint()
            : this(0.0, 0.0, 0)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="DcelPoint" /> class.
        /// </summary>
        /// <param name="x">The x coordinate.</param>
        /// <param name="y">The y coordinate.</param>
        public DcelPoint(double x, double y)
            : this(x, y, 0)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="DcelPoint" /> class.
        /// </summary>
        /// <param name="x">The x coordinate.</param>
        /// <param name="y">The y coordinate.</param>
        /// <param name="label">The point label.</param>
        public DcelPoint(double x, double y, int label)
        {
            this.x = x;
            this.y = y;
            this.label = label;
        }

        #region Public properties

        /// <summary>
        /// Gets or sets the vertex id.
        /// </summary>
        public int ID
        {
            get { return id; }
            set { id = value; }
        }

        /// <summary>
        /// Gets or sets the vertex x coordinate.
        /// </summary>
        public double X
        {
            get { return x; }
            set { x = value; }
        }

        /// <summary>
        /// Gets or sets the vertex y coordinate.
        /// </summary>
        public double Y
        {
            get { return y; }
            set { y = value; }
        }

#if USE_Z
        /// <summary>
        /// Gets or sets the vertex z coordinate.
        /// </summary>
        public double Z
        {
            get { return z; }
            set { z = value; }
        }
#endif

        /// <summary>
        /// Gets or sets a general-purpose label.
        /// </summary>
        /// <remarks>
        /// This is used for the vertex boundary mark.
        /// </remarks>
        public int Label
        {
            get { return label; }
            set { label = value; }
        }

        #endregion

        #region Overriding Equals() and == Operator

        /// <inheritdoc />
        public static bool operator ==(DcelPoint a, DcelPoint b)
        {
            if (a is null)
            {
                // If one is null, but not both, return false.
                return b is null;
            }

            // If both are same instance, return true.
            if (ReferenceEquals(a, b))
            {
                return true;
            }

            return a.Equals(b);
        }

        /// <inheritdoc />
        public static bool operator !=(DcelPoint a, DcelPoint b)
        {
            return !(a == b);
        }

        /// <inheritdoc />
        public override bool Equals(object obj) => Equals(obj as DcelPoint);

        /// <inheritdoc />
        public bool Equals(DcelPoint p)
        {
            // If object is null return false.
            if (p is null)
            {
                return false;
            }

            // Return true if the fields match.
            return (x == p.x) && (y == p.y);
        }

        #endregion

        /// <inheritdoc />
        public int CompareTo(DcelPoint other)
        {
            if (x == other.x && y == other.y)
            {
                return 0;
            }

            return (x < other.x || (x == other.x && y < other.y)) ? -1 : 1;
        }

        /// <inheritdoc />
        public override int GetHashCode()
        {
            int hash = 19;
            hash = hash * 31 + x.GetHashCode();
            hash = hash * 31 + y.GetHashCode();

            return hash;
        }
    }

}
