using System.Collections;

namespace SurfNet
{
    public class Array3Nulleable<T> : IEnumerable<T> where T : class
    {
        public Array3Nulleable()
        {
            array = new T?[3];
        }

        public Array3Nulleable(params T?[] values)
        {
#if !NO_PRECONDITION
            if (values.Length != 3) throw new ArgumentOutOfRangeException(nameof(values));
#endif
            this.array = values;
        }

        protected T?[] array;

        public T? this[int index]
        {
            get
            {
#if !NO_PRECONDITION
                if (index < 0 && index > 2) throw new ArgumentOutOfRangeException(nameof(index));
#endif
                return array[index];
            }
            set
            {
#if !NO_PRECONDITION
                if (index < 0 && index > 2) throw new ArgumentOutOfRangeException(nameof(index));
#endif
                array[index] = value;
            }
        }

        public IEnumerator<T> GetEnumerator()
        {
            return ((IEnumerable<T?>)array).GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return array.GetEnumerator();
        }
    }

    public class Array3<T> : IEnumerable<T> where T : class
    {
        public Array3()
        {
            array = new T[3];
        }

        public Array3(params T[] values)
        {
#if !NO_PRECONDITION
            if (values.Length != 3) throw new ArgumentOutOfRangeException(nameof(values));
#endif
            this.array = values;
        }

        protected T[] array;

        public T this[int index]
        {
            get
            {
#if !NO_PRECONDITION
                if (index < 0 && index > 2) throw new ArgumentOutOfRangeException(nameof(index));
#endif
                return array[index];
            }
            set
            {
#if !NO_PRECONDITION
                if (index < 0 && index > 2) throw new ArgumentOutOfRangeException(nameof(index));
#endif
                array[index] = value;
            }
        }

        public IEnumerator<T> GetEnumerator()
        {
            return ((IEnumerable<T>)array).GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return array.GetEnumerator();
        }
    }

    public class WavefrontVertex3 : Array3<WavefrontVertex>
    {
        public override string ToString()
        {
            return $"v:[{string.Join(", ", array.Select(e => (e != null ? e.Id.ToString() : " *").PadLeft(3)))}]";
        }
        public WavefrontVertex3()
        {
        }

        public WavefrontVertex3(params WavefrontVertex[] values) : base(values)
        {
        }
    }

    public class WavefrontEdge3 : Array3Nulleable<WavefrontEdge>
    {
        public override string ToString()
        {
            return $"e:[{string.Join(", ", array.Select(e => (e != null ? e.Id.ToString() : " *").PadLeft(3)))}]";
        }

        public WavefrontEdge3()
        {
        }

        public WavefrontEdge3(params WavefrontEdge[] values) : base(values)
        {
        }
    }

    public class KineticTriangle3 : Array3Nulleable<KineticTriangle>
    {
        public override string ToString()
        {
            return $"n:[{string.Join(", ", array.Select(e => (e != null ? e.Id.ToString() : " *").PadLeft(3)))}]";
        }
        public KineticTriangle3()
        {
        }

        public KineticTriangle3(params KineticTriangle[] values) : base(values)
        {
        }
    }

    public class HalfEdge3 : Array3<DcelHalfEdge>
    {
        public HalfEdge3(params DcelHalfEdge[] values) : base(values)
        {
        }
    }
}