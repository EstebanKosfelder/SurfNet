
using System.Diagnostics;


namespace SurfNet
{
    public static class Extensions
    {


        public static void Precondition(bool value, string? msg = null)
        {
            if (!value)
                throw new Exception(msg);
        }

        public static void postcondition(bool value, string? msg = null)
        {
            if (!value)
                throw new Exception(msg);
        }

        public static void assertion(bool value, string? msg = null)
        {
            if (!value)
                throw new Exception(msg);
        }

        public static void assert(bool value, string msg)
        {
            if (!value)
                throw new Exception(msg);
        }


        public static void push_back<T>(this List<T> list, T item) => list.Add(item);

        public static T Last<T>(this List<T> list) => list[list.Count - 1];

        public static T First<T>(this List<T> list) => list[0];





    }
    public static class DebuggerInfo
    {

        public static int CGAL_STRAIGHT_SKELETON_ENABLE_TRACE = 4;
        public static bool CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE = true;
        public static int CGAL_POLYGON_OFFSET_ENABLE_TRACE = 4;



        public static void ExceptionHandler(ArithmeticException ex)
        {
            Console.Error.WriteLine(ex);
            Debug.WriteLine(ex);
        }



        public static void CGAL_STSKEL_VALIDITY_TRACE(params string[] m) => CGAL_STSKEL_TRACE(string.Join("\n", m));

        public static void CGAL_STSKEL_VALIDITY_TRACE_IF(bool value, params string[] m)
        {
            if (value) CGAL_STSKEL_VALIDITY_TRACE(m);
        }

        public static void CGAL_STSKEL_TRAITS_TRACE(params string[] m) => CGAL_STSKEL_TRACE(string.Join("\n", m));
        public static void CGAL_STSKEL_BUILDER_TRACE(int l, params string[] m) => CGAL_STSKEL_BUILDER_TRACE(l, string.Join("\n", m));
        public static void CGAL_STSKEL_BUILDER_TRACE(int l, string m) { if (l <= CGAL_STRAIGHT_SKELETON_ENABLE_TRACE) CGAL_STSKEL_TRACE(m); }


        public static void CGAL_STSKEL_BUILDER_TRACE_IF(bool c, int l, string m)
        {
            if ((c) && l <= CGAL_STRAIGHT_SKELETON_ENABLE_TRACE) CGAL_STSKEL_TRACE(m);
        }

        public static void CGAL_STSKEL_TRACE(string m)
        {
            Console.WriteLine(m);
            Debug.WriteLine(m);
        }
        public static void CGAL_precondition(bool value, string? msg = null)
        {
            Extensions.Precondition(value, msg);
        }
        public static void CGAL_postcondition(bool value, string? msg = null) => Extensions.postcondition(value, msg);
        public static void CGAL_kernel_precondition(bool value, string? msg = null) => Extensions.Precondition(value, msg);

        public static void CGAL_assertion(bool value, string? msg = null) => Extensions.assertion(value, msg);

        public static void CGAL_kernel_assertion(bool value) => Extensions.assertion(value);



    }
}
