using System.Globalization;

namespace SurfNet
{
    public static class Extension
    {
        public static string Debug(this double num, int dec = 3) =>  num.ToString($"F{dec}",CultureInfo.InvariantCulture);
        public static T back<T>(this List<T> list) => list[list.Count - 1];
        public static int size<T>(this List<T> list) => list.Count;

        public static int ToInt(this bool value) => value ? 1 : 0;

    }
}
