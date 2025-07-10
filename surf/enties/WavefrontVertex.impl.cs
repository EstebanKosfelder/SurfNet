namespace SurfNet
{
    using static DebugLog;
    using static Mathex;

    public partial class WavefrontVertex
    {


        //std.ostream & operator <<(std.ostream& os, const WavefrontVertex.LineIntersectionType t)
        //{
        //    switch (t)
        //    {
        //        case WavefrontVertex.LineIntersectionType.ONE: return os << "ONE";
        //        case WavefrontVertex.LineIntersectionType.ALL: return os << "ALL";
        //        case WavefrontVertex.LineIntersectionType.NONE: return os << "NONE";
        //    }
        //    CANNOTHAPPEN_MSG << "Fell through switch which should cover all cases.";
        //    assert(false);
        //    abort();
        //}

        //std.ostream &
        //operator <<(std.ostream& os, const InfiniteSpeedType &a)
        //{
        //    switch (a)
        //    {
        //        case InfiniteSpeedType.NONE:
        //            os << "";
        //            break;
        //        case InfiniteSpeedType.OPPOSING:
        //            os << "^";
        //            break;
        //        case InfiniteSpeedType.WEIGHTED:
        //            os << "~";
        //            break;
        //    }
        //    return os;
        //}
        public bool Virtual { get; internal set; }
    }
}