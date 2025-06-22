using System.Runtime.Serialization;

namespace SurfNet
{
    public class GFLException : Exception
    {
        public GFLException()
        {
        }

        public GFLException(string? message) : base(message)
        {
        }

        public GFLException(string? message, Exception? innerException) : base(message, innerException)
        {
        }

        protected GFLException(SerializationInfo info, StreamingContext context) : base(info, context)
        {
        }
    }

    public class HasNotRelevantEdge : GFLException
    {
        public HasNotRelevantEdge()
        {
        }

        public HasNotRelevantEdge(string? message) : base(message)
        {
        }

        public HasNotRelevantEdge(string? message, Exception? innerException) : base(message, innerException)
        {
        }

        protected HasNotRelevantEdge(SerializationInfo info, StreamingContext context) : base(info, context)
        {
        }
    }
}
