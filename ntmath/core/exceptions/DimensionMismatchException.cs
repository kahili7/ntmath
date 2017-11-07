using System;
using System.Runtime.Serialization;

namespace nt.core
{
    [Serializable]
    public class DimensionMismatchException : ArgumentException
    {
        public DimensionMismatchException() { }

        public DimensionMismatchException(string paramName) :
            base(paramName, "Array dimensions must match.") { }

        public DimensionMismatchException(string paramName, string message) :
            base(message, paramName) { }

        public DimensionMismatchException(string message, Exception innerException) :
            base(message, innerException) { }

        protected DimensionMismatchException(SerializationInfo info, StreamingContext context) :
            base(info, context) { }
    }
}
