namespace nt.math
{
    using System;

    public abstract class MatrixFormatProviderBase : IMatrixFormatProvider
    {
        #region Formatting specification
        public string FormatMatrixStart { get; protected set; }

        public string FormatMatrixEnd { get; protected set; }

        public string FormatRowStart { get; protected set; }

        public string FormatRowEnd { get; protected set; }

        public string FormatColStart { get; protected set; }

        public string FormatColEnd { get; protected set; }

        public string FormatRowDelimiter { get; protected set; }

        public string FormatColDelimiter { get; protected set; }
        #endregion

        #region Parsing specification
        public string ParseMatrixStart { get; protected set; }

        public string ParseMatrixEnd { get; protected set; }

        public string ParseRowStart { get; protected set; }

        public string ParseRowEnd { get; protected set; }

        public string ParseColStart { get; protected set; }

        public string ParseColEnd { get; protected set; }

        public string ParseRowDelimiter { get; protected set; }

        public string ParseColDelimiter { get; protected set; }
        #endregion

        public IFormatProvider InnerProvider { get; protected set; }

        protected MatrixFormatProviderBase(IFormatProvider innerProvider)
        {
            this.InnerProvider = innerProvider;
        }

        public object GetFormat(Type formatType)
        {
            if (formatType == typeof(ICustomFormatter))
            {
                return new MatrixFormatter();
            }

            return null;
        }
    }
}
