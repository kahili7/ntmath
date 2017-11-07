namespace nt.math
{
    using System;
    using System.Globalization;

    public sealed class DefaultArrayFormatProvider : MatrixFormatProviderBase
    {
        public DefaultArrayFormatProvider(IFormatProvider innerProvider) : base(innerProvider)
        {
            FormatMatrixStart = String.Empty;
            FormatMatrixEnd = String.Empty;
            FormatRowStart = String.Empty;
            FormatRowEnd = String.Empty;
            FormatColStart = String.Empty;
            FormatColEnd = String.Empty;
            FormatRowDelimiter = " ";
            FormatColDelimiter = String.Empty;

            ParseMatrixStart = String.Empty;
            ParseMatrixEnd = String.Empty;
            ParseRowStart = String.Empty;
            ParseRowEnd = String.Empty;
            ParseColStart = String.Empty;
            ParseColEnd = String.Empty;
            ParseRowDelimiter = " ";
            ParseColDelimiter = String.Empty;
        }

        public static DefaultArrayFormatProvider CurrentCulture
        {
            get { return currentCulture; }
        }

        public static DefaultArrayFormatProvider InvariantCulture
        {
            get { return invariantCulture; }
        }

        private static readonly DefaultArrayFormatProvider currentCulture = new DefaultArrayFormatProvider(CultureInfo.CurrentCulture);
        private static readonly DefaultArrayFormatProvider invariantCulture = new DefaultArrayFormatProvider(CultureInfo.InvariantCulture);
    }
}
