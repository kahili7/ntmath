namespace nt.math
{
    using System;
    using System.Globalization;

    public sealed class DefaultMatrixFormatProvider : MatrixFormatProviderBase
    {
        public DefaultMatrixFormatProvider(IFormatProvider innerProvider) : base(innerProvider)
        {
            FormatMatrixStart = String.Empty;
            FormatMatrixEnd = String.Empty;
            FormatRowStart = String.Empty;
            FormatRowEnd = String.Empty;
            FormatColStart = String.Empty;
            FormatColEnd = String.Empty;
            FormatRowDelimiter = " \n";
            FormatColDelimiter = " ";

            ParseMatrixStart = String.Empty;
            ParseMatrixEnd = String.Empty;
            ParseRowStart = String.Empty;
            ParseRowEnd = String.Empty;
            ParseColStart = String.Empty;
            ParseColEnd = String.Empty;
            ParseRowDelimiter = "\n";
            ParseColDelimiter = " ";
        }

        public static DefaultMatrixFormatProvider CurrentCulture
        {
            get { return currentCulture; }
        }

        public static DefaultMatrixFormatProvider InvariantCulture
        {
            get { return invariantCulture; }
        }

        private static readonly DefaultMatrixFormatProvider currentCulture = new DefaultMatrixFormatProvider(CultureInfo.CurrentCulture);
        private static readonly DefaultMatrixFormatProvider invariantCulture = new DefaultMatrixFormatProvider(CultureInfo.InvariantCulture);
    }
}
