namespace nt.math
{
    using System;

    public interface IMatrixFormatProvider : IFormatProvider  
    {
        #region Formatting specification
        string FormatMatrixStart { get; }

        string FormatMatrixEnd { get; }

        string FormatRowStart { get; }

        string FormatRowEnd { get; }

        string FormatColStart { get; }

        string FormatColEnd { get; }

        string FormatRowDelimiter { get; }

        string FormatColDelimiter { get; }
        #endregion

        #region Parsing specification
        string ParseMatrixStart { get; }

        string ParseMatrixEnd { get; }

        string ParseRowStart { get; }

        string ParseRowEnd { get; }

        string ParseColStart { get; }

        string ParseColEnd { get; }

        string ParseRowDelimiter { get; }

        string ParseColDelimiter { get; }
        #endregion

        IFormatProvider InnerProvider { get; }
    }
}
