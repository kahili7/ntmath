namespace nt.math
{
    using System;

    public static partial class Matrix
    {
        #region ToString
        public static string ToString(this double[,] matrix)
        {
            return ToString(matrix, DefaultMatrixFormatProvider.CurrentCulture);
        }

        public static string ToString(this double[,] matrix, bool multiline, IMatrixFormatProvider provider)
        {
            if (multiline)
            {
                return ToString(matrix, Environment.NewLine, provider);
            }
            else
            {
                return ToString(matrix, null, provider);
            }
        }

        public static string ToString(this double[,] matrix, IMatrixFormatProvider provider)
        {
            return ToString(matrix, null, provider);
        }

        public static string ToString(this double[,] matrix, string format, IMatrixFormatProvider provider)
        {
            return MatrixFormatter.Format(format, matrix, provider);
        }

        public static string ToString(this double[,] matrix, string format)
        {
            return MatrixFormatter.Format(format, matrix, DefaultMatrixFormatProvider.CurrentCulture);
        }

        public static string ToString(this double[][] matrix)
        {
            return ToString(matrix, DefaultMatrixFormatProvider.CurrentCulture);
        }

        public static string ToString(this double[][] matrix, IMatrixFormatProvider provider)
        {
            return ToString(matrix, null, provider);
        }

        public static string ToString(this double[][] matrix, string format, IMatrixFormatProvider provider)
        {
            return MatrixFormatter.Format(format, matrix, provider);
        }

        public static string ToString(this double[][] matrix, string format)
        {
            return MatrixFormatter.Format(format, matrix, DefaultMatrixFormatProvider.CurrentCulture);
        }

        public static string ToString(this double[] array)
        {
            return ToString(array, DefaultArrayFormatProvider.CurrentCulture);
        }

        public static string ToString(this double[] array, IMatrixFormatProvider provider)
        {
            return ToString(array, null, provider);
        }

        public static string ToString(this double[] matrix, string format, IMatrixFormatProvider provider)
        {
            return MatrixFormatter.Format(format, matrix, provider);
        }

        public static string ToString(this double[] array, string format)
        {
            return MatrixFormatter.Format(format, array, DefaultArrayFormatProvider.CurrentCulture);
        }
        #endregion

        #region Parse
        public static double[,] Parse(string str)
        {
            return MatrixFormatter.ParseMultidimensional(str, DefaultMatrixFormatProvider.CurrentCulture);
        }

        public static double[,] Parse(string str, IMatrixFormatProvider provider)
        {
            return MatrixFormatter.ParseMultidimensional(str, provider);
        }

        public static double[][] ParseJagged(string s, IMatrixFormatProvider provider)
        {
            return MatrixFormatter.ParseJagged(s, provider);
        }

        public static bool TryParse(string s, IMatrixFormatProvider provider, out double[,] matrix)
        {
            try
            {
                matrix = Parse(s, provider);
            }
            catch (FormatException)
            {
                matrix = null;
            }
            catch (ArgumentNullException)
            {
                matrix = null;
            }

            return matrix != null;
        }

        public static bool TryParse(string s, IMatrixFormatProvider provider, out double[][] matrix)
        {
            try
            {
                matrix = ParseJagged(s, provider);
            }
            catch (FormatException)
            {
                matrix = null;
            }
            catch (ArgumentNullException)
            {
                matrix = null;
            }

            return matrix != null;
        }
        #endregion
    }
}
