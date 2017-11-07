namespace nt.math
{
    using System;

    public static partial class Matrix
    {
        #region ToMatrix
        public static T[,] ToMatrix<T>(this T[][] array)
        {
            return ToMatrix(array, false);
        }

        public static T[,] ToMatrix<T>(this T[][] array, bool transpose)
        {
            int rows = array.Length;
            if (rows == 0) return new T[0, rows];
            int cols = array[0].Length;
            T[,] m;

            if (transpose)
            {
                m = new T[cols, rows];

                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        m[j, i] = array[i][j];
            }
            else
            {
                m = new T[rows, cols];

                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        m[i, j] = array[i][j];
            }

            return m;
        }

        public static T[,] ToMatrix<T>(this T[] array)
        {
            T[,] m = new T[1, array.Length];

            for (int i = 0; i < array.Length; i++)
                m[0, i] = array[i];

            return m;
        }

        public static T[][] ToArray<T>(this T[] array, bool asColumnVector = true)
        {
            if (asColumnVector)
            {
                T[][] m = new T[array.Length][];

                for (int i = 0; i < array.Length; i++)
                    m[i] = new[] { array[i] };
                return m;
            }
            else
            {
                return new T[][] { array };
            }
        }

        public static T[,] ToMatrix<T>(this T[] array, bool asColumnVector = false)
        {
            if (asColumnVector)
            {
                T[,] m = new T[1, array.Length];

                for (int i = 0; i < array.Length; i++)
                    m[0, i] = array[i];
                return m;
            }
            else
            {
                T[,] m = new T[array.Length, 1];

                for (int i = 0; i < array.Length; i++)
                    m[i, 0] = array[i];
                return m;
            }
        }

        public static T[][] ToArray<T>(this T[,] matrix)
        {
            return ToArray(matrix, false);
        }

        public static T[][] ToArray<T>(this T[,] matrix, bool transpose)
        {
            T[][] array;

            if (transpose)
            {
                int cols = matrix.GetLength(1);

                array = new T[cols][];

                for (int i = 0; i < cols; i++)
                    array[i] = matrix.GetColumn(i);
            }
            else
            {
                int rows = matrix.GetLength(0);

                array = new T[rows][];

                for (int i = 0; i < rows; i++)
                    array[i] = matrix.GetRow(i);
            }

            return array;
        }
        #endregion

    }
}
