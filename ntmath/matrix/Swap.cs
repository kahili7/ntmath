namespace nt.math
{
    using System;
    using System.Collections.Generic;
    using nt.core;

    public static partial class Matrix
    {
        public static void Swap<T>(this T[] source, T[] destination)
        {
            if(source.GetLength(0) != destination.GetLength(0))
                throw new ArgumentException("Argument out of range.");

            for(int i = 0; i < source.GetLength(0); i++)
            {
                T tmp = source[i];

                source[i] = destination[i];
                destination[i] = tmp;
            }
        }

        public static void Swap<T>(this T[,] source, T[,] destination)
        {
            if (source.GetLength(0) != destination.GetLength(0))
                throw new ArgumentException("Argument out of range.");

            if (source.GetLength(1) != destination.GetLength(1))
                throw new ArgumentException("Argument out of range.");

            for (int i = 0; i < source.GetLength(0); i++)
            {
                for (int j = 0; j < source.GetLength(1); j++)
                {
                    T tmp = source[i, j];

                    source[i, j] = destination[i, j];
                    destination[i, j] = tmp;
                }
            }
        }
    }
}
