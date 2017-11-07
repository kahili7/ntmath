namespace nt.math
{
    using System;
    using System.Collections.Generic;

    public static class Combinatorics
    {
        public static int[][] TruthTable(int length)
        {
            return TruthTable(2, length);
        }

        public static int[][] TruthTable(int symbols, int length)
        {
            int[] sym = new int[length];

            for (int i = 0; i < sym.Length; i++)
                sym[i] = symbols;

            return TruthTable(sym);
        }

        public static int[][] TruthTable(int[] symbols)
        {
            int size = 1;

            for (int i = 0; i < symbols.Length; i++)
                size *= symbols[i];

            int[][] sequences = new int[size][];
            sequences[0] = new int[symbols.Length];

            for (int i = 1; i < sequences.Length; i++)
            {
                var row = sequences[i] = (int[])sequences[i - 1].Clone();

                for (int j = symbols.Length - 1; j >= 0; j--)
                {
                    if (row[j] < symbols[j] - 1)
                    {
                        row[j]++;
                        break;
                    }

                    row[j] = 0;
                }
            }

            return sequences;
        }

        public static IEnumerable<int[]> Sequences(int symbols, int length, bool inPlace = false)
        {
            int[] sym = new int[length];

            for (int i = 0; i < sym.Length; i++)
                sym[i] = symbols;

            return Sequences(sym, inPlace);
        }

        public static IEnumerable<int[]> Sequences(int[] symbols, bool inPlace = false)
        {
            var current = new int[symbols.Length];

            while (true)
            {
                yield return inPlace ? current : (int[])current;

                bool match = true;

                for (int i = 0; i < current.Length && match; i++)
                    if (current[i] != symbols[i] - 1)
                        match = false;

                if (match)
                    break;

                for (int j = symbols.Length - 1; j >= 0; j--)
                {
                    if (current[j] < symbols[j] - 1)
                    {
                        current[j]++;
                        break;
                    }

                    current[j] = 0;
                }
            }
        }

        public static IEnumerable<T[]> Combinations<T>(T[] values, int k, bool inPlace = false)
        {
            int n = values.Length;
            int t = k;
            int[] c = new int[t + 3];
            T[] current = new T[t];
            int j, x;

            for (j = 1; j <= t; j++)
                c[j] = j - 1;

            c[t + 1] = n;
            c[t + 2] = 0;
            j = t;

            do
            {
                for (int i = 0; i < current.Length; i++)
                    current[i] = values[c[i + 1]];

                yield return (inPlace ? current : (T[])current.Clone());

                if (j > 0)
                {
                    x = j;
                }
                else
                {
                    if (c[1] + 1 < c[2])
                    {
                        c[1]++;
                        continue;
                    }
                    else
                    {
                        j = 2;
                    }
                }

                while (true)
                {
                    c[j - 1] = j - 2;
                    x = c[j] + 1;

                    if (x == c[j + 1]) j++;
                    else break;
                }

                c[j] = x;
                j--;
            } while (j < t);
        }

        public static IEnumerable<T[]> Permutations<T>(T[] values, bool inPlace = false)
        {
            T[] current = new T[values.Length];

            yield return (inPlace ? values : (T[])values.Clone());

            int[] idx = Matrix.Indices(0, values.Length);
            int j, l;

            while (true)
            {
                for (j = values.Length - 2; j >= 0; j--)
                    if (idx[j + 1] > idx[j]) break;

                if (j == -1) yield break;

                for (l = values.Length - 1; l > j; l--)
                    if (idx[l] > idx[j]) break;

                int temp = idx[j];

                idx[j] = idx[l];
                idx[l] = temp;

                for (int i = j + 1; i < idx.Length; i++)
                {
                    if (i > idx.Length - i + j) break;

                    temp = idx[i];
                    idx[i] = idx[idx.Length - i + j];
                    idx[idx.Length - i + j] = temp;
                }

                for (int i = 0; i < values.Length; i++)
                    current[i] = values[idx[i]];

                yield return (inPlace ? current : (T[])current.Clone());
            }
        }
    }
}
