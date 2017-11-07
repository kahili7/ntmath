namespace nt.math.comparers
{
    using System;
    using System.Collections.Generic;

    public class CustomComparer<T> : IComparer<T>, IEqualityComparer<T>
    {
        private Func<T, T, int> comparer;

        public CustomComparer(Func<T, T, int> comparer)
        {
            this.comparer = comparer;
        }

        public int Compare(T x, T y)
        {
            return comparer(x, y);
        }

        public bool Equals(T x, T y)
        {
            return comparer(x, y) == 0;
        }

        public int GetHashCode(T obj)
        {
            return obj.GetHashCode();
        }
    }
}
