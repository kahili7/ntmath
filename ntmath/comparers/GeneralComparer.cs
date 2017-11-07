namespace nt.math.comparers
{
    using System;
    using System.Collections.Generic;

    public enum ComparerDirection
    {
        Ascending = +1,
        Descending = -1
    };

    public class GeneralComparer : IComparer<double>, IComparer<int>
    {
        Func<double, double> map;
        private int direction = 1;

        public ComparerDirection Direction
        {
            get { return (ComparerDirection)direction; }
            set { direction = (int)value; }
        }

        public GeneralComparer(ComparerDirection direction) : this(direction, false) { }

        public GeneralComparer(ComparerDirection direction, bool useAbsoluteValues)
        {
            if (useAbsoluteValues) this.map = Math.Abs;
            else this.map = (a) => a;

            this.direction = (int)direction;
        }

        public GeneralComparer(ComparerDirection direction, Func<double, double> map)
        {
            this.map = map;
            this.direction = (int)direction;
        }

        public int Compare(double x, double y)
        {
            return direction * (map(x).CompareTo(map(y)));
        }

        public int Compare(int x, int y)
        {
            return direction * (map(x).CompareTo(map(y)));
        }
    }
}
