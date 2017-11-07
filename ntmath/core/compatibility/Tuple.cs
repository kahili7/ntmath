namespace nt.core
{
    using System.Collections.Generic;

    public class Tuple<T1, T2>
    {

        /// 
        public T1 Item1 { get; set; }

        /// 
        public T2 Item2 { get; set; }

        /// 
        public Tuple(T1 item1, T2 item2)
        {
            Item1 = item1;
            Item2 = item2;
        }

        public override bool Equals(object obj)
        {
            var other = obj as Tuple<T1, T2>;

            if (other == null) return false;

            return EqualityComparer<T1>.Default.Equals(Item1, other.Item1)
                && EqualityComparer<T2>.Default.Equals(Item2, other.Item2);
        }

        public override int GetHashCode()
        {
            return (Item1 == null ? 1 : Item1.GetHashCode()) ^ (Item2 == null ? 2 : Item2.GetHashCode());
        }
    }

    [System.Diagnostics.CodeAnalysis.SuppressMessage("Microsoft.Design", "CA1005:AvoidExcessiveParametersOnGenericTypes")]
    public class Tuple<T1, T2, T3>
    {
        public T1 Item1 { get; set; }
 
        public T2 Item2 { get; set; }

        public T3 Item3 { get; set; }

        public Tuple(T1 item1, T2 item2, T3 item3)
        {
            Item1 = item1;
            Item2 = item2;
            Item3 = item3;
        }
 
        public override bool Equals(object obj)
        {
            var other = obj as Tuple<T1, T2, T3>;

            if (other == null) return false;

            return EqualityComparer<T1>.Default.Equals(Item1, other.Item1)
                && EqualityComparer<T2>.Default.Equals(Item2, other.Item2)
                && EqualityComparer<T3>.Default.Equals(Item3, other.Item3);
        }

        public override int GetHashCode()
        {
            return (Item1 == null ? 1 : Item1.GetHashCode())
                 ^ (Item2 == null ? 2 : Item2.GetHashCode())
                 ^ (Item3 == null ? 3 : Item3.GetHashCode());
        }
    }

    public static class Tuple
    {
        public static Tuple<T1, T2> Create<T1, T2>(T1 item1, T2 item2)
        {
            var tuple = new Tuple<T1, T2>(item1, item2);

            return tuple;
        }

        public static Tuple<T1, T2, T3> Create<T1, T2, T3>(T1 item1, T2 item2, T3 item3)
        {
            var tuple = new Tuple<T1, T2, T3>(item1, item2, item3);

            return tuple;
        }
    }
}