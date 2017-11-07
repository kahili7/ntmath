namespace nt.core
{
    using System;
    using System.Threading;

    internal class Lazy<T>
    {
        private readonly Func<T> valueFactory;
        private readonly object lockObj = new Object();

        private bool isValueCreated;
        private T value;

        public Lazy(Func<T> valueFactory) : this(valueFactory, true) { }

        public Lazy(Func<T> valueFactory, bool isThreadSafe)
        {
            if (valueFactory == null)
                throw new ArgumentNullException("valueFactory");

            if (isThreadSafe == false)
                throw new ArgumentException("This implementation only supports thread-safe instances.", "isThreadSafe");

            this.valueFactory = valueFactory;
        }

        public T Value
        {
            get
            {
                if (!isValueCreated)
                {
                    lock (lockObj)
                    {
                        if (!isValueCreated)
                        {
                            value = valueFactory();

                            Thread.MemoryBarrier();
                            isValueCreated = true;
                        }
                    }
                }

                Thread.MemoryBarrier();
                return value;
            }
        }

        public bool IsValueCreated
        {
            get
            {
                lock (lockObj)
                {
                    return isValueCreated;
                }
            }
        }

        public override string ToString()
        {
            return Value.ToString();
        }
    }
}