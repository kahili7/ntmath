namespace nt.math.decompositions
{
    public interface ISolverMatrixDecomposition<T>  where T : struct  
    {
        T[,] Solve(T[,] value);

        T[] Solve(T[] value);

        T[,] Inverse();
    }
}
