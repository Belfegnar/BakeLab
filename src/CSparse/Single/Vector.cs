// -----------------------------------------------------------------------
// <copyright file="Vector.cs">
// Copyright (c) 2012-2016, Christian Woltering
// </copyright>
// -----------------------------------------------------------------------

namespace CSparse.Single
{
    using System;

    /// <summary>
    /// Vector helper methods.
    /// </summary>
    public static class Vector
    {
        /// <summary>
        /// Copy one vector to another.
        /// </summary>
        /// <param name="src">The source array.</param>
        /// <param name="dst">The destination array.</param>
        public static void Copy(float[] src, float[] dst)
        {
            Buffer.BlockCopy(src, 0, dst, 0, src.Length * Constants.SizeOfSingle);
        }

        /// <summary>
        /// Copy one vector to another.
        /// </summary>
        /// <param name="src">The source array.</param>
        /// <param name="dst">The destination array.</param>
        /// <param name="n">Number of values to copy.</param>
        public static void Copy(float[] src, float[] dst, int n)
        {
            Buffer.BlockCopy(src, 0, dst, 0, n * Constants.SizeOfSingle);
        }

        /// <summary>
        /// Create a new vector.
        /// </summary>
        public static float[] Create(int length, float value)
        {
            float[] result = new float[length];

            for (int i = 0; i < length; i++)
            {
                result[i] = value;
            }

            return result;
        }

        /// <summary>
        /// Clone the given vector.
        /// </summary>
        public static float[] Clone(float[] src)
        {
            float[] result = new float[src.Length];

            Buffer.BlockCopy(src, 0, result, 0, src.Length * Constants.SizeOfSingle);

            return result;
        }

        /// <summary>
        /// Set vector values to zero.
        /// </summary>
        public static void Clear(float[] x)
        {
            Array.Clear(x, 0, x.Length);
        }

        /// <summary>
        /// Computes the dot product of two vectors.
        /// </summary>
        public static float DotProduct(float[] x, float[] y)
        {
            int length = x.Length;

            float result = 0.0f;

            for (int i = 0; i < length; i++)
            {
                result += x[i] * y[i];
            }

            return result;
        }

        /// <summary>
        /// Computes the pointwise product of two vectors.
        /// </summary>
        public static void PointwiseMultiply(float[] x, float[] y, float[] z)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                z[i] = x[i] * y[i];
            }
        }

        /// <summary>
        /// Computes the norm of a vector, sqrt( x' * x ).
        /// </summary>
        public static float Norm(float[] x)
        {
            int length = x.Length;

            float result = 0.0f;

            for (int i = 0; i < length; ++i)
            {
                result += x[i] * x[i];
            }

			return (float)Math.Sqrt(result);
        }

        /// <summary>
        /// Computes the norm of a vector avoiding overflow, sqrt( x' * x ).
        /// </summary>
        public static float NormRobust(float[] x)
        {
            int length = x.Length;

            float scale = 0.0f, ssq = 1.0f;

            for (int i = 0; i < length; ++i)
            {
                if (x[i] != 0.0f)
                {
                    float absxi = Math.Abs(x[i]);
                    if (scale < absxi)
                    {
                        ssq = 1.0f + ssq * (scale / absxi) * (scale / absxi);
                        scale = absxi;
                    }
                    else
                    {
                        ssq += (absxi / scale) * (absxi / scale);
                    }
                }
            }

			return scale * (float)Math.Sqrt(ssq);
        }

        /// <summary>
        /// Scales a vector by a given factor, x = alpha * x.
        /// </summary>
        public static void Scale(float alpha, float[] x)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                x[i] *= alpha;
            }
        }

        /// <summary>
        /// Add a scaled vector to another vector, y = y + alpha * x.
        /// </summary>
        public static void Axpy(float alpha, float[] x, float[] y)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                y[i] += alpha * x[i];
            }
        }

        /// <summary>
        /// Add two scaled vectors, z = alpha * x + beta * y.
        /// </summary>
        public static void Add(float alpha, float[] x, float beta, float[] y, float[] z)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                z[i] = alpha * x[i] + beta * y[i];
            }
        }
    }
}
